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

#include <vector>
#include <string>
#include <stdexcept>

#include <Wt/WWebWidget>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;

const int ReferenceLineInfo::sm_xmlSerializationVersion = 1;


namespace
{
  template<class T> struct index_compare_assend
  {
    index_compare_assend(const T arr, const T arr2) : m_arr(arr), m_arr2( arr2 ) {} //pass the actual values you want sorted into here
    bool operator()(const size_t a, const size_t b) const
    {
      if( m_arr[a] == m_arr[b] )
        return m_arr2[a] < m_arr2[b];
      return m_arr[a] < m_arr[b];
    }
    const T m_arr;
    const T m_arr2;
  };//struct index_compare
  
  std::string jsQuote( const std::string &str )
  {
    return Wt::WWebWidget::jsStringLiteral(str,'\"');
  }
  
}//namespace


//I need practice/help memorizing background lines, so I typed the following in
//  from "Practical gamma-ray spectrometry" pg 361
const OtherRefLine BackgroundReactionLines[28] =
{
  OtherRefLine( 53.44f,   0.1034f, "Ge(n,g)", OtherRefLineType::BackgroundReaction, "Ge72(n,g), Ge74(n,2n)" ),
  OtherRefLine( 68.75f,   0.001f,   "Ge(n,n)", OtherRefLineType::BackgroundReaction, "Ge73(n,n) broad antisymetric peak" ),
  OtherRefLine( 139.68f,  0.390f,  "Ge75m", OtherRefLineType::BackgroundReaction, "Ge74(n,g), Ge76(n,2n)" ),
  OtherRefLine( 159.7f,   0.1033f, "Ge(n,g)", OtherRefLineType::BackgroundReaction, "Ge76(n,g)" ),
  OtherRefLine( 174.95f,  0.0f,   "Ge(n,g)", OtherRefLineType::BackgroundReaction, "Ge70(n,g) activation" ),
  OtherRefLine( 198.39f,  0.0f,   "Ge71m", OtherRefLineType::BackgroundReaction, "Sum peak Ge70(n,g)" ),
  OtherRefLine( 278.26f,  0.0f,   "Cu64", OtherRefLineType::BackgroundReaction, "Cu63(n,g), Cu65(n,2n) prompt gamma" ),
  OtherRefLine( 336.24f,  0.459f,  "Cd115m,In115m", OtherRefLineType::BackgroundReaction, "Activation of Cd (descendant of Cd115)" ),
  OtherRefLine( 416.86f,  0.277f,  "In116m", OtherRefLineType::BackgroundReaction, "In115(n,g) activation of In metal seal" ),
  OtherRefLine( 527.90f,  0.275f,  "Cd115", OtherRefLineType::BackgroundReaction, "Cd114(n,g) activation" ),
  OtherRefLine( 558.46f,  0.0f,   "Cd114", OtherRefLineType::BackgroundReaction, "Cd113(n,g) prompt gamma" ),
  OtherRefLine( 569.7f,   0.9789f, "Pb207m", OtherRefLineType::BackgroundReaction, "Pb207(n,n)" ),
  OtherRefLine( 579.2f,   0.0f,   "Pb207", OtherRefLineType::BackgroundReaction, "Pb207(n,n) prompt gamma" ),
  OtherRefLine( 595.85f,  0.0f,   "Ge74", OtherRefLineType::BackgroundReaction, "Ge74(n,n) broad asymmetric peak" ),
  OtherRefLine( 669.62f,  0.0f,   "Cu63", OtherRefLineType::BackgroundReaction, "Cu63(n,n) prompt gamma" ),
  OtherRefLine( 689.6f,   0.0f,   "Ge72", OtherRefLineType::BackgroundReaction, "Ge72(n,n) broad asymetric peak" ),
  OtherRefLine( 803.06f,  0.0f,   "Pb206", OtherRefLineType::BackgroundReaction, "Pb206(n,n) prompt gamma" ),
  OtherRefLine( 843.76f,  0.718f,  "Mg27", OtherRefLineType::BackgroundReaction, "Mg26(n,g) or Al27(n,p) of encapsilation" ),
  OtherRefLine( 846.77f,  0.0f,   "Fe56", OtherRefLineType::BackgroundReaction, "Fe56(n,n)" ),
  OtherRefLine( 962.06f,  0.0f,   "Cu63", OtherRefLineType::BackgroundReaction, "Cu63(n,n) prompt gamma" ),
  OtherRefLine( 1014.44f, 0.280f,  "Mg27", OtherRefLineType::BackgroundReaction, "Mg26(n,g) or Al27(n,p) of encapsilation" ),
  OtherRefLine( 1063.66f, 0.885f,  "Pb207m", OtherRefLineType::BackgroundReaction, "Pb207(n,n)" ),
  OtherRefLine( 1097.3f,  0.562f,  "In116", OtherRefLineType::BackgroundReaction, "In115(n,g) activation of In metal seal" ),
  OtherRefLine( 1115.56f, 0.0f,   "Cu65", OtherRefLineType::BackgroundReaction, "Cu65(n,n)" ),
  OtherRefLine( 1173.23f, 0.9985f, "Co60", OtherRefLineType::BackgroundReaction, "Activation" ),
  OtherRefLine( 1293.54f, 0.844f,  "In116", OtherRefLineType::BackgroundReaction, "In115(n,g) activation of In metal seal" ),
  OtherRefLine( 1332.49f, 0.9998f, "Co60", OtherRefLineType::BackgroundReaction, "Activation" ),
  OtherRefLine( 2224.57f, 0.0f,   "H2", OtherRefLineType::BackgroundReaction, "H(n,g)" )
};//BackgroundReactionLines[]


//BackgroundLines looks to be taking up ~14 kb of executable size on Win7
const OtherRefLine BackgroundLines[89] =
{
  /*xrays below 46 keV not inserted*/
  OtherRefLine( 46.54f,   0.0425f,  "Pb210", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 53.23f,   0.01060f, "Pb214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 63.28f,   0.048f,   "Th234", OtherRefLineType::U238Series,         "" ),
  OtherRefLine( 72.81f,   0.277f,   "Pb xray", OtherRefLineType::BackgroundXRay,     "Flourescence and Tl208 decay" ),
  OtherRefLine( 74.82f,   0.277f,   "Bi xray", OtherRefLineType::BackgroundXRay,     "Pb212, Pb214 decay" ),
  OtherRefLine( 74.97f,   0.462f,   "Pb xray", OtherRefLineType::BackgroundXRay,     "Flourescence and Tl208 decay" ),
  OtherRefLine( 77.11f,   0.462f,   "Bi xray", OtherRefLineType::BackgroundXRay,     "Pb212, Pb214 decay" ),
  OtherRefLine( 79.29f,   0.461f,   "Po xray", OtherRefLineType::BackgroundXRay,     "Fluorescence and Bi212, Bi214 decay" ),
  OtherRefLine( 81.23f,   0.009f,   "Th231", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 84.94f,   0.107f,   "Pb xray", OtherRefLineType::BackgroundXRay,     "Flourescence and Tl208 decay" ),
  OtherRefLine( 87.3f,    0.0391f,  "Pb xray", OtherRefLineType::BackgroundXRay,     "Flourescence and Tl208 decay" ),
  OtherRefLine( 87.35f,   0.107f,   "Bi xray", OtherRefLineType::BackgroundXRay,     "Pb212, Pb214 decay" ),
  OtherRefLine( 89.78f,   0.0393f,  "Bi xray", OtherRefLineType::BackgroundXRay,     "Pb212, Pb214 decay" ),
  OtherRefLine( 89.96f,   0.281f,   "Th xray", OtherRefLineType::BackgroundXRay,     "U235 and Ac228 decay" ),
  OtherRefLine( 92.58f,   0.0558f,  "Th234", OtherRefLineType::U238Series,         " - doublet" ),
  OtherRefLine( 93.35f,   0.454f,   "Th xray", OtherRefLineType::BackgroundXRay,     "U235 and Ac228 decay" ),
  OtherRefLine( 105.6f,   0.107f,   "Th xray", OtherRefLineType::BackgroundXRay,     "U235 and Ac228 decay" ),
  OtherRefLine( 109.16f,  0.0154f,  "U235", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 112.81f,  0.0028f,  "Th234", OtherRefLineType::U238Series,         "" ),
  OtherRefLine( 122.32f,  0.01192f, "Ra223", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 129.06f,  0.0242f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 143.76f,  0.1096f,  "U235", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 163.33f,  0.0508f,  "U235", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 185.72f,  0.572f,   "U235", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 186.21f,  0.03555f, "Ra226", OtherRefLineType::U238Series,         "" ),
  OtherRefLine( 205.31f,  0.0501f,  "U235", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 209.26f,  0.0389f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 238.63f,  0.436f,   "Pb212", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 240.89f,  0.0412f,  "Ra224", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 242.0f,   0.07268f, "Pb214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 269.49f,  0.137f,   "Ra223", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 270.24f,  0.0346f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 277.37f,  0.0237f,  "Tl208", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 295.22f,  0.185f,   "Pb214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 299.98f,  0.0216f,  "Th227", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 300.07f,  0.0247f,  "Pa231", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 300.09f,  0.0318f,  "Pb212", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 328.0f,   0.0295f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 338.28f,  0.0279f,  "Ra223", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 338.32f,  0.1127f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 351.06f,  0.1291f,  "Bi211", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 351.93f,  0.3560f,  "Pb214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 409.46f,  0.0192f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 447.6f,   0.1044f,  "Be7", OtherRefLineType::OtherBackground,    "Cosmic" ),
  OtherRefLine( 462.0f,   0.00213f, "Pb214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 463.0f,   0.0440f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 510.7f,   0.0629f,  "Tl208", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 511.0f,   0.010f,   "", OtherRefLineType::OtherBackground,    "Annihilation radiation (beta+)" ),
  OtherRefLine( 570.82f,  0.00182f, "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 583.19f,  0.306f,   "Tl208", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 609.31f,  0.4549f,  "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 661.66f,  0.400f,   "Cs137", OtherRefLineType::OtherBackground,    "Fission" ),
  OtherRefLine( 726.86f,  0.0062f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 727.33f,  0.0674f,  "Bi212", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 755.31f,  0.010f,   "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 768.36f,  0.04891f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 794.95f,  0.0425f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 806.17f,  0.01262f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 832.01f,  0.0352f,  "Pb211", OtherRefLineType::U235Series,         "" ),
  OtherRefLine( 835.71f,  0.0161f,  "AC228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 839.04f,  0.00587f, "Pb214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 860.56f,  0.0448f,  "Tl208", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 911.20f,  0.258f,   "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 934.06f,  0.03096f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 964.77f,  0.0499f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 968.97f,  0.158f,   "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1001.03f, 0.01021f, "Pa234m", OtherRefLineType::U238Series,         "..." ),
  OtherRefLine( 1120.29f, 0.1491f,  "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1155.19f, 0.01635f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1238.11f, 0.05827f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1377.67f, 0.03967f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1407.98f, 0.02389f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1459.14f, 0.0083f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1460.82f, 0.1066f,  "K40", OtherRefLineType::K40Background,      "" ),
  OtherRefLine( 1588.2f,  0.0322f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1592.51f, 0.010f,   "Th232 D.E. 2614 keV", OtherRefLineType::OtherBackground,    "" ),
  OtherRefLine( 1620.74f, 0.0151f,  "Bi212", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1630.63f, 0.0151f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1686.09f, 0.00095f, "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1661.28f, 0.0112f,  "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1729.6f,  0.02843f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1764.49f, 0.1528f,  "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1838.36f, 0.00346f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1847.42f, 0.02023f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 2103.51f, 0.020f,   "Th232 S.E. 2614 keV", OtherRefLineType::OtherBackground,    "" ),
  OtherRefLine( 2118.55f, 0.00454f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 2204.21f, 0.04913f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 2447.86f, 0.01518f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 2614.51f, 0.3585f,  "Tl208", OtherRefLineType::Th232Series,        "Th232 series; Pb208(n,p)" ),
};//BackgroundLines



const char *to_str( const OtherRefLineType type )
{
  switch( type )
  {
    case OtherRefLineType::U238Series: return "U238Series";
    case OtherRefLineType::U235Series: return "U235Series";
    case OtherRefLineType::Th232Series: return "Th232Series";
    case OtherRefLineType::Ra226Series: return "Ra226Series";
    case OtherRefLineType::K40Background: return "K40Background";
    case OtherRefLineType::BackgroundXRay: return "BackgroundXRay";
    case OtherRefLineType::BackgroundReaction: return "BackgroundReaction";
    case OtherRefLineType::OtherBackground: return "Other";
  }//switch( type )

  return "InvalidOtherRefLineType";
}//to_str(...)


OtherRefLineType other_ref_line_type_from_str( const std::string &str )
{
  for( auto type = OtherRefLineType(0);
    type != OtherRefLineType::OtherBackground;
    type = OtherRefLineType( static_cast<int>(type) + 1) )
  {
    if( str == to_str( type ) )
      return type;
  }

  return OtherRefLineType::OtherBackground;
}


#if( DEV_REF_LINE_UPGRADE_20221212 )
ReferenceLineInfo::RefLineInput::RefLineInput()
  : m_input_txt(),
  m_age(),
  m_color(),
  m_lower_br_cutt_off( 0.0 ),
  m_promptLinesOnly( false ),
  m_showGammas( false ),
  m_showXrays( false ),
  m_showAlphas( false ),
  m_showBetas( false ),
  m_showCascades( false ),
  m_detector_name(),
  m_det_intrinsic_eff(),
  m_shielding_name(),
  m_shieldingThickness( 0.0 ),
  m_shielding_att()
{
}


ReferenceLineInfo::RefLine::RefLine()
  : m_energy( 0.0 ),
  m_normalized_intensity( 0.0 ),
  m_drf_factor( 1.0f ),
  m_shield_atten( 1.0f ),
  m_particle_sf_applied( 1.0f ),
  m_particlestr(),
  m_decaystr(),
  m_elementstr(),
  m_decay_intensity( 0.0 ),
  m_particle_type( ReferenceLineInfo::RefLine::Particle::Gamma ),
  m_parent_nuclide( nullptr ),
  m_transition( nullptr ),
  m_source_type( ReferenceLineInfo::RefLine::RefGammaType::Normal ),
  m_element( nullptr ),
  m_reaction( nullptr )
{
}
#endif //#if( DEV_REF_LINE_UPGRADE_20221212 )

ReferenceLineInfo::ReferenceLineInfo()
{
  reset();
}

void ReferenceLineInfo::reset()
{
  nuclide = NULL;
  element = NULL;
#if( DEV_REF_LINE_UPGRADE_20221212 )
  m_ref_lines.clear();
  m_input_warnings.clear();
  m_validity = InputValidity::Blank;
  m_has_coincidences = false;
  m_input = RefLineInput();
  m_source_type = ReferenceLineInfo::SourceType::None;
#endif
  energies.clear();
  intensities.clear();
  particlestrs.clear();
  decaystrs.clear();
  elementstrs.clear();
  particle_sf.clear();
  reactionGammas.clear();
  otherRefLines.clear();
  labelTxt.clear();
  reactionsTxt.clear();
  lineColor = Wt::WColor();
  showGammas = showXrays = showAlphas = showBetas = showCascades = false;
  promptLinesOnly = showLines = isOtherRef = displayLines = false;
  lowerBrCuttoff = age = 0.0;
}//void ReferenceLineInfo::reset()


bool ReferenceLineInfo::empty() const
{
  return (!nuclide && !element
          && reactionGammas.empty() && otherRefLines.empty());
}//bool empty() const


bool ReferenceLineInfo::operator==( const ReferenceLineInfo &rhs ) const
{
  return (nuclide==rhs.nuclide
          && element==rhs.element
          && labelTxt==rhs.labelTxt );
}//ReferenceLineInfo::operator==


void ReferenceLineInfo::toJson( string &json ) const
{
  // For reference, for Th232, the JSON returned is about 32 kb, and U238 is 70.5 kb (just gamma and xray).
  //  TODO: The "decay" for each line could be specified in a separate map, so like 'Ba133 to Cs133 via Electron Capture' isnt included in the JSON a bunch of times.  could also change "particle" to "p", "decay" to "d", and particle values from "gamma", "xray", etc, to "g", "x", etc.
  
  //const size_t color_index = (num % ns_numColors);
  
  std::stringstream jsonstrm;
  
  jsonstrm << "{\"color\":\"" << (lineColor.isDefault() ? "#0000FF" : lineColor.cssText(false) ) << "\","
  << "\"parent\":\"" << parentLabel() << "\",";
  
  if( promptLinesOnly )
    jsonstrm << "\"prompt\":true,";
  if( !otherRefLines.empty() )
    jsonstrm << "\"age\":\"Primordial\",";
  else
    jsonstrm << "\"age\":" << jsQuote( PhysicalUnits::printToBestTimeUnits(age) ) << ",";
  
  if( !detectorName.empty() )
    jsonstrm << "\"detector\":" << jsQuote( detectorName ) << ",";
  
  if( !shieldingName.empty() && (shieldingThickness > 0.1*PhysicalUnits::um ) )
  {
    const string name = jsQuote( shieldingName );
    string thickness = PhysicalUnits::printToBestLengthUnits( shieldingThickness );
    thickness = jsQuote( thickness );
    
    jsonstrm << "\"shielding\":" << name << ",";
    jsonstrm << "\"shieldingThickness\":" << thickness << ",";
  }//if( !shieldingName.empty() )
  
  //jsonstrm << "particleSf:{";
  //for( map<string,double>::const_iterator iter = particle_sf.begin();
  //     iter != particle_sf.end(); ++iter )
 // {
 //   jsonstrm << ((iter==particle_sf.begin()) ? "" : ",") << jsQuote(iter->first)
 //            << ":" << iter->second;
 // }
 // jsonstrm << "},";
  
  jsonstrm << "\"lines\":[";
  
  bool printed = false;
  char intensity_buffer[32] = { '\0' };
  
  // Round to the nearest 10 eV; probably the extent to which any data useful, or even good to
  const auto round_energy = []( const double e ) -> double { return std::round(100.0*e)/100.0; };
  
  for( size_t i = 0; i < energies.size(); ++i )
  {
    if( intensities[i] <= 0.0 )
      continue;
    
    const double energy = round_energy( energies[i] );
    
    // There are situations where two lines have either the exact same energies, or super-close
    //  energies, and the spectrum chart is not particularly smart about this, so we effectively
    //  lose some amplitude, so we'll combine them here.
    // However, this additional munging has a overhead for a rare edge-case, so we'll split
    //  the code-paths, even though this adds code...
    // For U238, there are 39 pairs of lines that get combined, and 1 triplet of lines combined.
    
    // We will assume entries are sorted by energy, which is only guaranteed when
    //  ReferencePhotopeakDisplay::updateDisplayChange() sets the data - but also, I think this
    //  is the only place that sets the data!

    // TODO: This current way of doing things wil sum a Cascade sum with a gamma (e.g., the 387.8 
    //       keV of U235), which is probably not the right thing to do because it then shows up as
    //       giant gamma on the chart, which is deceptive
    const bool next_gamma_close = (((i+1) < energies.size()) && (round_energy(energies[i+1]) == energy) );
    
    if( next_gamma_close )
    {
      double intensity = 0.0;
      size_t num_combined = 0;
      // TODO: be a little more efficient than allocating strings in these sets...
      set<string> particles, decays, elements;
      for( size_t index = i; index < energies.size(); ++index )
      {
        const double this_energy = round_energy(energies[index]);
        if( this_energy != energy )
          break;
        
        intensity += intensities[index];
        if( !particlestrs[index].empty() )
          particles.insert( particlestrs[index] );
        if( !decaystrs[index].empty() )
          decays.insert( decaystrs[index] );
        if( !elementstrs[index].empty() )
          elements.insert( elementstrs[index] );
        
        num_combined += 1;
      }//for( loop over to find all energies to cluster together )
      
//      if( particles.size() > 1 || decays.size() > 1 || elements.size() > 1 )
//        cout << "Combined " << num_combined << " lines for " << energy << " keV" << endl
//             << "\tAnd had " << particles.size() << " particles, " << decays.size()
//             << " decays, and " << elements.size() << " elements" << endl;
      
      assert( num_combined != 0 ); //could tighten this up to (num_combined > 1)
      
      if( IsNan(intensity) || IsInf(intensity) )
        snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
      else
        snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", intensity );
      
      auto combine_strs = []( const set<string> &strs ) -> string {
        if( strs.empty() )
          return jsQuote("");
        
        if( strs.size() == 1 )
          return jsQuote( *begin(strs) );
        
        string answer;
        for( const auto &s : strs )
          answer += (answer.empty() ? "" : ", ") + s;
        return jsQuote( answer );
      };//combine_strs lambda
      
      jsonstrm << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
      if( !particles.empty() )
        jsonstrm << ",\"particle\":" << combine_strs(particles);
      if( !decays.empty() )
        jsonstrm << ",\"decay\":" << combine_strs(decays);
      if( !elements.empty() )
        jsonstrm << ",\"el\":" << combine_strs(elements);
      
      // Now increment 'i' so we'll skip over these lines we've already covered.
      i += (num_combined >= 1) ? (num_combined - 1) : size_t(0);
    }else
    {
      if( IsNan(intensities[i]) || IsInf(intensities[i]) )
        snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
      else
        snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", intensities[i] );
      
      jsonstrm << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
      if( !particlestrs[i].empty() )
        jsonstrm << ",\"particle\":" << jsQuote(particlestrs[i]);
      if( !decaystrs[i].empty() )
        jsonstrm << ",\"decay\":" << jsQuote(decaystrs[i]);
      if( !elementstrs[i].empty() )
        jsonstrm << ",\"el\":" << jsQuote(elementstrs[i]);
    }//if( next gamma line is close ) / else
    
    jsonstrm << "}";
    
    printed = true;
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  jsonstrm <<"]}";
  
  json += jsonstrm.str();
  
#if( DEV_REF_LINE_UPGRADE_20221212 )
  {
    std::stringstream jsons;
    
    jsons << "{\"color\":\"" << (lineColor.isDefault() ? "#0000FF" : lineColor.cssText(false)) << "\","
    << "\"parent\":\"" << parentLabel() << "\",";
    
    if( m_input.m_promptLinesOnly )
      jsons << "\"prompt\":true,";
    if( m_source_type == ReferenceLineInfo::SourceType::Background )
      jsons << "\"age\":\"Primordial\",";
    else // DEV_REF_LINE_UPGRADE_20221212 TODO: after comparison to old, only due age for nuclides
      jsons << "\"age\":" << jsQuote( PhysicalUnits::printToBestTimeUnits(age) ) << ",";
    
    if( !detectorName.empty() )
      jsons << "\"detector\":" << jsQuote( m_input.m_detector_name ) << ",";
    
    if( !m_input.m_shielding_name.empty() && (m_input.m_shieldingThickness > 0.1*PhysicalUnits::um ) )
    {
      const string name = jsQuote( m_input.m_shielding_name );
      string thickness = PhysicalUnits::printToBestLengthUnits( m_input.m_shieldingThickness );
      thickness = jsQuote( thickness );
      
      jsons << "\"shielding\":" << name << ",";
      jsons << "\"shieldingThickness\":" << thickness << ",";
    }//if( !shieldingName.empty() )
    
    //jsons << "particleSf:{";
    //for( map<string,double>::const_iterator iter = particle_sf.begin();
    //    iter != particle_sf.end(); ++iter )
    //{
    //  jsons << ((iter==particle_sf.begin()) ? "" : ",") << jsQuote(iter->first)
    //  << ":" << iter->second;
    //}
    //jsons << "},"
    
    jsons << "\"lines\":[";
    
    bool printed = false;
    char intensity_buffer[32] = { '\0' };
    
    // Round to the nearest 10 eV; probably the extent to which any data useful, or even good to
    const auto round_energy = []( const double e ) -> double { return std::round(100.0*e)/100.0; };
    
    for( size_t index = 0; index < m_ref_lines.size(); ++index )
    {
      const RefLine &line = m_ref_lines[index];
      if( line.m_normalized_intensity <= 0.0 )
        continue;
      
      const double energy = round_energy( line.m_energy );
      
      // There are situations where two lines have either the exact same energies, or super-close
      //  energies, and the spectrum chart is not particularly smart about this, so we effectively
      //  lose some amplitude, so we'll combine them here.
      // However, this additional munging has a overhead for a rare edge-case, so we'll split
      //  the code-paths, even though this adds code...
      // For U238, there are 39 pairs of lines that get combined, and 1 triplet of lines combined.
      
      // We will assume entries are sorted by energy, which is only guaranteed when
      //  ReferencePhotopeakDisplay::updateDisplayChange() sets the data - but also, I think this
      //  is the only place that sets the data!
      
      // TODO: This current way of doing things will sum a Cascade sum with a gamma (e.g., the 387.8
      //       keV of U235), which is probably not the right thing to do because it then shows up as
      //       giant gamma on the chart, which is deceptive
      const bool next_gamma_close = (((index+1) < m_ref_lines.size())
                                     && (round_energy(m_ref_lines[index+1].m_energy) == energy)
                                     //&& (m_ref_lines[index+1].m_source_type == line.m_source_type)
                                     //&& (m_ref_lines[index+1].m_particle_type == line.m_particle_type)
                                     //&& (m_ref_lines[index+1].m_particle_type == RefLine::Particle::Gamma)
                                     );
      
      if( next_gamma_close )
      {
        double intensity = 0.0;
        size_t num_combined = 0;
        // TODO: be a little more efficient than allocating strings in these sets...
        set<string> particles, decays, elements;
        for( size_t inner_index = index; inner_index < m_ref_lines.size(); ++inner_index )
        {
          const RefLine &inner_line = m_ref_lines[inner_index];
          
          const double this_energy = round_energy(inner_line.m_energy);
          if( this_energy != energy )
            break;
          
          // TODO: skip over zero intensity lines - leaving commented out for comparison
          //if( inner_line.m_normalized_intensity <= 0.0 )
          //{
            //num_combined += 1;
            //continue;
          //}
            
          intensity += inner_line.m_normalized_intensity;
          if( !inner_line.m_particlestr.empty() )
            particles.insert( inner_line.m_particlestr );
          if( !inner_line.m_decaystr.empty() )
            decays.insert( inner_line.m_decaystr );
          if( !inner_line.m_elementstr.empty() )
            elements.insert( inner_line.m_elementstr );
          
          num_combined += 1;
        }//for( loop over inner_index to find all energies to cluster together )
        
        assert( num_combined != 0 ); //could tighten this up to (num_combined > 1)
        
        if( IsNan(intensity) || IsInf(intensity) )
          snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
        else
          snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", intensity );
        
        auto combine_strs = []( const set<string> &strs ) -> string {
          string answer;
          for( const auto &s : strs )
            answer += (answer.empty() ? "" : ", ") + s;
          return answer;
        };//combine_strs lambda
        
        jsons << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
        if( !particles.empty() )
          jsons << ",\"particle\":" << jsQuote(combine_strs(particles));
        if( !decays.empty() )
          jsons << ",\"decay\":" << jsQuote(combine_strs(decays));
        if( !elements.empty() )
          jsons << ",\"el\":" << jsQuote(combine_strs(elements));
        
        // Now increment 'i' so we'll skip over these lines we've already covered.
        index += (num_combined >= 1) ? (num_combined - 1) : size_t(0);
      }else
      {
        if( IsNan(line.m_normalized_intensity) || IsInf(line.m_normalized_intensity) )
          snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
        else
          snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", line.m_normalized_intensity );
        
        jsons << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
        if( !line.m_particlestr.empty() )
          jsons << ",\"particle\":" << jsQuote(line.m_particlestr);
        if( !line.m_decaystr.empty() )
          jsons << ",\"decay\":" << jsQuote(line.m_decaystr);
        if( !line.m_elementstr.empty() )
          jsons << ",\"el\":" << jsQuote(line.m_elementstr);
      }//if( next gamma line is close ) / else
      
      jsons << "}";
      
      printed = true;
    }//for( size_t index = 0; index < m_ref_lines.size(); ++index )
    
    jsons <<"]}";
    
    if( jsons.str() != jsonstrm.str() )
    {
      cout << "New and old JSON string dont match:\n\nNew: " << jsons.str()
      << "\n\nOld: " << jsonstrm.str() << endl;
    }else
    {
      cout << "\n\nNew and old JSON match! Len=" << jsonstrm.str().size() << "\n\n";
    }
    
  }
#endif
}//std::string toJson( const ReferenceLineInfo &displnuc )


void ReferenceLineInfo::sortByEnergy()
{
  const size_t len = energies.size();
  
  if( len != energies.size() || len != intensities.size()
     || len != particlestrs.size() || len != decaystrs.size()
     || len != elementstrs.size()
#if( DEV_REF_LINE_UPGRADE_20221212 )
 //   || len != m_ref_lines.size() 
#endif
    )
    throw runtime_error( "ReferenceLineInfo::sortByEnergy(): inconsistent input size" );
  
  vector<size_t> sort_indices( energies.size() );
  for( size_t i = 0; i < len; ++i )
    sort_indices[i] = i;
  
#if( DEV_REF_LINE_UPGRADE_20221212 )
// THis is only necassary for comparing new and old values
#if _MSC_VER
#pragma message("Doing sort by energy and intensity for comparison only")
#else
  #warning "Doing sort by energy and intensity for comparison only"
#endif
#endif
  std::sort( sort_indices.begin(), sort_indices.end(),
            index_compare_assend<vector<double>&>(energies, intensities) );
  
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

#if( DEV_REF_LINE_UPGRADE_20221212 )
  std::sort( begin( m_ref_lines ), end( m_ref_lines ), 
    []( const RefLine &lhs, const RefLine &rhs ) -> bool {
    if( lhs.m_energy == rhs.m_energy )
      return lhs.m_normalized_intensity < rhs.m_normalized_intensity;
    return lhs.m_energy < rhs.m_energy;
    } );
#endif
  
#if( PERFORM_DEVELOPER_CHECKS )
  for( size_t i = 1; i < len; ++i )
  {
    if( energies[i-1] > energies[i] )
    {
      log_developer_error( __func__, "Error sorting reference photopeaks" );
      break;
    }
  }//for( size_t i = 0; i < len; ++i )
#endif
}//void sortByEnergy()


std::string ReferenceLineInfo::parentLabel() const
{
  if( nuclide )
    return nuclide->symbol;
  
  if( element )
    return element->symbol;
  
  if( !reactionGammas.empty() )
    return reactionsTxt;
  
  if( !otherRefLines.empty() )
  {
    if( (labelTxt.size() >= 4)
      && ((labelTxt[0] == 'b') || (labelTxt[0] == 'B'))
      && ((labelTxt[0] == 'a') || (labelTxt[0] == 'A'))
      && ((labelTxt[0] == 'c') || (labelTxt[0] == 'C'))
      && ((labelTxt[0] == 'k') || (labelTxt[0] == 'K')) )
    {
      return "Background";
    }
  }//if( !otherRefLines.empty() )

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
  
  node = base_node->first_node( "LineColor", 9 );
  if( node && node->value_size() )
  {
    try{
      lineColor = Wt::WColor( node->value() );
    }catch(...){
      lineColor = Wt::WColor();
    }
  }else
  {
    lineColor = Wt::WColor();
  }
  
  node = base_node->first_node( "DisplayLines", 12 );
  if( node && node->value_size() )
    displayLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowGammas", 10 );
  if( node && node->value_size() )
    showGammas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowXrays", 9 );
  if( node && node->value_size() )
    showXrays = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowAlphas", 10 );
  if( node && node->value_size() )
    showAlphas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowBetas", 9 );
  if( node && node->value_size() )
    showBetas = (node->value()[0] == '1');

  node = base_node->first_node("ShowCascades", 12);
  if (node && node->value_size())
    showCascades = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowLines", 9 );
  if( node && node->value_size() )
    showLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "PromptLinesOnly", 15 );
  if( node && node->value_size() )
    promptLinesOnly = (node->value()[0] == '1');
  
  node = base_node->first_node( "IsBackground", 12 );  // IsBackground is holdover from pre 20221209 code
  if( node && node->value_size() )
    isOtherRef = (node->value()[0] == '1');

  node = base_node->first_node( "IsReaction", 10 );
  if( node && node->value_size() )
    isReaction = (node->value()[0] == '1');
  
  node = base_node->first_node( "Age", 3 );
  if( node && node->value_size() && !(stringstream(node->value()) >> age) )
    throw runtime_error( "Ivalid age parameter" );
  
  node = base_node->first_node( "LowestBranchRatio", 17 );
  if( node && node->value_size() && !(stringstream(node->value())>> lowerBrCuttoff) )
    throw runtime_error( "Ivalid LowerBrCuttoff parameter" );
  
  node = base_node->first_node( "ShieldingName", 13 );
  if( node && node->value_size() && node->value_size() )
    shieldingName = node->value();

  node = base_node->first_node( "ShieldingThickness", 18 );
  if( node && node->value_size() && !(stringstream(node->value())>> shieldingThickness) )
    throw runtime_error( "Ivalid ShieldingThickness parameter" );

  node = base_node->first_node( "DetectorName", 12 );
  if( node && node->value_size() )
    detectorName = node->value();
  
  node = base_node->first_node( "Nuclide", 7 );
  if( node && node->value_size() )
    labelTxt = node->value();
  

  if( isOtherRef )
  {
    node = base_node->first_node( "OtherRefLines", 13 );
    if( !node )
    {
      // Kept here for backwards compatibility with states saved by pre 20221209 code
      for( const OtherRefLine &bl : BackgroundLines )
      {
        if( (std::get<3>( bl ) != OtherRefLineType::BackgroundXRay) || showXrays )
          otherRefLines.push_back( bl );
      }//for( const BackgroundLine &bl : BackgroundLines )
    }else
    {
      for( const rapidxml::xml_node<char> *ref_node = node->first_node( "RefLine", 7 );
        ref_node; ref_node = ref_node->next_sibling( "RefLine", 7 ) )
      {
        const rapidxml::xml_node<char> *energy_node = ref_node->first_node( "Energy", 6 );
        const rapidxml::xml_node<char> *br_node = ref_node->first_node( "Intensity", 9 );
        const rapidxml::xml_node<char> *symbol_node = ref_node->first_node( "Symbol", 6 );
        const rapidxml::xml_node<char> *type_node = ref_node->first_node( "Type", 4 );
        const rapidxml::xml_node<char> *desc_node = ref_node->first_node( "Desc", 4 );

        if( !energy_node || !br_node || !type_node 
          || !energy_node->value() || !br_node->value()  || !type_node->value() )
          throw runtime_error( "Invalid other ref line element" );
        
        OtherRefLine line;
        if( !(stringstream( energy_node->value() ) >> get<0>(line)) )
          throw runtime_error( "Invalid other ref line energy" );
        if( !(stringstream( br_node->value() ) >> get<1>( line )) )
          throw runtime_error( "Invalid other ref line B.R." );

        get<2>( line ) = (symbol_node && symbol_node->value_size()) ? symbol_node->value() : "";
        get<3>( line ) = other_ref_line_type_from_str( type_node->value() );
        get<4>( line ) = (desc_node && desc_node->value_size()) ? desc_node->value() : "";

        otherRefLines.push_back( line );
      }
    }//if( !node )
  }//if( isOtherRef )
  
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
      double energy = 0.0, intensity = 0.0;
      energy_node = node->first_node( "Energy", 6 );
      br_node = node->first_node( "Intensity", 9 );
      
      if( !energy_node || !br_node || !energy_node->value() || !br_node->value() )
        throw runtime_error( "Invalid line element" );
      
      if( !(stringstream(energy_node->value()) >> energy) )
        throw runtime_error( "Invalid line element energy" );
        
      // Some intensities can (what I presume to be) non-normal (i.e, like 1.0E-312), so we'll
      //  just set it to zero
      if( !(stringstream(br_node->value()) >> intensity) )
        intensity = 0.0;
      
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
      {
        //got an error on iPhone 20190716 with an error for str=="xray", but couldnt follow through on diagnosing.
        //if( str == SandiaDecay::to_str(SandiaDecay::ProductType::XrayParticle)
        //   && particle_sf.count(SandiaDecay::to_str(SandiaDecay::ProductType::GammaParticle)) )
        //{
        //  particle_sf[str] = particle_sf[SandiaDecay::to_str(SandiaDecay::ProductType::GammaParticle)];
        //}else
        //{
          throw runtime_error( "Particle type does not have a SF: " + str );
        //}
      }
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
  
  name = "Nuclide";
  value = doc->allocate_string( labelTxt.c_str() );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  if( !lineColor.isDefault() )
  {
    name = "LineColor";
    value = doc->allocate_string( lineColor.cssText(false).c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }
  
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
  
  name = "ShowCascades";
  value = (showCascades ? "1" : "0");
  node = doc->allocate_node(rapidxml::node_element, name, value);
  base_node->append_node(node);

  name = "PromptLinesOnly";
  value = (promptLinesOnly ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowLines";
  value = (showLines ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "IsBackground";
  value = (isOtherRef ? "1" : "0");
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
  
  
  if( !otherRefLines.empty() )
  {
    // Pre 20221209 code will not read in the <OtherRefLines> node, and instead load the 
    //  BackgroundLines into this->otherRefLines - I guess this will be fine...
    name = "OtherRefLines";
    rapidxml::xml_node<char> *other_ref_lines = doc->allocate_node( rapidxml::node_element, name );
    base_node->append_node( other_ref_lines );

    for( const OtherRefLine &bl : otherRefLines )
    {
      using rapidxml::node_element;

      name = "RefLine";
      rapidxml::xml_node<char> *line_node = doc->allocate_node( rapidxml::node_element, name );
      other_ref_lines->append_node( line_node );

      snprintf( buffer, sizeof( buffer ), "%1.7g", get<0>(bl) );
      value = doc->allocate_string( buffer );
      node = doc->allocate_node( node_element, "Energy", value );
      line_node->append_node( node );

      snprintf( buffer, sizeof( buffer ), "%1.7g", get<1>( bl ) );
      value = doc->allocate_string( buffer );
      node = doc->allocate_node( node_element, "Intensity", value );
      line_node->append_node( node );

      value = doc->allocate_string( get<2>( bl ).c_str() );
      node = doc->allocate_node( node_element, "Symbol", value );
      line_node->append_node( node );

      value = to_str( get<3>( bl ) );
      node = doc->allocate_node( node_element, "Type", value );
      line_node->append_node( node );

      value = doc->allocate_string( get<4>( bl ).c_str() );
      node = doc->allocate_node( node_element, "Desc", value );
      line_node->append_node( node );
    }//for( const OtherRefLine &bl : otherRefLines )
  }//if( !otherRefLines.empty() )


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



