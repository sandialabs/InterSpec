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

#include <Wt/WWebWidget>  //For quoting strings only

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GammaInteractionCalc.h"


using namespace std;

/**
 version 0 and 1: serialized the ReferenceLineInfo struct.
 version 2: implemented 20221220 serializes just RefLineInput - you can then use this to generate ReferenceLineInfo
 */
const int RefLineInput::sm_xmlSerializationVersion = 2;


namespace
{
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
  //TODO: do some better way of indicating D.E. than a string match
  OtherRefLine( 1592.51f, 0.010f,   "Th232 D.E. 2614 keV", OtherRefLineType::OtherBackground,    "" ),
  OtherRefLine( 1620.74f, 0.0151f,  "Bi212", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1630.63f, 0.0151f,  "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1686.09f, 0.00095f, "Ac228", OtherRefLineType::Th232Series,        "" ),
  OtherRefLine( 1661.28f, 0.0112f,  "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1729.6f,  0.02843f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1764.49f, 0.1528f,  "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1838.36f, 0.00346f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  OtherRefLine( 1847.42f, 0.02023f, "Bi214", OtherRefLineType::Ra226Series,        "" ),
  //TODO: do some better way of indicating S.E. than a string match
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



RefLineInput::RefLineInput()
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
  m_shielding_thickness(),
  m_shielding_an(),
  m_shielding_ad(),
  m_shielding_att()
{
}

void RefLineInput::reset()
{
  *this = RefLineInput{};
}

ReferenceLineInfo::RefLine::RefLine()
  : m_energy( 0.0 ),
  m_normalized_intensity( 0.0 ),
  m_drf_factor( 1.0f ),
  m_shield_atten( 1.0f ),
  m_particle_sf_applied( 1.0f ),
  m_decaystr(),
  m_decay_intensity( 0.0 ),
  m_particle_type( ReferenceLineInfo::RefLine::Particle::Gamma ),
  m_parent_nuclide( nullptr ),
  m_transition( nullptr ),
  m_source_type( ReferenceLineInfo::RefLine::RefGammaType::Normal ),
  m_element( nullptr ),
  m_reaction( nullptr )
{
}


const std::string &ReferenceLineInfo::RefLine::particlestr() const
{
  static const string cascade_sum = "cascade-sum";
  static const string gamma_sum = "sum-gamma";
  static const string single_escape = "S.E.";
  static const string double_escape = "D.E.";
  
  static const string alpha = "alpha";
  static const string beta = "beta";
  static const string gamma = "gamma";
  static const string xray = "xray";
  static const string other = "invalid";
  
  
  switch( m_source_type )
  {
    case RefGammaType::Normal:
    case RefGammaType::Annihilation:
      break;
      
    case RefGammaType::SingleEscape:        return single_escape;
    case RefGammaType::DoubleEscape:        return double_escape;
    case RefGammaType::CoincidenceSumPeak:  return cascade_sum;
    case RefGammaType::SumGammaPeak:        return gamma_sum;
  }//switch( m_source_type )
  
  
  switch( m_particle_type )
  {
    case ReferenceLineInfo::RefLine::Particle::Alpha: return alpha;
    case ReferenceLineInfo::RefLine::Particle::Beta:  return beta;
    case ReferenceLineInfo::RefLine::Particle::Gamma: return gamma;
    case ReferenceLineInfo::RefLine::Particle::Xray:  return xray;
  }//switch( line.m_particle_type )
  
  assert( 0 );
  return other;
}//particlestr(...)

ReferenceLineInfo::ReferenceLineInfo()
{
  reset();
}

void ReferenceLineInfo::reset()
{
  m_ref_lines.clear();
  m_input_warnings.clear();
  m_validity = InputValidity::Blank;
  m_has_coincidences = false;
  m_input = RefLineInput();
  m_source_type = ReferenceLineInfo::SourceType::None;

  m_nuclide = nullptr;
  m_element = nullptr;
  m_reactions.clear();
}//void ReferenceLineInfo::reset()	  


void ReferenceLineInfo::toJson( string &json ) const
{
  // For reference, for Th232, the JSON returned is about 32 kb, and U238 is 70.5 kb (just gamma and xray).
  //  TODO: The "decay" for each line could be specified in a separate map, so like 'Ba133 to Cs133 via Electron Capture' isnt included in the JSON a bunch of times.  could also change "particle" to "p", "decay" to "d", and particle values from "gamma", "xray", etc, to "g", "x", etc.
  
  std::stringstream jsons;
  
  // We will put individual lines descriptions in a separate array, and give the index into
  //  the array that line should go; this is because decays like U238 have many gammas for
  //  many of the transitions, so we can save a decent about of space.
  //  For {U238,Th232,Ba133}, the JSON size goes from {71, 37, 1} kb, to {52,27,1} kb
  vector<string> decay_strs;
  map<string,size_t> decay_str_indexs;
  
  jsons << "{\"color\":\"" << (m_input.m_color.isDefault() ? "#0000FF" : m_input.m_color.cssText(false)) << "\","
  << "\"parent\":\"" << m_input.m_input_txt << "\",";
  
  if( m_input.m_promptLinesOnly )
    jsons << "\"prompt\":true,";
  if( m_source_type == ReferenceLineInfo::SourceType::Background )
    jsons << "\"age\":\"Primordial\",";
  else if( !m_input.m_age.empty() )
    jsons << "\"age\":" << jsQuote( m_input.m_age ) << ",";
  
  if( !m_input.m_detector_name.empty() )
    jsons << "\"detector\":" << jsQuote( m_input.m_detector_name ) << ",";
  
  if( !m_input.m_shielding_name.empty() && !m_input.m_shielding_thickness.empty() )
  {
    assert( m_input.m_shielding_an.empty() );
    assert( m_input.m_shielding_ad.empty() );
#ifndef NDEBUG
    PhysicalUnits::stringToDistance(m_input.m_shielding_thickness);
#endif

    jsons << "\"shielding\":" << jsQuote( m_input.m_shielding_name ) << ",";
    jsons << "\"shieldingThickness\":" << jsQuote( m_input.m_shielding_thickness ) << ",";
  }else if( !m_input.m_shielding_an.empty() && !m_input.m_shielding_ad.empty() )
  {
    assert( m_input.m_shielding_name.empty() );
    assert( m_input.m_shielding_thickness.empty() );
#ifndef NDEBUG
    double dummy;
    assert( (stringstream(m_input.m_shielding_an) >> dummy) );
    assert( (stringstream(m_input.m_shielding_ad) >> dummy) );
#endif
    
    const string shield = "AN=" + m_input.m_shielding_an + ", AD=" + m_input.m_shielding_ad + " g/cm2";
    jsons << "\"shielding\":" << jsQuote( shield ) << ",";
  }

  
  jsons << "\"lines\":[";
  
  bool printed = false;
  char intensity_buffer[32] = { '\0' };
  
  // Round to the nearest 10 eV; probably the extent to which any data useful, or even good to
  const auto round_energy = []( const double e ) -> double { return std::round(100.0*e)/100.0; };
  
  for( size_t index = 0; index < m_ref_lines.size(); ++index )
  {
    const RefLine &line = m_ref_lines[index];
    if( line.m_normalized_intensity <= std::numeric_limits<float>::min() )  //numeric_limits<float>::min()==1.17549e-38
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
      set<string> particles, decays;
      for( size_t inner_index = index; inner_index < m_ref_lines.size(); ++inner_index )
      {
        const RefLine &inner_line = m_ref_lines[inner_index];
        
        const double this_energy = round_energy(inner_line.m_energy);
        if( this_energy != energy )
          break;
        
        if( inner_line.m_normalized_intensity <= 0.0 )
        {
          num_combined += 1;
          continue;
        }
        
        intensity += inner_line.m_normalized_intensity;
        particles.insert( inner_line.particlestr() );
        if( !inner_line.m_decaystr.empty() )
          decays.insert( inner_line.m_decaystr );
        
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
      {
        const string decay_str = combine_strs(decays);
        auto decay_str_iter = decay_str_indexs.find(decay_str);
        if( decay_str_iter == end(decay_str_indexs) )
        {
          decay_str_iter = decay_str_indexs.insert( {decay_str, decay_strs.size()} ).first;
          decay_strs.push_back(decay_str);
        }
        
        jsons << ",\"decayind\":" << decay_str_iter->second;
      }
      
      // Now increment 'i' so we'll skip over these lines we've already covered.
      index += (num_combined >= 1) ? (num_combined - 1) : size_t(0);
    }else
    {
      if( IsNan(line.m_normalized_intensity) || IsInf(line.m_normalized_intensity) )
        snprintf( intensity_buffer, sizeof(intensity_buffer), "0" );
      else
        snprintf( intensity_buffer, sizeof(intensity_buffer), "%.3g", line.m_normalized_intensity );
      
      jsons << (printed ? "," : "") << "{\"e\":" << energy << ",\"h\":" << intensity_buffer;
      jsons << ",\"particle\":" << jsQuote( line.particlestr() );
      
      if( !line.m_decaystr.empty() )
      {
        auto decay_str_iter = decay_str_indexs.find(line.m_decaystr);
        if( decay_str_iter == end(decay_str_indexs) )
        {
          decay_str_iter = decay_str_indexs.insert( {line.m_decaystr, decay_strs.size()} ).first;
          decay_strs.push_back(line.m_decaystr);
        }
        
        jsons << ",\"desc_ind\":" << decay_str_iter->second;
      }//if( !line.m_decaystr.empty() )
    }//if( next gamma line is close ) / else
    
    jsons << "}";
    
    printed = true;
  }//for( size_t index = 0; index < m_ref_lines.size(); ++index )
  jsons <<"]";

  
  jsons <<", \"desc_strs\":[";
  for( size_t i = 0; i < decay_strs.size(); ++i )
    jsons << (i ? "," : "") << jsQuote( decay_strs[i] );
  jsons <<"]}";
  
  json += jsons.str();
}//std::string toJson( const ReferenceLineInfo &displnuc )


void ReferenceLineInfo::sortByEnergy()
{
  std::sort( begin( m_ref_lines ), end( m_ref_lines ), 
    []( const RefLine &lhs, const RefLine &rhs ) -> bool {
    if( lhs.m_energy == rhs.m_energy )
      return lhs.m_normalized_intensity < rhs.m_normalized_intensity;
    return lhs.m_energy < rhs.m_energy;
    } );
}//void sortByEnergy()


void RefLineInput::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  // The schema used was for #ReferenceLineInfo, before 20221220, so its
  //  not exactly what you would write from scratch, but good-enough
  //  considering it gives us backward compatibility, mostly.
  // <DisplayedSource> is for pre-20221220, and <RefLineInput> is
  // post-20221220 (so that old-code wont even try to deserialize the stuff).
  
  if( !base_node
     || ((base_node->name() != string("DisplayedSource"))
         && (base_node->name() != string("RefLineInput"))) )
    throw runtime_error( "Invalid base node for DisplayedSource" );
  
  static_assert( sm_xmlSerializationVersion == 2, "You need to update xml serialization" );
  
  int version;
  const rapidxml::xml_attribute<char> * const version_attr = base_node->first_attribute( "version", 7 );
  if( !version_attr || !version_attr->value()
     || !(stringstream(version_attr->value()) >> version)
     || (version != 0 && version != 1 && version != 2) )
    throw runtime_error( "Missing or invalid DisplayedSource version" );
 
  reset();
  RefLineInput &input = *this;
  
  const rapidxml::xml_node<char> *node = base_node->first_node( "Nuclide", 7 );
  if( node && node->value_size() )
    input.m_input_txt = node->value();
  
  node = base_node->first_node( "Element", 7 );
  if( node && node->value_size() && input.m_input_txt.empty() )
    input.m_input_txt = node->value();
  
  node = base_node->first_node( "LineColor", 9 );
  if( node && node->value_size() )
    try{ input.m_color = Wt::WColor( node->value() ); }catch(...){ }
  
  //node = base_node->first_node( "DisplayLines", 12 );
  //if( node && node->value_size() )
  //  displayLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowGammas", 10 );
  if( node && node->value_size() )
    input.m_showGammas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowXrays", 9 );
  if( node && node->value_size() )
    input.m_showXrays = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowAlphas", 10 );
  if( node && node->value_size() )
    input.m_showAlphas = (node->value()[0] == '1');
  
  node = base_node->first_node( "ShowBetas", 9 );
  if( node && node->value_size() )
    input.m_showBetas = (node->value()[0] == '1');
  
  node = base_node->first_node("ShowCascades", 12);
  if (node && node->value_size())
    input.m_showCascades = (node->value()[0] == '1');
  
  //node = base_node->first_node( "ShowLines", 9 );
  //if( node && node->value_size() )
  //  showLines = (node->value()[0] == '1');
  
  node = base_node->first_node( "PromptLinesOnly", 15 );
  if( node && node->value_size() )
    input.m_promptLinesOnly = (node->value()[0] == '1');
  
  //node = base_node->first_node( "IsBackground", 12 );
  //if( node && node->value_size() )
  //  isOtherRef = (node->value()[0] == '1');
  
  //node = base_node->first_node( "IsReaction", 10 );
  //if( node && node->value_size() )
  //  isReaction = (node->value()[0] == '1');
  
  node = base_node->first_node( "Age", 3 );
  if( node && node->value_size() )
    input.m_age = node->value();
  
  node = base_node->first_node( "LowestBranchRatio", 17 );
  if( node && node->value_size() && !(stringstream(node->value()) >> input.m_lower_br_cutt_off) )
    throw runtime_error( "Invalid LowerBrCuttoff parameter" );
  
  node = base_node->first_node( "ShieldingName", 13 );
  if( node && node->value_size() )
  {
    input.m_shielding_name = node->value();
    
    // Pre 20221220, generic shielding would look like "26, 12.1 g/cm2", so we'll convert
    //  this to AN and AD.
    if( (version < 2) && SpecUtils::icontains(input.m_shielding_name, "g/cm2") )
    {
      vector<string> parts;
      SpecUtils::split( parts, input.m_shielding_name, "," );
      if( parts.size() == 2 )
      {
        input.m_shielding_an = SpecUtils::trim_copy(parts[0]);
        input.m_shielding_ad = SpecUtils::trim_copy(parts[1]);
        double dummy;
        if( !(stringstream(input.m_shielding_an) >> dummy)
           || !(stringstream(input.m_shielding_ad) >> dummy) )
        {
          input.m_shielding_an = input.m_shielding_ad = "";
          cerr << "Old RefLineInput Shielding '" << input.m_shielding_name << "' was not valid AD, AN." << endl;
        }else
        {
          input.m_shielding_name = "";
        }
      }else
      {
        cerr << "\n\ninput.m_shielding_name_='" << input.m_shielding_name << "', isnt formed as expected" << endl;
      }
    }
    
    node = base_node->first_node( "ShieldingThickness", 18 );
    
    if( node && node->value_size() )
    {
      input.m_shielding_thickness = node->value();
      try
      {
        // Make sure a proper distance
        PhysicalUnits::stringToDistance(input.m_shielding_thickness);
      }catch( std::exception &e )
      {
        input.m_shielding_thickness = "";
        cerr << "Unexpected invalid shielding thickness when deserializing: " << e.what() << endl;
      }
    }//if( node && node->value_size() )
  }//if( node && node->value_size() && node->value_size() )
  
  node = base_node->first_node( "ShieldingAN", 11 );
  if( node && node->value_size() )
  {
    input.m_shielding_an = node->value();
  
    double dummy;
    if( !(stringstream(input.m_shielding_an) >> dummy) )
    {
      cerr << "Unexpected invalid AN (" << input.m_shielding_an << ")" << endl;
      input.m_shielding_an = input.m_shielding_ad = "";
    }else
    {
      assert( input.m_shielding_name.empty() );
      assert( input.m_shielding_thickness.empty() );
      
      input.m_shielding_name = input.m_shielding_thickness = "";
    }//
  }//if( node && node->value_size() )
  
  
  node = base_node->first_node( "ShieldingAD", 11 );
  if( node && node->value_size() )
  {
    input.m_shielding_ad = node->value();
    double dummy;
    if( !(stringstream(input.m_shielding_ad) >> dummy) )
    {
        cerr << "Unexpected invalid AD (" << input.m_shielding_ad << ")" << endl;
        input.m_shielding_an = input.m_shielding_ad = "";
    }else
    {
      assert( input.m_shielding_name.empty() );
      assert( input.m_shielding_thickness.empty() );
      input.m_shielding_name = input.m_shielding_thickness = "";
    }//
  }//if( node && node->value_size() )
  
  
  if( !input.m_shielding_an.empty() || !input.m_shielding_ad.empty() )
  {
    assert( input.m_shielding_name.empty() );
    assert( input.m_shielding_thickness.empty() );
    input.m_shielding_name.clear();
    input.m_shielding_thickness.clear();
  }
  
  if( !input.m_shielding_name.empty() || !input.m_shielding_thickness.empty() )
  {
    assert( input.m_shielding_an.empty() );
    assert( input.m_shielding_ad.empty() );
    input.m_shielding_an.clear();
    input.m_shielding_ad.clear();
  }
  
  node = base_node->first_node( "DetectorName", 12 );
  if( node && node->value_size() )
    input.m_detector_name = node->value();
}//void deSerialize( const rapidxml::xml_node<char> *node );


void RefLineInput::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  // See notes at top of #RefLineInput::deSerialize for history of the function.
  
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  char buffer[64];
  
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node, *lines_node, *sf_node,
  *part_sf_node, *energy_node, *br_node, *str_node;
  rapidxml::xml_attribute<char> *attr;
  
  name = "RefLineInput";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  static_assert( sm_xmlSerializationVersion == 2, "You need to update xml serialization" );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  snprintf( buffer, sizeof(buffer), "%i", sm_xmlSerializationVersion );
  value = doc->allocate_string( buffer );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  name = "Nuclide";
  value = doc->allocate_string( m_input_txt.c_str() );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  if( !m_color.isDefault() )
  {
    name = "LineColor";
    value = doc->allocate_string( m_color.cssText(false).c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }
  
  if( !m_age.empty() )
  {
    name = "Age";
    value = doc->allocate_string( m_age.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( age >= 0.0 )
  
  name = "LowestBranchRatio";
  snprintf( buffer, sizeof(buffer), "%g", m_lower_br_cutt_off );
  value = doc->allocate_string( buffer );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  //name = "DisplayLines";
  //value = (displayLines ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  name = "ShowGammas";
  value = (m_showGammas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowXrays";
  value = (m_showXrays ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowAlphas";
  value = (m_showAlphas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowBetas";
  value = (m_showBetas ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ShowCascades";
  value = (m_showCascades ? "1" : "0");
  node = doc->allocate_node(rapidxml::node_element, name, value);
  base_node->append_node(node);

  name = "PromptLinesOnly";
  value = (m_promptLinesOnly ? "1" : "0");
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  //name = "ShowLines";
  //value = (showLines ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  //name = "IsBackground";
  //value = (isOtherRef ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  //name = "IsReaction";
  //value = (isReaction ? "1" : "0");
  //node = doc->allocate_node( rapidxml::node_element, name, value );
  //base_node->append_node( node );
  
  if( !m_shielding_name.empty() )
  {
    assert( m_shielding_an.empty() );
    assert( m_shielding_ad.empty() );
    
    name = "ShieldingName";
    value = doc->allocate_string( m_shielding_name.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( shieldingName.size() )
  
  if( !m_shielding_thickness.empty() )
  {
    assert( m_shielding_an.empty() );
    assert( m_shielding_ad.empty() );
    
    name = "ShieldingThickness";
    value = doc->allocate_string( m_shielding_thickness.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }
  
  
  if( !m_shielding_an.empty() || !m_shielding_ad.empty() )
  {
    assert( m_shielding_name.empty() );
    assert( m_shielding_thickness.empty() );
    
#ifndef NDEBUG
    double dummy;
    assert( m_shielding_an.empty() || (stringstream(m_shielding_an) >> dummy) );
    assert( m_shielding_ad.empty() || (stringstream(m_shielding_ad) >> dummy) );
#endif
    
    if( !m_shielding_an.empty() )
    {
      name = "ShieldingAN";
      value = doc->allocate_string( m_shielding_an.c_str() );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      base_node->append_node( node );
    }
    
    if( !m_shielding_ad.empty() )
    {
      name = "ShieldingAD";
      value = doc->allocate_string( m_shielding_ad.c_str() );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      base_node->append_node( node );
    }
  }//if( !m_shielding_an.empty() && !m_shielding_ad.empty() )
  
  
  if( !m_detector_name.empty() )
  {
    name = "DetectorName";
    value = doc->allocate_string( m_detector_name.c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );
  }//if( detectorName.size() )
}//void ReferenceLineInfo::serialize( rapidxml::xml_node<char> *parent_node )


void RefLineInput::setShieldingAttFcn( const MaterialDB *db )
{
  // Check for generic shielding
  if( !m_shielding_an.empty() && !m_shielding_ad.empty() )
  {
    assert( m_shielding_name.empty() );
    assert( m_shielding_thickness.empty() );
    
    double atomic_number, areal_density;
    if( !(stringstream(m_shielding_an) >> atomic_number) )
      throw runtime_error( "Couldnt convert atomic number string '"
                           + m_shielding_an + "' to a floating point." );
    
    if( !(stringstream(m_shielding_ad) >> areal_density) )
      throw runtime_error( "Couldnt convert areal density string '"
                          + m_shielding_ad + "' to a floating point." );
    
    if( (atomic_number < MassAttenuation::sm_min_xs_atomic_number)
       || (atomic_number > MassAttenuation::sm_max_xs_atomic_number) )
      throw runtime_error( "Invalid atomic value '" + m_shielding_an + "'." );
    
    if( areal_density < 0.0 )
      throw runtime_error( "Invalid areal density value '" + m_shielding_ad + "'." );
    
    areal_density *= (PhysicalUnits::gram / PhysicalUnits::cm2);
    
    m_shielding_att = [atomic_number, areal_density]( float energy ) -> double {
      const double att_coef = GammaInteractionCalc::transmition_coefficient_generic( atomic_number, areal_density, energy );
      return exp( -1.0 * att_coef );
    };
    
    return;
  }//if( !m_shielding_an.empty() && !m_shielding_ad.empty() )
  
  assert( m_shielding_an.empty() );
  assert( m_shielding_ad.empty() );
  
  if( m_shielding_name.empty() || m_shielding_thickness.empty() )
  {
    m_shielding_att = nullptr;
    return;
  }
  
  if( m_shielding_name.empty() )
  {
    assert( m_shielding_thickness.empty() );
    m_shielding_att = nullptr;
    return;
  }//if( m_shielding_name.empty() )
  
  if( !db )
    throw runtime_error( "No MaterialDB passed in." );
  
  const Material * const material_raw = db->material(m_shielding_name);
  if( !material_raw )
    throw runtime_error( "No material named '" + m_shielding_name + "' available." );
  
  const auto material = make_shared<Material>(*material_raw);
  const double thickness = PhysicalUnits::stringToDistance( m_shielding_thickness );
  if( thickness < 0.0 )
    throw runtime_error( "Distance '" + m_shielding_thickness + "' is negative." );
  
  m_shielding_att = [material, thickness]( float energy ) -> double {
    const double att_coef = GammaInteractionCalc::transmition_coefficient_material( material.get(), energy, thickness );
    return exp( -1.0 * att_coef );
  };
}//void ReferenceLineInfo::setShieldingAttFcn( const MaterialDB *db )


std::shared_ptr<ReferenceLineInfo> ReferenceLineInfo::generateRefLineInfo( RefLineInput input )
{
  // The gamma or xray energy below which we wont show lines for.
  //  x-rays for nuclides were limited at above 10 keV, so we'll just impose this
  //  as a lower limit to show to be consistent.
  const float lower_photon_energy = 10.0f;
  
  auto answer_ptr = make_shared<ReferenceLineInfo>();
  ReferenceLineInfo &answer = *answer_ptr;
  
  
  // We want to set the final _modified_ version of input to the answer, before
  // returning, so we'll just use a helper for this.
  //on_scope_exit on_exit( [&answer_ptr, &input](){
  //  answer_ptr->m_input = input;
  //  answer_ptr->lineColor = input.m_color;
  //  } );
  answer_ptr->m_input = input;
  
  if( input.m_input_txt.empty() )
  {
    answer.m_validity = ReferenceLineInfo::InputValidity::Blank; //Should already be this value, but being explicit
    
    return answer_ptr;
  }
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  double age = 0.0;
  const SandiaDecay::Nuclide * const nuc = db->nuclide( input.m_input_txt );
  
  if( !nuc )
  {
    input.m_age = "";
  }else
  {
    input.m_input_txt = nuc->symbol;
    input.m_promptLinesOnly = (input.m_promptLinesOnly && nuc->canObtainPromptEquilibrium());
    
    answer.m_nuclide = nuc;
    
    if( input.m_promptLinesOnly )
    {
      age = 5.0*nuc->promptEquilibriumHalfLife();
      input.m_age = PhysicalUnits::printToBestTimeUnits(age, 2);
    }else if( input.m_age == "" )
    {
      age = PeakDef::defaultDecayTime( nuc, &input.m_age );
    }else if( nuc->isStable() )
    {
      age = 0;
      input.m_age = "";
      answer.m_input_warnings.push_back( nuc->symbol + " is a stable isotope." );
      answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
      answer_ptr->m_input = input;
      
      return answer_ptr;
    }else
    {
      try
      {
        age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( input.m_age, nuc->halfLife );
      } catch( std::exception & )
      {
        answer.m_input_warnings.push_back( "Invalid nuclide age input." );
        answer.m_validity = ReferenceLineInfo::InputValidity::InvalidAge;
        
        answer_ptr->m_input = input;
        
        return answer_ptr;
      }//try /catch to get the age
      
      if( age > 100.0 * nuc->halfLife || age < 0.0 )
      {
        const string old_age_str = input.m_age;
        age = PeakDef::defaultDecayTime( nuc, &input.m_age );
        answer.m_input_warnings.push_back( "Changed age to a more reasonable value for "
                                          + nuc->symbol + " from '" + old_age_str + "' to " + input.m_age );
      }
    }//if( prompt only ) / else
  }//if( nuc )
  
  
  const bool check_element = (!nuc && (input.m_input_txt.find_first_of( "0123456789" ) == string::npos));
  const SandiaDecay::Element * const el = check_element ? db->element( input.m_input_txt ) : nullptr;
  
  if( el )
  {
    input.m_input_txt = el->symbol;
    input.m_age = ""; //JIC
    input.m_showXrays = true;
    
    answer.m_element = el;
  }//if( !nuc )
  
  
  string reaction_txt; //CSV list of reactions - I think for our context, only ever a single reaction
  vector<ReactionGamma::ReactionPhotopeak> rctn_gammas;
  if( !nuc && !el )
  {
    const size_t open_paren = input.m_input_txt.find( "(" );
    const size_t close_paren = (open_paren == string::npos) ? string::npos
    : input.m_input_txt.find( ")", open_paren );
    
    if( close_paren != string::npos )
    {
      try
      {
        const ReactionGamma *rctnDb = ReactionGammaServer::database();
        if( rctnDb )
        {
          reaction_txt = rctnDb->gammas( input.m_input_txt, rctn_gammas );
          SpecUtils::ireplace_all( reaction_txt, "'", "" );
          
          // We will fill in answer.m_reactions later on.
          
          // Note: we are-not setting the input to the reaction_txt, as we dont need to bother
          //       the user that the underlying data is for the isotopics of the element, and also
          //       we dont want to change the user input, because we want to keep contributions
          //       normalized to natural abundance.
          //input.m_input_txt = reaction_txt;
          
          input.m_age = "";
          input.m_showGammas = true;
        }
      }catch( std::exception &e )
      {
        // Not a reaction
      }//try / catch
    }//if( (open_paren != string::npos) && (close_paren != string::npos) )
  }//if( !nuc && !el )
  
  
  vector<OtherRefLine> otherRefLinesToShow;
  const bool is_background = (nuc || el || !rctn_gammas.empty()) ? false
  : SpecUtils::icontains( input.m_input_txt, "background" );
  if( is_background )
  {
    input.m_input_txt = "background";
    
    for( const OtherRefLine &bl : BackgroundLines )
    {
      const bool isXray = (std::get<3>(bl) == OtherRefLineType::BackgroundXRay);
      if( (isXray && input.m_showXrays) || (!isXray && input.m_showGammas) )
        otherRefLinesToShow.push_back( bl );
    }//for( const BackgroundLine &bl : BackgroundLines )
    
    input.m_age = "";
  }//if( is_background )
  
  bool is_custom_energy = false;
  if( !nuc && !el && rctn_gammas.empty() && !is_background )
  {
    try
    {
      const float energy = static_cast<float>(PhysicalUnits::stringToEnergy( input.m_input_txt ));
      
      //BackgroundLine<Energy, RelBranchRatio, "Symbol", OtherRefLineType, "Description">
      OtherRefLine line{energy, 1.0f, "", OtherRefLineType::OtherBackground, input.m_input_txt};
      otherRefLinesToShow.push_back( line );
      is_custom_energy = true;
      input.m_age = "";
    }catch( std::exception & )
    {
    }
  }//if( !nuc && !el && rctnGammas.empty() )
  
  
  if( !nuc && !el && rctn_gammas.empty() && !is_background && !is_custom_energy )
  {
    answer.m_validity = ReferenceLineInfo::InputValidity::InvalidSource;
    answer.m_input_warnings.push_back( input.m_input_txt + " is not a valid isotope, element, reaction, or energy." );
    
    answer_ptr->m_input = input;
    
    return answer_ptr;
  }//if( we couldnt match input text to a source )
  
  answer.m_validity = ReferenceLineInfo::InputValidity::Valid;
  
  input.m_showGammas = input.m_showGammas;
  input.m_showXrays = input.m_showXrays;
  input.m_showAlphas = input.m_showAlphas;
  input.m_showBetas  = input.m_showBetas;
  // We will also update showing cascades later based on answer.m_has_coincidences
  input.m_showCascades = input.m_showCascades;
  
  answer_ptr->m_input = input;
  
  
  if( nuc )
    answer.m_source_type = ReferenceLineInfo::SourceType::Nuclide;
  else if( el )
    answer.m_source_type = ReferenceLineInfo::SourceType::FluorescenceXray;
  else if( !rctn_gammas.empty() )
    answer.m_source_type = ReferenceLineInfo::SourceType::Reaction;
  else if( is_background )
    answer.m_source_type = ReferenceLineInfo::SourceType::Background;
  else if( is_custom_energy )
    answer.m_source_type = ReferenceLineInfo::SourceType::CustomEnergy;
  else
    answer.m_source_type = ReferenceLineInfo::SourceType::None;
  
  bool use_particle[SandiaDecay::ProductType::XrayParticle + 1] = {false};
  
  // We'll loop + switch over SandiaDecay::ProductType so we'll at least get
  //  a compiler warning if SandiaDecay::ProductType changes.
  for( auto type = SandiaDecay::ProductType(0);
      type <= SandiaDecay::ProductType::XrayParticle;
      type = SandiaDecay::ProductType(type+1) )
  {
    switch( type )
    {
      case SandiaDecay::ProductType::BetaParticle:
        use_particle[type] = input.m_showBetas;
        break;
        
      case SandiaDecay::ProductType::GammaParticle:
        use_particle[type] = (input.m_showGammas || input.m_showCascades);
        break;
        
      case SandiaDecay::ProductType::AlphaParticle:
        use_particle[type] = input.m_showAlphas;
        break;
        
      case SandiaDecay::ProductType::PositronParticle:
        use_particle[type] = input.m_showGammas;
        break;
        
      case SandiaDecay::ProductType::CaptureElectronParticle:
        break;
        
      case SandiaDecay::ProductType::XrayParticle:
        use_particle[type] = (el || input.m_showXrays);
        break;
    }//switch(type)
  }//for( loop over SandiaDecay::ProductType )
  
  
  vector<ReferenceLineInfo::RefLine> lines;
  
  //transition, first gamma BR, first gamma energy, second gamma energy, coincidence fraction, second gamma BR (just for debug)
  vector<tuple<const SandiaDecay::Transition *, double, float, float, float, double>> gamma_coincidences;
  
  if( nuc )
  {
    SandiaDecay::NuclideMixture mixture;
    
    if( input.m_promptLinesOnly )
    {
      age = 0.0;
      mixture.addNuclideInPromptEquilibrium( nuc, 1.0E-3 * SandiaDecay::curie );
    } else
    {
      mixture.addNuclideByActivity( nuc, 1.0E-3 * SandiaDecay::curie );
    }//if( we want promt only ) / else
    
    
    const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( age );
    
    // x-rays are slightly problematic - we can *almost* treat them like gammas, but for
    //  some decays they essentually get duplicated - so instead we'll be a little ineffient
    //  and track them seperately.
    vector<SandiaDecay::EnergyRatePair> xrays = use_particle[SandiaDecay::ProductType::XrayParticle]
    ? mixture.xrays( age )
    : vector<SandiaDecay::EnergyRatePair>{};
    
    const double parent_activity = nuc ? mixture.activity( age, nuc ) : 0.0;
    
    // We will accumulate positrons as a single line, and just assign the transition
    //  as the first one we run into
    //  TODO: should make sure the first transition we run into is the heaviest one.
    ReferenceLineInfo::RefLine positron_line;
    positron_line.m_energy = 510.9989 * PhysicalUnits::keV;
    positron_line.m_decay_intensity = 0.0;
    positron_line.m_parent_nuclide = nuc;
    positron_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
    positron_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Annihilation;
    
    
    for( const SandiaDecay::NuclideActivityPair &nap : activities )
    {
      const SandiaDecay::Nuclide *nuclide = nap.nuclide;
      const double activity = nap.activity;
      
      for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
      {
        for( const SandiaDecay::RadParticle &particle : transition->products )
        {
          assert( particle.type <= SandiaDecay::ProductType::XrayParticle );
          
          if( !use_particle[particle.type] )
            continue;
          
          if( ((particle.type == SandiaDecay::GammaParticle)
               || (particle.type == SandiaDecay::XrayParticle))
             && (particle.energy < lower_photon_energy) )
          {
            continue;
          }
          
          if( particle.type == SandiaDecay::PositronParticle )
          {
            if( !positron_line.m_transition )
              positron_line.m_transition = transition;
            
            const double br = activity * particle.intensity
            * transition->branchRatio / parent_activity;
            
            positron_line.m_decay_intensity += 2.0 * br;
            
            continue;
          }//if( particle.type == SandiaDecay::PositronParticle )
          
          
          ReferenceLineInfo::RefLine line;
          line.m_parent_nuclide = nuc;
          line.m_energy = particle.energy;
          line.m_transition = transition;
          
          if( particle.type == SandiaDecay::XrayParticle )
          {
            size_t index = 0;
            for( ; index < xrays.size(); ++index )
            {
              if( fabs( xrays[index].energy - particle.energy ) < 1.0E-6 )
                break;
            }
            
            if( index < xrays.size() )
            {
              line.m_decay_intensity = xrays[index].numPerSecond / parent_activity;
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
              line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
              
              lines.push_back( line );
              
              // Erase this x-ray so we dont double-count it
              xrays.erase( begin( xrays ) + index );
            } else
            {
              // We've already accounted for this energy.
            }
            
            continue;
          }//if( particle.type == SandiaDecay::XrayParticle )
          
          
          if( !answer.m_has_coincidences && (particle.type == SandiaDecay::GammaParticle) )
            answer.m_has_coincidences = !particle.coincidences.empty();
          
          const double br = activity * particle.intensity
          * transition->branchRatio / parent_activity;
          
          if( input.m_showCascades && (particle.type == SandiaDecay::GammaParticle) )
          {
            for( size_t coinc_index = 0; coinc_index < particle.coincidences.size(); ++coinc_index )
            {
              const unsigned short int part_ind = particle.coincidences[coinc_index].first;
              const float fraction = particle.coincidences[coinc_index].second;
              assert( part_ind < transition->products.size() );
              if( part_ind < transition->products.size() )
              {
                const SandiaDecay::RadParticle &coinc_part = transition->products[part_ind];
                
                // The BR of second gamma is just for debugging
                const double second_br = activity * coinc_part.intensity
                * transition->branchRatio / parent_activity;
                
                if( coinc_part.type == SandiaDecay::ProductType::GammaParticle )
                  gamma_coincidences.emplace_back( transition, br, particle.energy, coinc_part.energy, fraction, second_br );
              }//if (part_ind < transition->products.size())
            }//for( loop over coincidences )
          }//if( show cascade gammas )
          
          
          // If type is GammaParticle, we could be here if user selected to show cascade, but
          //  not actual gammas, so we need to check for this, and if so not add this gamma
          //  in to be shown
          if( (particle.type == SandiaDecay::GammaParticle) && !input.m_showGammas )
            continue;
          
          line.m_decay_intensity = br;
          line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
          
          switch( particle.type )
          {
            case SandiaDecay::ProductType::BetaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Beta;
              break;
              
            case SandiaDecay::ProductType::GammaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
              break;
              
            case SandiaDecay::ProductType::AlphaParticle:
              line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Alpha;
              break;
              
            case SandiaDecay::ProductType::PositronParticle:
            case SandiaDecay::ProductType::CaptureElectronParticle:
            case SandiaDecay::ProductType::XrayParticle:
              assert( 0 );
              continue;
              break;
          }//switch( particle.type )
          
          lines.push_back( line );
        }//for( const SandiaDecay::RadParticle &particle : transition->products )
      }//for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
    }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
    
    if( positron_line.m_decay_intensity > 0.0 )
      lines.push_back( positron_line );
  }//if( nuc )
  
  // Update showing cascades based on if there are actually any present
  input.m_showCascades = (answer.m_has_coincidences && input.m_showCascades);
  answer.m_input.m_showCascades = (answer.m_has_coincidences && input.m_showCascades);
  
  
  
  if( el )
  {
    for( const SandiaDecay::EnergyIntensityPair &eip : el->xrays )
    {
      if( eip.energy < lower_photon_energy )
        continue;
      
      ReferenceLineInfo::RefLine line;
      line.m_element = el;
      line.m_energy = eip.energy;
      line.m_decay_intensity = eip.intensity;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      
      lines.push_back( line );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( m_showXrays->isChecked() )
  
  
  if( !rctn_gammas.empty() )
  {
    for( const ReactionGamma::ReactionPhotopeak &eip : rctn_gammas )
    {
      if( eip.reaction )
        answer.m_reactions.insert( eip.reaction );
      
      if( eip.energy < lower_photon_energy )
        continue;
      
      ReferenceLineInfo::RefLine line;
      line.m_reaction = eip.reaction;
      line.m_energy = eip.energy;
      line.m_decay_intensity = eip.abundance;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      
      lines.push_back( line );
    }//for( const SandiaDecay::EnergyIntensityPair &eip : element->xrays )
  }//if( !rctn_gammas.empty() )
  
  
  if( !otherRefLinesToShow.empty() )
  {
    assert( is_background || is_custom_energy );
    assert( !nuc && !el && rctn_gammas.empty() );
    
    const SandiaDecay::Nuclide *u238 = is_background ? db->nuclide( "U238" ) : nullptr;
    const SandiaDecay::Nuclide *u235 = is_background ? db->nuclide( "U235" ) : nullptr;
    const SandiaDecay::Nuclide *th232 = is_background ? db->nuclide( "Th232" ) : nullptr;
    const SandiaDecay::Nuclide *ra226 = is_background ? db->nuclide( "Ra226" ) : nullptr;
    const SandiaDecay::Nuclide *k40 = is_background ? db->nuclide( "K40" ) : nullptr;
    
    for( const OtherRefLine &bl : otherRefLinesToShow )
    {
      if( std::get<0>( bl ) < lower_photon_energy )
        continue;
      
      ReferenceLineInfo::RefLine line;
      line.m_energy = std::get<0>( bl );
      line.m_decay_intensity = std::get<1>( bl );
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      
      if( !std::get<2>(bl).empty() )
        line.m_decaystr = std::get<2>(bl) + ", ";
      
      switch( std::get<3>( bl ) )
      {
        case OtherRefLineType::U238Series:
          line.m_parent_nuclide = u238;
          line.m_decaystr += "U238 series";
          break;
          
        case OtherRefLineType::U235Series:
          line.m_parent_nuclide = u235;
          line.m_decaystr += "U235 series";
          break;
          
        case OtherRefLineType::Th232Series:
          line.m_parent_nuclide = th232;
          line.m_decaystr += "Th232 series";
          break;
          
        case OtherRefLineType::Ra226Series:
          line.m_parent_nuclide = ra226;
          line.m_decaystr += "U238 (Ra226) series";
          break;
          
        case OtherRefLineType::K40Background:
          line.m_parent_nuclide = k40;
          line.m_decaystr += "Primordial";
          break;
          
        case OtherRefLineType::BackgroundXRay:
        {
          //std::get<2>( bl ) will be like "Pb xray"
          line.m_decaystr = std::get<2>( bl );
          if( !std::get<4>( bl ).empty() )
            line.m_decaystr += (line.m_decaystr.empty() ? "" : ", ") + std::get<4>( bl );
          line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Xray;
          vector<string> parts;
          SpecUtils::split( parts, std::get<2>( bl ), " " );
          if( !parts.empty() )
            line.m_element = db->element( parts[0] );
          break;
        }//case OtherRefLineType::BackgroundXRay:
          
        case OtherRefLineType::OtherBackground:
        case OtherRefLineType::BackgroundReaction:
        {
          line.m_decaystr = std::get<2>( bl );
          
          line.m_parent_nuclide = db->nuclide( std::get<2>( bl ) );
          // TODO: try to get reaction if didnt get nuclide - also, nuclide list may be CSV, could split that
          //const ReactionGamma *rctnDb = ReactionGammaServer::database();
          //if( rctnDb )
          //{
          //  rctnDb->gammas( input.m_input_txt, rctn_gammas );
          //...
          //}
          
          if( !std::get<4>( bl ).empty() )
            line.m_decaystr += (line.m_decaystr.empty() ? "" : ", ") + std::get<4>( bl );
          
          // TODO: We should probably do a little better than string matching to define as a single/double escape peak.
          if( line.m_decaystr.find("S.E.") != string::npos )
            line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::SingleEscape;
          else if( line.m_decaystr.find("D.E.") != string::npos )
            line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape;
          
          break;
        }
      }//switch( get<3>(*bl) )
      
      
      lines.push_back( line );
    }//for( otherRefLinesToShow )
  }//if( !otherRefLinesToShow.empty() )
  
  
  // Now calc detector response and shielding
  //  Up to now, we shouldnt have any escape or sum gammas in answer.m_ref_lines
  double max_alpha_br = 0.0, max_beta_br = 0.0, max_photon_br = 0.0;
  for( ReferenceLineInfo::RefLine &line : lines )
  {
    assert( (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Normal)
           || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::Annihilation)
           || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::SingleEscape)
           || (line.m_source_type == ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape) );
    
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
        max_alpha_br = std::max( max_alpha_br, line.m_decay_intensity );
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Beta:
        max_beta_br = std::max( max_beta_br, line.m_decay_intensity );
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
      {
        double energy = line.m_energy;
        switch( line.m_source_type )
        {
          case ReferenceLineInfo::RefLine::RefGammaType::Normal:
          case ReferenceLineInfo::RefLine::RefGammaType::Annihilation:
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::SingleEscape:
            // TODO: Need to put in S.E. DRF factor here
            energy += 510.998950;
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::DoubleEscape:
            // TODO: Need to put in D.E. DRF factor here
            energy += 2.0*510.998950;
            break;
            
          case ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak:
          case ReferenceLineInfo::RefLine::RefGammaType::SumGammaPeak:
            assert( 0 );
            break;
        }//switch( line.m_source_type )
        
        if( input.m_det_intrinsic_eff )
          line.m_drf_factor = input.m_det_intrinsic_eff( energy );
        
        if( input.m_shielding_att )
          line.m_shield_atten = input.m_shielding_att( energy );
        
        max_photon_br = std::max( max_photon_br,
                                 line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten );
        break;
      }
    }//switch( line.m_particle_type )
  }//for( ReferenceLineInfo::RefLine &line : lines )
  
  
  const double alpha_sf = ((max_alpha_br > 0.0) && !IsNan(max_alpha_br)) ? (1.0 / max_alpha_br) : 1.0;
  const double beta_sf = ((max_beta_br > 0.0) && !IsNan(max_beta_br)) ? (1.0 / max_beta_br) : 1.0;
  const double photon_sf = ((max_photon_br > 0.0) && !IsNan(max_photon_br)) ? (1.0 / max_photon_br) : 1.0;
  
  for( ReferenceLineInfo::RefLine &line : lines )
  {
    switch( line.m_particle_type )
    {
      case ReferenceLineInfo::RefLine::Particle::Alpha:
        line.m_particle_sf_applied = alpha_sf;
        line.m_normalized_intensity = line.m_decay_intensity * alpha_sf;
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Beta:
        line.m_particle_sf_applied = beta_sf;
        line.m_normalized_intensity = line.m_decay_intensity * beta_sf;
        break;
        
      case ReferenceLineInfo::RefLine::Particle::Gamma:
      case ReferenceLineInfo::RefLine::Particle::Xray:
        line.m_particle_sf_applied = photon_sf;
        line.m_normalized_intensity = photon_sf * line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten;
        break;
    }//switch( line.m_particle_type )
    
    // We wont filter out lines smaller than wanted here
    if( (line.m_transition || is_background)
       && (line.m_normalized_intensity <= input.m_lower_br_cutt_off
           || IsInf( line.m_normalized_intensity )
           || IsNan( line.m_normalized_intensity )) )
    {
      continue;
    }
    
    
    // Now lets fill out line.m_decaystr
    if( line.m_decaystr.empty() )
    {
      if( line.m_transition )
      {
        if( line.m_transition->parent )
          line.m_decaystr = line.m_transition->parent->symbol;
        if( line.m_transition->child )
          line.m_decaystr += " to " + line.m_transition->child->symbol;
        
        // TODO: for alphas and betas its pretty redundant to have this next line (I guess its redundant no matter what actually)
        line.m_decaystr += string( " via " ) + SandiaDecay::to_str( line.m_transition->mode );
        
        switch( line.m_transition->mode )
        {
          case SandiaDecay::AlphaDecay:
          case SandiaDecay::BetaDecay:
          case SandiaDecay::BetaPlusDecay:
          case SandiaDecay::ProtonDecay:
            line.m_decaystr += " decay";
            break;
            
          case SandiaDecay::IsometricTransitionDecay:
            line.m_decaystr += " transition";
            break;
          
          default:
            break;
        }//switch( line.m_transition->mode )
      }else if( line.m_reaction )
      {
        // TODO: we can probably come up with a better way to describe reactions
        if( line.m_reaction->targetNuclide && line.m_reaction->targetElement )
        {
          switch( line.m_reaction->type )
          {
            case AlphaNeutron:   line.m_decaystr = "Alphas on "; break;
            case NeutronAlpha:   line.m_decaystr = "Neutrons on "; break;
            case AlphaProton:    line.m_decaystr = "Alphas on "; break;
            case NeutronCapture: line.m_decaystr = "Neutron capture by "; break;
            case NeutronInelasticScatter: line.m_decaystr = "Neutron inelastic scatter on "; break;
            case AnnihilationReaction: break;
            case NumReactionType:      break;
          }//switch( line.m_reaction->type )
          
          line.m_decaystr += line.m_reaction->targetNuclide->symbol;
          //if( line.m_reaction->productNuclide )
          //  line.m_decaystr += " to give " + line.m_reaction->productNuclide->symbol;
        }
      }else if( line.m_element )
      {
        line.m_decaystr = line.m_element->name + " fluorescence";
      }
    }//if( line.m_decaystr.empty() )
    
    const double &amp = line.m_normalized_intensity;
    if( !IsNan( amp ) && !IsInf( amp )
       && (amp >= std::numeric_limits<float>::min()) // numeric_limits<float>::min()==1.17549e-38
       && (amp > input.m_lower_br_cutt_off)
       )
    {
      answer.m_ref_lines.push_back( line );
    }
  }//for( ReferenceLineInfo::RefLine &line : answer.m_ref_lines )
  
  
  // If we add in escape peaks - we could put them in here
  
  // Add in coincident gammas
  if( !gamma_coincidences.empty() )
  {
    double max_coincidence_br = 0.0;
    vector<ReferenceLineInfo::RefLine> coinc_ref_lines;
    for( const auto &casc : gamma_coincidences )
    {
      const SandiaDecay::Transition *const &trans = get<0>( casc );
      const double &first_br = get<1>( casc );
      const float &first_energy = get<2>( casc );
      const float &second_energy = get<3>( casc );
      const float &coinc_frac = get<4>( casc );
      const double &second_br = get<5>( casc );
      
      const float energy = first_energy + second_energy;
      
      ReferenceLineInfo::RefLine line;
      line.m_energy = energy;
      line.m_decay_intensity = first_br * coinc_frac;
      line.m_parent_nuclide = nuc;
      line.m_transition = trans;
      line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::CoincidenceSumPeak;
      
      if( input.m_det_intrinsic_eff )
        line.m_drf_factor = input.m_det_intrinsic_eff( first_energy ) * input.m_det_intrinsic_eff( second_energy );
      
      if( input.m_shielding_att )
        line.m_shield_atten = input.m_shielding_att( first_energy ) * input.m_shielding_att( second_energy );
      
      const double amp = line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten;
      assert( !IsNan( amp ) && !IsInf( amp ) );
      if( IsNan( amp ) || IsInf( amp ) )
      {
        cerr << "Unexpected NaN or Inf coincidence amp." << endl;
        continue;
      }
      
      line.m_decaystr = "Cascade sum";
      if( trans && trans->parent )
        line.m_decaystr += " " + trans->parent->symbol;
      if( trans && trans->child )
        line.m_decaystr += " to " + trans->child->symbol;
      
      char buffer[128];
      snprintf( buffer, sizeof( buffer ),
               " (%.1f + %.1f keV, coinc=%.3g)",
               first_energy, second_energy, coinc_frac );
      
      line.m_decaystr += buffer;
      
      coinc_ref_lines.push_back( std::move( line ) );
      
      max_coincidence_br = std::max( max_coincidence_br, amp );
    }//for( loop over cascades )
    
    assert( coinc_ref_lines.empty()
           || ((max_coincidence_br > 0.0) && !IsNan( max_coincidence_br )) );
    
    // Scale the coincidence line amplitudes to be between 0
    for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    {
      const double sf = 1.0 / max_coincidence_br;
      line.m_particle_sf_applied = sf;
      const double amp = line.m_decay_intensity * line.m_drf_factor * line.m_shield_atten * sf;
      line.m_normalized_intensity = amp;
    }//for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    
    // There can be tons of cascade sums (4834 for U238), we'll limit the number
    //   we draw to an arbitrary 350, because this is even more than I expect to
    //   be relevant (although I didnt actually check this).
    //  TODO: limit based on importance, and not a flat limit, e.g., use something like
    //        yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))
    const size_t max_cascade_sums = 350;
    if( coinc_ref_lines.size() > max_cascade_sums )
    {
      std::sort( begin( coinc_ref_lines ), end( coinc_ref_lines ),
                []( const ReferenceLineInfo::RefLine &lhs, const ReferenceLineInfo::RefLine &rhs ) -> bool {
        if( lhs.m_normalized_intensity == rhs.m_normalized_intensity )
          return lhs.m_energy > rhs.m_energy;
        return lhs.m_normalized_intensity > rhs.m_normalized_intensity;
      } );
      
      cout << "Resizing cascade sums from " << coinc_ref_lines.size() << " to " << max_cascade_sums << endl;
      coinc_ref_lines.resize( max_cascade_sums );
    }//if( coinc_ref_lines.size() > 350 )
    
    answer.m_ref_lines.reserve( answer.m_ref_lines.size() + coinc_ref_lines.size() );
    
    for( const ReferenceLineInfo::RefLine &line : coinc_ref_lines )
    {
      const double &amp = line.m_normalized_intensity;
      if( !IsNan( amp ) && !IsInf( amp )
         && (amp >= std::numeric_limits<float>::min()) // numeric_limits<float>::min()==1.17549e-38
         && (amp > input.m_lower_br_cutt_off)
         )
        answer.m_ref_lines.push_back( line );
    }//for( ReferenceLineInfo::RefLine &line : coinc_ref_lines )
  }//if( !gamma_coincidences.empty() )
  
  
  //Clientside javascript currently doesnt know about this garuntee that gamma
  //  lines will be sorted by energy.
  answer.sortByEnergy();
  
  return answer_ptr;
}//std::shared_ptr<ReferenceLineInfo> generateRefLineInfo()
