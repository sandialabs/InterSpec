#ifndef ReferenceLineInfo_h
#define ReferenceLineInfo_h
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

#include <tuple>
#include <vector>
#include <string>

#include <Wt/WColor>

#include "InterSpec/ReactionGamma.h"

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml


enum class OtherRefLineType
{
  U238Series,    //"U238 series"
  U235Series,    //"U235 series"
  Th232Series,   //"Th232 series"
  Ra226Series,   //"U238 (Ra226) series"
  K40Background, //"Primordial"
  BackgroundXRay,
  BackgroundReaction,
  OtherBackground
};//enum OtherRefLineType

/* Returns simple string representation of OtherRefLineType enum. */
const char *to_str( const OtherRefLineType type );

/* Returns OtherRefLineType from string; input must exactly match what is 
  return by #to_str(OtherRefLineType), or else will return 
  #OtherRefLineType::OtherBackground.
*/
OtherRefLineType other_ref_line_type_from_str( const std::string &str );

/** Used to describe a custom-source reference line.

  Currently used to represent background lines (e.g., user enters "background"), and 
  if the user enters a single energy (e.x., "185 keV").  May be expanded to other use
  cases in the future.

  <Energy, RelBranchRatio, "Symbol", OtherRefLineType, "Description">

  \sa otherRefLines  
*/
typedef std::tuple<float,float,std::string, OtherRefLineType,std::string> OtherRefLine;

//BackgroundLines and BackgroundReactionLines need to be renamed; however, I'm
//  expecting to switch to a different mechanism to suppliment the nuclide and
//  reaction databases, so I'll leave these alone for now.
extern const OtherRefLine BackgroundLines[89];
extern const OtherRefLine BackgroundReactionLines[28];

struct ReferenceLineInfo
{
  //A simple structure to keep in server (c++) memmorry what is being
  //  displayed on the client side, and to allow serializing to the database
  //  what is being displayed, and more-or-less the current state of the
  //  widget (when there is a current line).
  //  Note: some of the variables are probably vestigial by now, and could
  //        stand to be removed.
  
  //nuclide displayed:
  const SandiaDecay::Nuclide *nuclide;
  
  //If (the xrays for) an element is displayed, and no nuclides, then element
  //  will be non-null
  const SandiaDecay::Element *element;
  
  //If reactionGammas is defined, then energies/intensities will be
  //  redundant
  std::vector<ReactionGamma::ReactionPhotopeak> reactionGammas;
  
  /* either background, or custom energy, reference lines displayed.
  
  Used by PeakSearchGuiUtils when setting nuc ID of peaks, when showing 
  background or other ref lines.
  */
  std::vector<OtherRefLine> otherRefLines;
  
  //TODO: place energies, intensities, particlestrs, decaystrs, and
  //      elementstrs into a tuple; also include pointers to Element, Nuclide, 
  //      or Reaction that the peak should be assigned to.  And also unify
  //      this struct to be shared by DecayParticleModel::RowData.  And once
  //      this happens we can get rid of #reactionGammas and #otherRefLines.
  // 
  //energy and intesities of displayed lines, should always be same size.
  //intensities is normalied to be between 0 and 1 for each particle type
  //  (e.g. gammas get a different SF than alphas, etc.  Xrays and gammas share
  //  a scale factor).
  std::vector<double> energies, intensities;
  
  //Description of the displayed lines, should always be same size as
  //  energies and intensities vectors
  std::vector<std::string> particlestrs, decaystrs, elementstrs;
  
  //Scale factors applied to intensities vector, for each particle type.
  //To extract the original BR for a particle typ, you could do:
  //  double br = this->intensities[i] * this->particle_sf[this->particlestrs[i]];
  std::map<std::string,double> particle_sf;

  
  //The following variables are necassarry to help serialize the state of the
  //  widget for loading later
  bool showGammas, showXrays, showAlphas, showBetas, showCascades;
  bool showLines, promptLinesOnly, isOtherRef, isReaction, displayLines;
  double age, lowerBrCuttoff;
  
  //labelTxt: what the user enetered to get this set of gammas
  std::string labelTxt;
  
  //reactionsTxt: a CSV list of reactions giving rise to reactionGammas
  std::string reactionsTxt;
  
  //shieldingName: name of shieding applied to amplitudes of lines; empty
  //  if none.  Will have format "AN=23, AD=53.1 g/cm2" if generic shielding.
  std::string shieldingName;
  
  //Thickness of shielding.  Will be zero if no shielding or generic shielding
  double shieldingThickness;
  
  //detectorName: Name of detector the amplitude of lines has been modulated
  //  with.
  std::string detectorName;
  
  /** Color to draw the line with.  Alpha not currently supported. */
  Wt::WColor lineColor;
  
  ReferenceLineInfo();
  bool operator==( const ReferenceLineInfo &rhs ) const;
  void reset();
  bool empty() const;
  
  void serialize( rapidxml::xml_node<char> *parent_node ) const;
  void deSerialize( const rapidxml::xml_node<char> *node );
  
  //toJson(...): appends the JS object necessary to draw this to 'json'.
  //  number_lines cooresponds to the which color series to draw.
  //  TODO: look into using the Wt::Json library to form the JSON instead of
  //        doing it by hand
  void toJson( std::string &json ) const;
  
  std::string parentLabel() const;
  
  void sortByEnergy();
  
  static const int sm_xmlSerializationVersion;
};//struct ReferenceLineInfo


#endif
