#ifndef ReferenceLinePredef_h
#define ReferenceLinePredef_h
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

#include <map>
#include <string>
#include <vector>
#include <optional>

#include <Wt/WColor>

#include "InterSpec/ReactionGamma.h"

namespace SandiaDecay
{
  struct Element;
  struct Nuclide;
  struct Transition;
}

namespace ReferenceLinePredef
{
  /** Component of a nuclide mixture with age offset and relative activity. */
  struct NucMixComp
  {
    /** An age offset, relative to the rest of the nuclides (i.e., t=0).
     
     A negative value will increase the age of the nuclide, and positive value reduce it.
     
     If a positive value is used, and it is larger than the requested mixture age in the
     "Reference Photopeak" tool, than an age of 0 will be used.
     
     This value corresponds to the "age-offset" attribute of the <Nuc /> element.
     
     TODO: as of 20240217, the use of "age-offset" is not well tested
     */
    double m_age_offset;
    
    /** The relative activity for this nuclide - the total activity fractions dont need to add up to 1, they will be normalized. */
    double m_rel_act;
    
    /** The nuclide - must not be nullptr. */
    const SandiaDecay::Nuclide *m_nuclide;
    
    /** Optional color specified in `add_ref_line.xml` for this nuclide; will override user selection if specified. */
    Wt::WColor m_color;
  };//struct NucMixComp


  /** A mixture of nuclides with specified relative activities and optional age settings. */
  struct NucMix
  {
    std::string m_name;
    double m_default_age;
    std::string m_default_age_str;
    
    /** If the activity fractions are fixed, at whatever time is specified (e.g., no matter the age, the relative
     activity fractions of the parent nuclides will always be the same).
     
     If the XML provides a "reference-age", then this will be false, and changing the age will change the relative
     ratio of the activities.
     
     TODO: as of 20240217, the use of "reference-age" is not well tested
     */
    bool m_fixed_act_fractions;
    
    /** Only applicable to dynamic reference lines.
     The user supplied bias weight for this nuclide - if multiple nuclides may be choosen between for a particular energy,
     then a larger value of this weight will make it more likely this nuclide is shown. Default is 1.0.
     */
    float m_weight = 1.0f;
    
    std::vector<NucMixComp> m_components;
  };//struct NucMix


  /** A custom reference line with energy, branch ratio, and optional nuclide association. */
  struct CustomLine
  {
    float m_energy;
    float m_branch_ratio;
    std::string m_info;
    Wt::WColor m_color;
    
    /** The nuclide the user may have specified.
     If user specifies the nuclide, than the transition will be picked out, assuming a weighting window of 0.25 keV (i.e., if nuclide
     is valid, so will transition), and that gamma was intended.
     */
    const SandiaDecay::Nuclide *m_nuclide;
    const SandiaDecay::Transition *m_transition;
    
    /** If detector efficiency and/or shielding apply to this line, for line height display purposes only. */
    bool m_atten_applies;
    
    CustomLine( float ener, float br, std::string &&inf,
               const SandiaDecay::Nuclide *nuc, const SandiaDecay::Transition *trans,
               const bool atten, Wt::WColor &&color )
    : m_energy(ener),
    m_branch_ratio(br),
    m_info( std::move(inf) ),
    m_color( std::move(color) ),
    m_nuclide( nuc ),
    m_transition( trans ),
    m_atten_applies( atten )
    {
    }
  };//struct CustomLine


  /** A collection of custom source lines loaded from XML configuration. */
  struct CustomSrcLines
  {
    std::string m_name;
    /** Only applicable to dynamic reference lines.
     The user supplied bias weight for this nuclide - if multiple nuclides may be choosen between for a particular energy,
     then a larger value of this weight will make it more likely this nuclide is shown. Default is 1.0.
     */
    float m_weight;
    float m_max_branch_ratio;
    std::vector<CustomLine> m_lines;
  };//struct CustomSrcLines

  /** An individual source specified in the XML file.
   This is only applicable to dynamic reference lines.
   */
  struct IndividualSource
  {
    std::string m_name;
    
    /** Only applicable to dynamic reference lines.
     The user supplied bias weight for this nuclide - if multiple nuclides may be choosen between for a particular energy,
     then a larger value of this weight will make it more likely this nuclide is shown. Default is 1.0.
     */
    float m_weight = 1.0f;
    
    /** Exactly one of the three will be non-null. */
    const SandiaDecay::Nuclide *m_nuclide = nullptr;
    const SandiaDecay::Element *m_element = nullptr;
    const ReactionGamma::Reaction *m_reaction = nullptr;
    
    bool is_background = false;
    std::optional<std::string> shielding_material;
    std::optional<double> shielding_thickness;
    
    std::optional<double> m_age;
    std::optional<Wt::WColor> m_color;
  };//struct IndividualSource

  /** Loads reference line information from an XML file.

   @param filepath Path to the XML file to load
   @param nuc_mixes Map to store loaded nuclide mixtures
   @param custom_lines Map to store loaded custom source lines
   @param individual_sources Optional vector of individual sources in the file; these are used for RefLineDynamic;
          cooresponds to <IndividualSource> XML elements.
   */
  void load_ref_line_file( const std::string& filepath, 
                          std::map<std::string,NucMix>& nuc_mixes, 
                          std::map<std::string,CustomSrcLines>& custom_lines,
                          std::vector<IndividualSource> *individual_sources );

} //namespace ReferenceLinePredef

#endif //ReferenceLinePredef_h
