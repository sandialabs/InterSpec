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

#include <set>
#include <tuple>
#include <vector>
#include <string>

#include <Wt/WColor>

#include "InterSpec/ReactionGamma.h"

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
  struct Transition;
}

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml

class MaterialDB;
class DetectorPeakResponse;

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


/** A struct that contains all the information needed to generate reference photopeak
 lines (i.e. the #ReferenceLineInfo information) from.
 
 This information is filled out either by the ReferencePhotopeakDisplay GUI state, or
 from XML.
 */
struct RefLineInput
{
  RefLineInput();
  void reset();
  
  bool operator==(const RefLineInput &rhs) const;
  
  /** The text that defines the source, ex., "U238", "Pb", "background", "511 keV".
   */
  std::string m_input_txt;
  
  /** The age of the nuclide.
   
   If it is a nuclide source, and the age string is empty, then the default age will
   be used.
   If it is a fluorescence x-ray, reaction, background, or custom energy, then this
   string should be blank.
   
   If non-empty, this string may look like "1.2 y" (1.2 years), "5 HL" (5 half-lives), etc.
   */
  std::string m_age;  //
  
  /** Color to draw the line with.  Alpha not currently supported. */
  Wt::WColor m_color;
  
  double m_lower_br_cutt_off;
  bool m_promptLinesOnly;
  
  bool m_showGammas;
  bool m_showXrays;
  bool m_showAlphas;
  bool m_showBetas;
  bool m_showCascades;
  bool m_showEscapes;
  
  
  /** The name of the detector to display as having had the reference
   lines modulated with.
   */
  std::string m_detector_name;
  
  /** The intrinsic efficiency of the detector, as a function of energy.
   
   Maybe be nullptr.
   
   I'm not hugely crazy about having this defined here, because this
   function will not be serialized to/from XML, which then leads to
   a potentially inconsistent state.
   */
  std::function<float(float)> m_det_intrinsic_eff;
  
  // TODO: I'm not super happy with shielding definition; I'm torn how to represent it everywhere...
  
  /** Shielding name, as can be retrieved from #MaterialDB::material(string).
   Will be blank string if no shielding or a generic shielding.
   
   If this is non-empty, then #m_shielding_an and #m_shielding_ad must be empty.
   */
  std::string m_shielding_name;
  
  /** Thickness of shielding.
   Will be empty if no shielding or generic shielding.
   If non empty, will be a valid distance (like "1.2 cm") that can be interpreted  by
   #PhysicalUnits::
   
   If this is non-empty, then #m_shielding_an and #m_shielding_ad must be empty.
   
   This is a string instead of a float or double to keep user input exact.
   */
  std::string m_shielding_thickness;
  
  /** The atomic number of the generic shielding.
   Will be empty if no shielding or a material shielding.
   
   If this is non-empty, then #m_shielding_name and #m_shielding_thickness will must be empty.
   
   This is a string instead of a float or double to keep user input exact.
   */
  std::string m_shielding_an;
  
  /** The areal density, in units of g/cm2, of the generic shielding.
   This string holds just a number (e.x, "7.2" or "1"), and does not contain the "g/cm2" units.
   
   Will be empty if no shielding or a material shielding.
   
   If this is non-empty, then #m_shielding_name and #m_shielding_thickness will must be empty.
   
   This is a string instead of a float or double to keep user input exact.
   */
  std::string m_shielding_ad;
  
  /** Function to return the shielding attenuation factor, as a function of energy.
   
   I.e., returns values between 0 (all gammas stopped), and 1.0 (no attenuation).
   
   May be nullptr.
   */
  std::function<double( float )> m_shielding_att;
  
  /** Serializes member variables to xml, appending a <RefLineInput> element onto the
   passed in node.
   
   Note: #m_det_intrinsic_eff and #m_shielding_att do not get serialized.
   */
  void serialize( rapidxml::xml_node<char> *parent_node ) const;
  
  /** De-serializes member variables from xml.
   
   Note: #m_det_intrinsic_eff and #m_shielding_att do not get de-serialized, so you must
         manually set these variables.  See
   */
  void deSerialize( const rapidxml::xml_node<char> *node );
  
  static const int sm_xmlSerializationVersion;
  
  /** Uses (#m_shielding_name and m_shielding_thickness) OR (#m_shielding_an and #m_shielding_ad)
   to set #m_shielding_att.
   
   Throw std::exception on failure, in which case #m_shielding_att wont be changed. 
   */
  void setShieldingAttFcn( const MaterialDB *db );
};//struct RefLineInput


/** A struct to represent reference gamma/alpha/beta lines displayed on the spectrum.
 
 It also defines the information the information the table in the "Reference Photopeaks"
 tab displays; currently the #DecayParticleModel class uses its own #DecayParticleModel::RowData
 struct to hold its information, but eventually it will use this same class.
 */
struct ReferenceLineInfo
{
  /** Generates the reference lines from the given input.
   
   The returned answer will always be non-null, and the function shouldnt
   throw exception, however, you should check #ReferenceLineInfo::m_validity
   to determine success status of the call.
   */
  static std::shared_ptr<ReferenceLineInfo> generateRefLineInfo( RefLineInput input );
  
  /** The additional nuclide mixtures and one-off sources defined in `data/add_ref_line.xml` */
  static std::vector<std::string> additional_ref_line_sources();
  
  /** Loads nuclide mixtures  defined in `data/add_ref_line.xml` into memory.
   
   You do not need to call this function, but you can if you want to avoid the ~2 ms to load it
   the first time you create a reference photopeak widget.
   */
  static void load_nuclide_mixtures();
  
  /** A struct to represent the information for a single reference line.  */
  struct RefLine
  {
    RefLine();

    /** Energy, in keV of the particle.
     
     For gammas and x-rays, this is the energy of the photopeak that would be detected.
     I.e., if #RefLine::m_source_type is #RefGammaType::SingleEscape, #RefGammaType::DoubleEscape,
     then you would need to add 511 keV or 1022 keV to this energy to compute the attenuation
     by the shielding.  If RefLine::m_source_type is #RefGammaType::CoincidenceSumPeak, or
     #RefGammaType::SumGammaPeak, then this is the sum energy of the peaks.
     
     For alphas and betas, this is the end-point energy.
     */
    double m_energy;

    /** The relative intensity of this line, after shielding attenuation, DRF,
     scaling for particle type, and overall normalization.
     I.e., between 0 and 1, and what to use for displaying ref line on chart.
     
     Lines less than or equal to this value will be filtered out (i.e., if set
     to zero, then lines with zero amplitude will be filtered out).
    */
    double m_normalized_intensity;

    /** The detectors intrinsic efficiency.
     Should be between in range [0,1], with it being 1.0 if no DRF is being used.
     */
    float m_drf_factor;
    
    /** The efficiency of the gamma/x-ray to make it through the shielding with
     interacting.
     Should be between in range [0,1], with it being 1.0 if no shielding present.
     */
    float m_shield_atten;

    /** The scale factor applied so that the largest amplitude particle of this 
    type will have a #m_normalized_intensity value of 1.0.  Note decay x-rays 
    are grouped together with gammas to get this factor.
    Ex, if largest gamma BR is 0.7, and no DRF or shielding, then this SF 
      would be 1.0/0.7 = 1.428
    */
    float m_particle_sf_applied;

    /** This is a string that describes to the user where this particle was generated.
     Examples:
     - "Cascade sum Pa234 to U234 (131.3 + 456.7 leV, coinc=0.0514) decay"
     - "Th234 to Pa234m via Beta decay"
     */
    std::string m_decaystr;
    
    /** If not default, the color for this line.
     
     Ordinarily, the line color will be given by `ReferenceLineInfo::m_input.m_color`, but
     if this variable is not default, then the line color will be overridden.  
     So far, this is only used for `CustomLine` sources in the `add_ref_line.xml` file.
     */
    Wt::WColor m_color;
    
    /** Returns one of: "cascade-sum", "sum-gamma", "S.E.", "D.E.", "alpha", "beta", "gamma", "xray"
     */
    const std::string &particlestr() const;

    /** The original branching ratio of the process, to this energy. 
    This is what the "Reference Photopeaks" table displays as the intensity.
    */
    double m_decay_intensity;

    /* Particle type enum so we dont have to include SandiaDecay here... */
    enum class Particle : int
    {
      Alpha,
      Beta, 
      Gamma,
      Xray
    };//enum Particle
    
    /** Decay particle type.  I.e., Gamma, Beta, alpha, capture-e, x-ray. */
    Particle m_particle_type;

    /** This is the ultimate parent nuclide to assign the peak to. 
    Will be nullptr if source is not a nuclide. 

    At most one of {m_parent_nuclide, m_element, m_reaction} will be non-null
    */
    const SandiaDecay::Nuclide *m_parent_nuclide;

    /* The transition actually responsible for the gamma.
     Will be nullptr if source is not a nuclide. 
    */
    const SandiaDecay::Transition *m_transition;

    enum class RefGammaType : int 
    {
      Normal,
      Annihilation,
      SingleEscape,
      DoubleEscape,
      CoincidenceSumPeak,
      SumGammaPeak
    };//enum class RefGammaType : int
    
    /** The gamma type (Normal, Annih., S.E., D.E., x-ray).
    
    Used both for nuclides, as well as reactions.
    */
    RefGammaType m_source_type;

    /** If detector efficiency and/or shielding attenuation should be applied to this lines amplitude.
     This is only used for `CustomLine` sources in the `add_ref_line.xml` file.
     */
    bool m_attenuation_applies;
    
    /** Element giving rise to this fluorescent x-ray.
    Will be element peak is assigned to, if non-null. 
    */
    const SandiaDecay::Element *m_element;

    /** The reaction giving rise to this gamma.
    
    Will be reaction peak is assigned to, if non-null. 
    */
    const ReactionGamma::Reaction *m_reaction;
  };//struct RefLine

  /** Enum to provide status of the validity of the #RefLineInput in generating reference lines */
  enum class InputValidity : int
  {
    Blank,
    InvalidSource,
    InvalidAge,
    Valid
  };//enum class InputValidity : int
  
  /** Enum to specify what type of input was detected for creating reference lines from.
   A few places in the code treat things a little differently, depending on the type of
   reference lines the user would like displayed.
   */
  enum class SourceType : int
  {
    Nuclide,
    FluorescenceXray,
    Reaction,
    Background,
    CustomEnergy,
    NuclideMixture,
    OneOffSrcLines,
    FissionRefLines,
    None
  };//enum class SourceType : int

  /** The status of generating reference lines from the given input. */
  InputValidity m_validity;
  
  /** The input used to generate the reference lines.
   Note: this may not be just a copy of the input given to the function that generated the
   reference lines (i.e. #generateRefLineInfo), but fields may have gotten modified.  For example,
   #RefLineInput::m_input_txt would get changed from "U-238" to "U238", and if #RefLineInput::m_age
   was blank, then it would be set to "20 y".
   That is, the GUI state should be updated based on values in this struct that were changed from
   the input.
   */
  RefLineInput m_input;
  
  /** The reference lines to use for display. */
  std::vector<RefLine> m_ref_lines;
  
  /** If the source had cascade gamma coincidences in it, even if there were not selected to be
   displayed.
   This variable is needed to decide if we should show this cascade sum checkbox in the GUI.
   */
  bool m_has_coincidences;
  
  /** Warnings to give the user about the #RefLineInput.  For example, if the age given was
   totally inappropriate and changed, a warning will be included here.
   */
  std::vector<std::string> m_input_warnings;
  
  /** What type of source was requested. The input source definition is essentially just a string,
   so you can use this field to determine if the input is a nuclide, fluorescence x-ray, reaction,
   etc.
   */
  SourceType m_source_type;

  /** The nuclide the source corresponds to.
   Non-nullptr only if m_source_type==SourceType::Nuclide.
   */
  const SandiaDecay::Nuclide *m_nuclide;
  
  /** The element the source fluorescence x-rays correspond to.
   Non-nullptr only if m_source_type==SourceType::FluorescenceXray.
   */
  const SandiaDecay::Element *m_element;
  
  /** The reaction(s) the source corresponds to.
   There may be multiple reactions present for an input if, for example, the user
   input an element, but the underlying data is available for isotopes of the element,
   then the reference lines will normalized to natural abundance of the element.
   
   non-empty only if m_source_type==SourceType::Reaction.
   */
  std::set<const ReactionGamma::Reaction *> m_reactions;
  
  
  ReferenceLineInfo();
  void reset();
  
  bool operator==(const ReferenceLineInfo &rhs) const;
  
  //toJson(...): appends the JS object necessary to draw this to 'json'.
  //  number_lines cooresponds to the which color series to draw.
  //  TODO: look into using the Wt::Json library to form the JSON instead of
  //        doing it by hand
  void toJson( std::string &json ) const;
  
  
  void sortByEnergy();
};//struct ReferenceLineInfo


#endif
