#ifndef ShieldingSourceFitCalc_h
#define ShieldingSourceFitCalc_h
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
#include <memory>
#include <vector>


#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
#include <boost/optional.hpp>
#endif


//Forward declarations
struct Material;

class PeakDef;
class MaterialDB;
class PopupDivMenu;
class DetectorPeakResponse;

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}//namespace SandiaDecay

namespace ROOT
{
  namespace Minuit2
  {
    class MnUserParameters;
  }//namespace Minuit2
}//namespace ROOT

namespace GammaInteractionCalc
{
  enum class GeometryType : int;
  enum class TraceActivityType : int;
};


/** This namespace is structs that represent the data users input in the `ShieldingSelect` class, and the inputs of
 ShieldingSourceDisplay.h/.cpp
 */
namespace ShieldingSourceFitCalc
{
  /** Struct holding information cooresponding to the `TraceSrcDisplay` class defined in ShieldingSelect.cpp;
   represents information about a trace-source in a shielding (e.g., a volumetric source distributed unifrmly in a shielding
   material, but does not effect the attenuation or density of that material.
   */
  struct TraceSourceInfo
  {
    GammaInteractionCalc::TraceActivityType m_type;
    bool m_fitActivity;
    const SandiaDecay::Nuclide *m_nuclide;
    double m_activity; //units are according to #m_type
    float m_relaxationDistance; //only applicable to #TraceActivityType::ExponentialDistribution
    
    TraceSourceInfo();
    void serialize( rapidxml::xml_node<char> *parent_node ) const;
    void deSerialize( const rapidxml::xml_node<char> *shield_node );
    
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    static void equalEnough( const TraceSourceInfo &lhs, const TraceSourceInfo &rhs );
#endif
  };//struct TraceSourceInfo

  
  /** Holds the information cooresponding to the `ShieldingSelect` widget. */
  struct ShieldingInfo
  {
    /** The geometry of the shielding.
     If a generic material, then must be #GammaInteractionCalc::GeometryType::NumGeometryType.
     */
    GammaInteractionCalc::GeometryType m_geometry;
    
    /** Wether is a generic shieelding (i.e., no physical extent, but specified by atomic number and areal density),
     or a physical shielding, that should have a valid material.
     */
    bool m_isGenericMaterial;
    
    /** A kinda vestigial variable that indicates if the user input that made this info was Shielding/Source fit tool,
     or somewhere else where fitting the material dimesnions or AN/AD wasnt intended.
     */
    bool m_forFitting;
    
    /** Material shielding is made out of; will be nullptr if generic material. */
    std::shared_ptr<const Material> m_material;
    
    /** Dimesnisons of this shielding; the meaning of the entries differs depending on the geometry,
     or if a generic material.
     
     Spherical: ['Thickness', n/a, n/a]
     Cylinder:  ['Radius','Length',n/a]
     Rectangle: ['Width','Height','Depth']
     Generic:   ['AtomicNumber','ArealDensity',n/a]
     */
    double m_dimensions[3];
    bool m_fitDimensions[3];
    
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    boost::optional<double> m_truthDimensions[3];
    boost::optional<double> m_truthDimensionsTolerances[3];
#endif
    
    // Self-atten source stuff
    bool m_fitMassFrac;
    
    /** Nuclide mass-fractions are only fit within the same element, and they are constrained to be the sum fraction
     of all the nuclides, for that element, that are being fit.
     
     That is, if you have {{I131, 0.1}, {I124,0.01},{Cs131,0.01},{Cs137,0.02}} (and 0.89 stable I127, and 0.97 Cs133), then
     if you fit mass fractions, I131+I124 will always sum to be 0.11 fraction of the Iodine, and Cs137+Cs131 will always
     sum to be 0.03 of the Cesium.
     */
    std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractions;
    
    // Trace-source stuff
    std::vector<TraceSourceInfo> m_traceSources;
    
    
    ShieldingInfo();
    
    rapidxml::xml_node<char> *serialize( rapidxml::xml_node<char> *parent_node ) const;
    void deSerialize( const rapidxml::xml_node<char> *shield_node, MaterialDB *materialDb );
    
    /** Encodes current tool state to app-url format.  Returned string is just the query portion of URL;
     so will look something like "V=1&G=S&D1=1.2cm", and it will not be url-encoded.
     
     TODO: Currently does not encode self-attenuating source, trace-source information, or "truth" values.
     */
    std::string encodeStateToUrl() const;
    
    void handleAppUrl( std::string query_str, MaterialDB *materialDb );
    
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    static void equalEnough( const ShieldingInfo &lhs, const ShieldingInfo &rhs );
#endif
    
    static const int sm_xmlSerializationMajorVersion;
    static const int sm_xmlSerializationMinorVersion;
  };//struct ShieldingInfo
  
}//namespace ShieldingSourceFitCalc

#endif //ShieldingSourceFitCalc_h
