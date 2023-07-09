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

#include <memory>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"


#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"


using namespace std;

namespace ShieldingSourceFitCalc
{
  const int ShieldingInfo::sm_xmlSerializationMajorVersion = 0;

  /** Change log:
   - 20230705: no version changed, but converted from being member of ShieldingSelect to ShieldingInfo
   */
  const int ShieldingInfo::sm_xmlSerializationMinorVersion = 1;

  
  TraceSourceInfo::TraceSourceInfo()
   : m_type( GammaInteractionCalc::TraceActivityType::NumTraceActivityType ),
    m_fitActivity( false ),
    m_nuclide( nullptr ),
    m_activity( 0.0 ),
    m_relaxationDistance( 0.0f )
  {
  }
   
  
  void TraceSourceInfo::serialize( rapidxml::xml_node<char> *parent_node ) const
  {
    rapidxml::xml_document<char> *doc = parent_node->document();
    
    rapidxml::xml_node<char> * const base_node
                                      = doc->allocate_node( rapidxml::node_element, "TraceSource" );
    parent_node->append_node( base_node );
    
    const char *value = m_nuclide ? m_nuclide->symbol.c_str() : "None";
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, "Nuclide", value );
    base_node->append_node( node );
    
    const string actstr = PhysicalUnits::printToBestActivityUnits(m_activity,12);
    value = doc->allocate_string( actstr.c_str() );
    node = doc->allocate_node( rapidxml::node_element, "DisplayActivity", value );
    base_node->append_node( node );
    
    //// Put the total activity as an attribute, just to make a little more human readable/checkable
    //const string total_act = PhysicalUnits::printToBestActivityUnits( m_currentTotalActivity );
    //value = doc->allocate_string( total_act.c_str() );
    //rapidxml::xml_attribute<char> *attr = doc->allocate_attribute( "TotalActivity", value );
    //node->append_attribute( attr );
    
    value = GammaInteractionCalc::to_str( m_type );
    node = doc->allocate_node( rapidxml::node_element, "TraceActivityType", value );
    base_node->append_node( node );
      
    value = m_fitActivity ? "1" : "0";
    node = doc->allocate_node( rapidxml::node_element, "AllowFitting", value );
    base_node->append_node( node );
    
    if( m_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
    {
      const string relax_dist = PhysicalUnits::printToBestLengthUnits(m_relaxationDistance,7);
      value = doc->allocate_string( relax_dist.c_str() );
      node = doc->allocate_node( rapidxml::node_element, "RelaxationDistance", value );
      base_node->append_node( node );
    }
  }//void serialize( rapidxml::xml_node<char> *parent_node ) const
  
  
  void TraceSourceInfo::deSerialize( const rapidxml::xml_node<char> *trace_node )
  {
    using GammaInteractionCalc::TraceActivityType;
    
    if( !trace_node || !trace_node->value()
       || rapidxml::internal::compare(trace_node->value(), trace_node->value_size(),"TraceSource", 11, true) )
      throw runtime_error( "TraceSourceInfo::deSerialize: called with invalid node" );
      
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    
    const rapidxml::xml_node<char> *nuc_node = trace_node->first_node( "Nuclide" );
    const rapidxml::xml_node<char> *disp_act_node = trace_node->first_node( "DisplayActivity" );
    const rapidxml::xml_node<char> *act_type_node = trace_node->first_node( "TraceActivityType" );
    const rapidxml::xml_node<char> *allow_fit_node = trace_node->first_node( "AllowFitting" );
    const rapidxml::xml_node<char> *relax_node = trace_node->first_node( "RelaxationDistance" );
    
    if( !nuc_node || !disp_act_node || !act_type_node || !allow_fit_node )
      throw runtime_error( "TraceSrcDisplay::deSerialize: missing node" );

    const string act_type = SpecUtils::xml_value_str(act_type_node);
    
    m_type = TraceActivityType::NumTraceActivityType;
    for( TraceActivityType t = TraceActivityType(0);
        t != TraceActivityType::NumTraceActivityType;
        t = TraceActivityType(static_cast<int>(t) + 1) )
    {
      if( act_type == GammaInteractionCalc::to_str(t) )
      {
        m_type = t;
        break;
      }
    }//for( loop over TraceActivityTypes )
    
    if( m_type == TraceActivityType::NumTraceActivityType )
      throw runtime_error( "TraceSrcDisplay::deSerialize: Trace activity type in XML"
                          " ('" + act_type + "'), is invalid." );
    
    const string nuclidestr = SpecUtils::xml_value_str(nuc_node);
    m_nuclide = db->nuclide( nuclidestr );
    
    string acttxt = SpecUtils::xml_value_str(disp_act_node);
    
    // Check activity text is actually valid, and convert to users current prefered units
    try
    {
      m_activity = PhysicalUnits::stringToActivity(acttxt);
    }catch(...)
    {
      throw runtime_error( "TraceSrcDisplay::deSerialize: invalid activity: " + acttxt );
    }
    
    if( m_type == TraceActivityType::ExponentialDistribution )
    {
      if( !relax_node || !relax_node->value_size() )
        throw runtime_error( "TraceSrcDisplay::deSerialize: relaxation distance not specified." );
      
      try
      {
        const string relaxstr = SpecUtils::xml_value_str(relax_node);
        m_relaxationDistance = PhysicalUnits::stringToDistance( relaxstr );
      }catch( std::exception & )
      {
        throw runtime_error( "TraceSrcDisplay::deSerialize: relaxation distance invalid distance." );
      }
    }//if( m_type == TraceActivityType::ExponentialDistribution )
    
    const string do_fit_xml = SpecUtils::xml_value_str(allow_fit_node);
    m_fitActivity = ((do_fit_xml == "1") || SpecUtils::iequals_ascii(do_fit_xml,"true"));
  }//void TraceSourceInfo::deSerialize( const rapidxml::xml_node<char> *trace_node )
    
  
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  void TraceSourceInfo::equalEnough( const TraceSourceInfo &lhs, const TraceSourceInfo &rhs )
  {
    if( lhs.m_type != rhs.m_type )
      throw runtime_error( "TraceSourceInfo LHS Type != RHS Type" );
    
    if( lhs.m_fitActivity != rhs.m_fitActivity )
      throw runtime_error( "TraceSourceInfo LHS FitActivity != RHS FitActivity" );
    
    if( lhs.m_nuclide != rhs.m_nuclide )
      throw runtime_error( "TraceSourceInfo LHS Nuclide != RHS Nuclide" );
    
    const double act_diff = fabs( lhs.m_activity - rhs.m_activity );
    if( (act_diff > 1.0E-14) && (act_diff > 1.0E-8*max(fabs(lhs.m_activity), fabs(rhs.m_activity)) ) )
      throw runtime_error( "TraceSourceInfo LHS Activity != RHS Activity" );
    
    if( lhs.m_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
    {
      const float act_diff = fabs( lhs.m_relaxationDistance - rhs.m_relaxationDistance );
      if( (act_diff > 1.0E-7) && (act_diff > 1.0E-6*max(fabs(lhs.m_activity), fabs(rhs.m_activity)) ) )
        throw runtime_error( "TraceSourceInfo LHS RelaxationDistance != RHS RelaxationDistance" );
    }
  }//void TraceSourceInfo::equalEnough( const TraceSourceInfo &lhs, const TraceSourceInfo &rhs )
#endif //#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  

  
  ShieldingInfo::ShieldingInfo()
  : m_geometry( GammaInteractionCalc::GeometryType::NumGeometryType ),
    m_isGenericMaterial( false ),
    m_forFitting( false ),
    m_material( nullptr ),
    m_dimensions{ 0.0f, 0.0f, 0.0f },
    m_fitDimensions{ false, false, false },
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    m_truthDimensions{},
    m_truthDimensionsTolerances{},
#endif
    m_fitMassFrac( false ),
    m_nuclideFractions{},
    m_traceSources{}
  {
  }
    
  
  rapidxml::xml_node<char> *ShieldingInfo::serialize( rapidxml::xml_node<char> *parent_node ) const
  {
    rapidxml::xml_document<char> *doc = parent_node->document();
    
    const char *name, *value;
    rapidxml::xml_node<char> *base_node, *node;
    rapidxml::xml_attribute<char> *attr;
    
    name = "Shielding";
    base_node = doc->allocate_node( rapidxml::node_element, name );
    parent_node->append_node( base_node );
    
    //If you change the available options or formatting or whatever, increment the
    //  version field of the XML!
    string versionstr = std::to_string(ShieldingInfo::sm_xmlSerializationMajorVersion)
                        + "." + std::to_string(ShieldingInfo::sm_xmlSerializationMinorVersion);
    value = doc->allocate_string( versionstr.c_str() );
    attr = doc->allocate_attribute( "version", value );
    base_node->append_attribute( attr );
    
    name = "ForFitting";
    value = m_forFitting ? "1" : "0";
    node = doc->allocate_node( rapidxml::node_element, name, value );
    base_node->append_node( node );

    if( m_isGenericMaterial )
    {
      char buffer[32] = { '\0' };
      
      rapidxml::xml_node<> *generic_node;
      
      name = "Generic";
      generic_node = doc->allocate_node( rapidxml::node_element, name );
      base_node->append_node( generic_node );
      
      name = "ArealDensity";
      snprintf( buffer, sizeof(buffer), "%.9g", m_dimensions[1] );
      value = doc->allocate_string( buffer );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      generic_node->append_node( node );
      if( m_forFitting )
      {
        value = m_fitDimensions[1] ? "1" : "0";
        attr = doc->allocate_attribute( "Fit", value );
        node->append_attribute( attr );
      }//if( m_fitArealDensityCB )
      
      name = "AtomicNumber";
      snprintf( buffer, sizeof(buffer), "%.9g", m_dimensions[0] );
      value = doc->allocate_string( buffer );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      generic_node->append_node( node );
      if( m_forFitting )
      {
        value = m_fitDimensions[0] ? "1" : "0";
        attr = doc->allocate_attribute( "Fit", value );
        node->append_attribute( attr );
      }//if( m_fitAtomicNumberCB )
      
  #if( INCLUDE_ANALYSIS_TEST_SUITE )
      auto addTruth = [doc,generic_node]( const char *truthName, const boost::optional<double> &value ){
        if( value )
        {
          const string strval = std::to_string(*value);
          const char *value = doc->allocate_string( strval.c_str() );
          rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, truthName, value );
          generic_node->append_node( node );
        }
      };//addTruth(...)
      addTruth( "TruthAD", m_truthDimensions[1] );
      addTruth( "TruthADTolerance", m_truthDimensionsTolerances[1] );
      addTruth( "TruthAN", m_truthDimensions[0] );
      addTruth( "TruthANTolerance", m_truthDimensionsTolerances[0] );
  #endif
    }else
    {
      name = "Geometry";
      value = GammaInteractionCalc::to_str(m_geometry);
      node = doc->allocate_node( rapidxml::node_element, name, value );
      base_node->append_node( node );
      
      
      rapidxml::xml_node<> *material_node, *mass_frac_node, *iso_node;
      
      name = "Material";
      material_node = doc->allocate_node( rapidxml::node_element, name );
      base_node->append_node( material_node );
      
      name = "Name";
      string mat = m_material ? m_material->name : string("");
      value = doc->allocate_string( mat.c_str(), mat.size() + 1 );
      node = doc->allocate_node( rapidxml::node_element, name, value );
      material_node->append_node( node );
      
      //Lambda to add in dimension elements
      auto addDimensionNode = [this,doc,material_node]( const char *name, double val, bool fit )
                                                       -> rapidxml::xml_node<char> * {
        const string valstr = PhysicalUnits::printToBestLengthUnits( val, 9 );
                                                         
        const char *value = doc->allocate_string( valstr.c_str(), valstr.size()+1 );
        auto node = doc->allocate_node( rapidxml::node_element, name, value );
        material_node->append_node( node );
        if( m_forFitting )
        {
          value = fit ? "1" : "0";
          auto attr = doc->allocate_attribute( "Fit", value );
          node->append_attribute( attr );
        }//if( m_forFitting )
        
        return node;
      };//addDimensionNode lambda
      
      
      // For backward compatibility with XML serialization version 0.0, we will add in a <Thickness>
      //  element, corresponding to dimension along detector axis
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          node = addDimensionNode( "Thickness", m_dimensions[0], m_fitDimensions[0] );
          node->append_attribute( doc->allocate_attribute( "Remark", "SphericalThickness" ) );
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
          node = addDimensionNode( "Thickness", m_dimensions[1], m_fitDimensions[1] );
          node->append_attribute( doc->allocate_attribute( "Remark", "CylinderLengthThickness" ) );
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          node = addDimensionNode( "Thickness", m_dimensions[0], m_fitDimensions[0] );
          node->append_attribute( doc->allocate_attribute( "Remark", "CylinderRadiusThickness" ) );
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          node = addDimensionNode( "Thickness", m_dimensions[2], m_fitDimensions[2] );
          node->append_attribute( doc->allocate_attribute( "Remark", "RectangularDepthThickness" ) );
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      // We could write dimensions from all the WLineEdits into the XML, which would kinda keep state
      //  across changing geometries, but maybe for the moment we'll just write the relevant
      //  dimensions.  TODO: think about this a little more
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          addDimensionNode( "SphericalThickness", m_dimensions[0], m_fitDimensions[0] );
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          addDimensionNode( "CylinderRadiusThickness", m_dimensions[0], m_fitDimensions[0] );
          addDimensionNode( "CylinderLengthThickness", m_dimensions[1], m_fitDimensions[1] );
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          addDimensionNode( "RectangularWidthThickness",  m_dimensions[0], m_fitDimensions[0] );
          addDimensionNode( "RectangularHeightThickness", m_dimensions[1], m_fitDimensions[1] );
          addDimensionNode( "RectangularDepthThickness",  m_dimensions[2], m_fitDimensions[2] );
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      
      if( m_forFitting && m_fitMassFrac && !m_nuclideFractions.empty() )
      {
        name = "FitMassFraction";
        value = m_fitMassFrac ? "1" : "0";
        mass_frac_node = doc->allocate_node( rapidxml::node_element, name, value );
        material_node->append_node( mass_frac_node );
      }//if( m_forFitting )
      
      for( const auto &nuc_to_frac : m_nuclideFractions )
      {
        const SandiaDecay::Nuclide * const nuc = nuc_to_frac.first;
        const float frac = nuc_to_frac.second;
        
        if( nuc )
        {
          iso_node = doc->allocate_node( rapidxml::node_element, "Nuclide" );
          material_node->append_node( iso_node );
          
          value = doc->allocate_string( nuc->symbol.c_str() );
          node = doc->allocate_node( rapidxml::node_element, "Name", value );
          iso_node->append_node( node );
          
          const string fracstr = PhysicalUnits::printCompact(frac, 12);
          value = doc->allocate_string( fracstr.c_str(), fracstr.size() + 1 );
          node = doc->allocate_node( rapidxml::node_element, "MassFrac", value );
          iso_node->append_node( node );
        }//if( nuc )
      }//for( const auto &nuc_to_frac : nuc_to_fractions )
      
      for( const TraceSourceInfo &trace : m_traceSources )
      {
        trace.serialize( material_node );
      }//for( const TraceSourceInfo &trace : m_traceSources )
      
      
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
      auto addTruth = [doc,material_node]( const char *truthName, const boost::optional<double> &value ){
        if( value )
        {
          const string strval = PhysicalUnits::printToBestLengthUnits(*value,6);
          const char *value = doc->allocate_string( strval.c_str() );
          rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, truthName, value );
          material_node->append_node( node );
        }
      };//addTruth(...)
      
      addTruth( "TruthThickness", m_truthDimensions[0] );
      addTruth( "TruthThicknessTolerance", m_truthDimensionsTolerances[0] );
      addTruth( "TruthThicknessD2", m_truthDimensions[1] );
      addTruth( "TruthThicknessD2Tolerance", m_truthDimensionsTolerances[1] );
      addTruth( "TruthThicknessD3", m_truthDimensions[2] );
      addTruth( "TruthThicknessD3Tolerance", m_truthDimensionsTolerances[2] );
  #endif
    }//if( m_isGenericMaterial ) / else
    
    return base_node;
  }//void serialize( rapidxml::xml_node<char> *parent_node ) const
  
  
  void ShieldingInfo::deSerialize( const rapidxml::xml_node<char> *shield_node, MaterialDB *materialDb )
  {
    using GammaInteractionCalc::GeometryType;
    
    rapidxml::xml_attribute<char> *attr;
    rapidxml::xml_node<char> *node, *geom_node, *generic_node, *material_node;
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
   
    if( shield_node->name() != string("Shielding") )
      throw runtime_error( "ShieldingSelects XML node should be 'Shielding'" );
    
    attr = shield_node->first_attribute( "version", 7 );
    int version;
    if( !attr || !attr->value() || !(stringstream(attr->value())>>version) )
      throw runtime_error( "ShieldingSelects should be versioned" );
    
    if( version != sm_xmlSerializationMajorVersion )
      throw runtime_error( "Invalid XML version for ShieldingSelect" );
    
    // Note that version is either "0" (implied minor version 0), or "0.1" or something; we could read
    //  in the minor version to sm_xmlSerializationMinorVersion and use it, but theres really no need.
    
    // Determine if
    node = shield_node->first_node( "ForFitting", 10 );
    if( !node || !node->value() || !(stringstream(node->value()) >> m_forFitting) )
      throw runtime_error( "Missing/invalid for fitting node" );
    
    generic_node = shield_node->first_node( "Generic", 7 );
    material_node = shield_node->first_node( "Material", 8 );
    
    if( (!generic_node) == (!material_node) )
      throw runtime_error( "Shielding element must have exactly one Generic or Material node." );
      
    if( generic_node )
    {
      m_isGenericMaterial = true;
      m_geometry = GeometryType::NumGeometryType;
      m_fitMassFrac = false;
      m_nuclideFractions.clear();
      m_traceSources.clear();
      
      const rapidxml::xml_node<char> *ad_node, *an_node;
      ad_node = generic_node->first_node( "ArealDensity", 12 );
      an_node = generic_node->first_node( "AtomicNumber", 12 );
      
      if( !ad_node || !ad_node->value() || !an_node || !an_node->value() )
        throw runtime_error( "Generic material must have ArealDensity and"
                             " AtomicNumber nodes" );
      
      double ad, an;
      if( !(std::stringstream(an_node->value()) >> an) )
        throw runtime_error( "Invalid atomic number specified - not a floating point" );
      
      if( !(std::stringstream(ad_node->value()) >> ad) )
        throw runtime_error( "Invalid areal density - not a floating point" );
      
      m_dimensions[0] = an;
      m_dimensions[1] = ad;
      m_dimensions[2] = 0.0;
      
      m_fitDimensions[0] = m_fitDimensions[1] = m_fitDimensions[2] = false;
      
      if( m_forFitting  )
      {
        bool fit;
        attr = an_node->first_attribute( "Fit", 3 );
        if( attr && attr->value() && (stringstream(attr->value())>>fit) )
          m_fitDimensions[0] = fit;
        
        attr = ad_node->first_attribute( "Fit", 3 );
        if( attr && attr->value() && (stringstream(attr->value())>>fit) )
          m_fitDimensions[1] = fit;
      }//if( m_fitArealDensityCB )
      
  #if( INCLUDE_ANALYSIS_TEST_SUITE )
      auto getTruth = [generic_node]( const char *truthName, boost::optional<double> &value ){
        value.reset();
        auto node = generic_node->first_node( truthName );
        if( !node || !node->value() )
          return;
        
        double dblvalue;
        if( (stringstream(node->value()) >> dblvalue) )
          value = dblvalue;
        else
          cerr << "\n\nFailed to deserialize shielding " << truthName << " from " << node->value()
          << "\n\n" << endl;
      };//getTruth(...)
      
      getTruth( "TruthAN", m_truthDimensions[0] );
      getTruth( "TruthANTolerance", m_truthDimensionsTolerances[0] );
      getTruth( "TruthAD", m_truthDimensions[1] );
      getTruth( "TruthADTolerance", m_truthDimensionsTolerances[1] );
  #endif
    }//if( generic_node )
    
    if( material_node )
    {
      if( !materialDb )
        throw runtime_error( "De-serializing shielding requires valid Material DB" );
        
      m_isGenericMaterial = false;
      m_fitMassFrac = false;
      m_nuclideFractions.clear();
      m_traceSources.clear();
      
      //To be very compatible (pre non-spherical days), if no <Geometry> node, then we
      //  will assume spherical geometry.
      m_geometry = GeometryType::Spherical;
      geom_node = shield_node->first_node( "Geometry", 8 );
      if( geom_node && geom_node->value() )
      {
        bool matched = false;
        const string geomstr( geom_node->value(), geom_node->value() + geom_node->value_size() );
        for( GeometryType type = GeometryType(0);
            !matched && (type != GeometryType::NumGeometryType);
            type = GeometryType( static_cast<int>(type) + 1 ) )
        {
          matched = (geomstr == GammaInteractionCalc::to_str(type));
          if( matched )
            m_geometry = type;
        }//for( loop over GeometryTypes )
        
        if( !matched )
          throw runtime_error( "Invalid geometry type: '" + geomstr + "'" );
      }//if( geom_node )
      
      assert( m_geometry != GeometryType::NumGeometryType );
      
      const rapidxml::xml_node<> *dim_nodes[3] = { nullptr, nullptr, nullptr };

      const rapidxml::xml_node<> *name_node = material_node->first_node( "Name", 4 );
      if( !name_node || !name_node->value() )
        throw runtime_error( "Material node didnt have name node" );
      
      const string material_name( name_node->value(), name_node->value() + name_node->value_size() );
      m_material.reset();
      if( !material_name.empty() )
      {
        const Material *mat = materialDb->material( material_name );
        if( !mat )
          throw runtime_error( "Invalid shielding material: '" + material_name + "'" );
        m_material = make_shared<Material>( *mat );
      }//if( !material_name.empty() )
      
      
      int required_dim = 0;
      switch( m_geometry )
      {
        case GeometryType::Spherical:
          required_dim = 1;
          dim_nodes[0] = XML_FIRST_NODE(material_node, "SphericalThickness");
          if( !dim_nodes[0] )
            dim_nodes[0] = XML_FIRST_NODE(material_node, "Thickness");
          break;
          
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn:
          required_dim = 2;
          dim_nodes[0] = XML_FIRST_NODE(material_node, "CylinderRadiusThickness");
          dim_nodes[1] = XML_FIRST_NODE(material_node, "CylinderLengthThickness");
          break;
          
        case GeometryType::Rectangular:
          required_dim = 3;
          dim_nodes[0] = XML_FIRST_NODE(material_node, "RectangularWidthThickness");
          dim_nodes[1] = XML_FIRST_NODE(material_node, "RectangularHeightThickness");
          dim_nodes[2] = XML_FIRST_NODE(material_node, "RectangularDepthThickness");
          break;
          
        case GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      
      for( size_t i = 0; i < 3; ++i )
      {
        m_dimensions[i] = 0.0;
        m_fitDimensions[i] = false;
      }//for( size_t i = 0; i < 3; ++i )
      
      
      for( int i = 0; i < required_dim; ++i )
      {
        // Check and make sure all the dimension node values are valid
        if( !dim_nodes[i] || !dim_nodes[i]->value() )
          throw runtime_error( "Missing required dimension node" );
        
        const string val( dim_nodes[i]->value(), dim_nodes[i]->value() + dim_nodes[i]->value_size() );
        try
        {
          m_dimensions[i] = PhysicalUnits::stringToDistance(val);
          if( m_dimensions[i] < 0.0 )
            throw runtime_error( "" );
        }catch( std::exception & )
        {
          throw runtime_error( "Invalid dimension given '" + val + "'" );
        }
        
        if( m_forFitting )
        {
          const auto attr = dim_nodes[i]->first_attribute( "Fit", 3 );
          if( !attr || !(stringstream(attr->value()) >> m_fitDimensions[i]) )
            throw runtime_error( "Material node expected thickness Fit attribute" );
        }
      }//for( int i = 0; i < required_dim; ++i )
      
#if( INCLUDE_ANALYSIS_TEST_SUITE )
      auto getTruth = [material_node]( const char *truthName, boost::optional<double> &value ){
        value.reset();
        auto node = material_node->first_node( truthName );
        if( !node || !node->value() )
          return;
      
        const string valstr = SpecUtils::xml_value_str( node );
        try
        {
          value = PhysicalUnits::stringToDistance( valstr );
        }catch( std::exception &e )
        {
          throw runtime_error( "Shielding truth-thickness is invalid distance: '" + valstr + "'" );
        }
      };//getTruth(...)
    
      getTruth( "TruthThickness", m_truthDimensions[0] );
      getTruth( "TruthThicknessTolerance", m_truthDimensionsTolerances[0] );
      getTruth( "TruthThicknessD2", m_truthDimensions[1] );
      getTruth( "TruthThicknessD2Tolerance", m_truthDimensionsTolerances[1] );
      getTruth( "TruthThicknessD3", m_truthDimensions[2] );
      getTruth( "TruthThicknessD3Tolerance", m_truthDimensionsTolerances[2] );
#endif //INCLUDE_ANALYSIS_TEST_SUITE
      
      if( m_forFitting )
      {
        const rapidxml::xml_node<> *fitmassfrac_node = XML_FIRST_NODE(material_node, "FitMassFraction");
        if( fitmassfrac_node && fitmassfrac_node->value() )
        {
          const string valstr = SpecUtils::xml_value_str( fitmassfrac_node );
          if( !(stringstream(valstr) >> m_fitMassFrac) )
            throw runtime_error( "Invalid FitMassFraction value '" + valstr + "'" );
        }//if( m_forFitting )
      }//if( m_forFitting )

      double last_frac = 0.0;
      const SandiaDecay::Nuclide *last_nuc = NULL;
      
      for( const rapidxml::xml_node<> *iso_node = XML_FIRST_NODE(material_node, "Nuclide");
          iso_node; iso_node = XML_NEXT_TWIN(iso_node) )
      {
        name_node = XML_FIRST_NODE( iso_node, "Name");
        const rapidxml::xml_node<> * const frac_node = XML_FIRST_NODE( iso_node, "MassFrac" );
        
        const string namestr = SpecUtils::xml_value_str( name_node );
        const string fractionstr = SpecUtils::xml_value_str( frac_node );
        if( namestr.empty() || fractionstr.empty() )
          throw runtime_error( "Missing invalid name/mass frac node form iso" );
        
        const SandiaDecay::Nuclide *nuc = db->nuclide( namestr );
        const SandiaDecay::Element *el = nuc ? db->element( nuc->atomicNumber ) : nullptr;
        if( !nuc || !el )
          throw runtime_error( namestr + " is not a valid isotope" );
        
        double fraction;
        if( !(stringstream(fractionstr) >> fraction) )
          throw runtime_error( "Invalid mass fraction: '" + fractionstr + "'" );
        
        if( m_nuclideFractions.count(nuc) )
          throw runtime_error( "Nuclide '" + namestr + "' specified multiple times." );
        
        if( !m_material )
          throw runtime_error( "Mass-fraction specified, without a valid shielding material" );
        
        bool valid_element = false;
        for( size_t i = 0; !valid_element && (i < m_material->elements.size()); ++i )
          valid_element = (m_material->elements[i].first == el);
        
        // I dont recal if m_material->nuclides are a subset of elements that are in
        //  m_material->elements so we'll check, just to be sure, for the moment
        for( size_t i = 0; !valid_element && (i < m_material->nuclides.size()); ++i )
          valid_element = (m_material->nuclides[i].first == nuc);
        
        if( !valid_element )
          throw runtime_error( "A mass-fraction for '" + nuc->symbol
                              + "' was specified, but material '" + m_material->name
                              + "' doesnt have that element." );
        
        m_nuclideFractions[nuc] = fraction;
      }//for( loop over isotope nodes )
     
        
      for( auto trace_node = XML_FIRST_NODE(material_node, "TraceSource");
            trace_node; trace_node = XML_NEXT_TWIN(trace_node) )
      {
        TraceSourceInfo trace;
        trace.deSerialize( trace_node );
        m_traceSources.push_back( std::move(trace) );
      }//for( loop over trace nodes )
    }//if( material_node )
  }//void deSerialize( const rapidxml::xml_node<char> *shield_node, MaterialDB *materialDb )
    

  std::string ShieldingInfo::encodeStateToUrl() const
  {
    // "V=1&G=S&D1=1.2cm&N=Fe&NTRACE=3&TRACEN=1&V=1&N=U238&A=1.2uCi&T=total&F=1"
    // "V": version
    // "G": geometry
    // "F": for fitting; if not specified than false
    // "D1": "D2": Thickness, depth, etc
    // "FD1": "FD2": fit the cooresponding dimensions
    // "N": material name
    // "AN": atomic number
    // "FAN": fit atomic number - if not specified than false
    // "AD": areal density
    // "FAD": fit areal density - if not specified than false
    // ...
    
    string answer = "V=1";
    
    if( m_forFitting )
      answer += "&F=1";
    
    if( m_isGenericMaterial )
    {
      answer += "&AD=" + PhysicalUnits::printCompact(m_dimensions[1], 7);
      answer += "&AN=" + PhysicalUnits::printCompact(m_dimensions[0], 7);
      
      if( m_forFitting && m_fitDimensions[0] )
        answer += "FAN=1";
      if( m_forFitting && m_fitDimensions[1] )
        answer += "FAD=1";
    }else
    {
      std::string material = m_material ? m_material->name : string();
      SpecUtils::ireplace_all(material, "#", "%23" );
      SpecUtils::ireplace_all(material, "&", "%26" );
      
      answer += "&N=" + material;
      
      using PhysicalUnits::cm;
        
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          answer += "&G=S&D1=" + PhysicalUnits::printCompact(m_dimensions[0]/cm, 7) + "cm";
          if( m_forFitting && m_fitDimensions[0] )
            answer += "&FD1=1";
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          if( m_geometry == GammaInteractionCalc::GeometryType::CylinderEndOn )
            answer += "&G=CE";
          else
            answer += "&G=CS";
          answer += "&D1=" + PhysicalUnits::printCompact(m_dimensions[0]/cm, 7) + "cm";
          answer += "&D2=" + PhysicalUnits::printCompact(m_dimensions[1]/cm, 7) + "cm";
          if( m_forFitting && m_fitDimensions[0] )
            answer += "&FD1=1";
          if( m_forFitting && m_fitDimensions[1] )
            answer += "&FD2=1";
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          answer += "&G=R";
          answer += "&D1=" + PhysicalUnits::printCompact(m_dimensions[0]/cm, 7) + "cm";
          answer += "&D2=" + PhysicalUnits::printCompact(m_dimensions[1]/cm, 7) + "cm";
          answer += "&D3=" + PhysicalUnits::printCompact(m_dimensions[2]/cm, 7) + "cm";
          if( m_forFitting && m_fitDimensions[0] )
            answer += "&FD1=1";
          if( m_forFitting && m_fitDimensions[1] )
            answer += "&FD2=1";
          if( m_forFitting && m_fitDimensions[2] )
            answer += "&FD3=1";
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      // TODO: encode self-attenuating and trace sources, and maybe "truth" value; with a string like "NTRACE=3&TRACEN=1&V=1&N=U238&A=1.2uCi&T=total&F=1"
    }//if( m_isGenericMaterial ) / else
    
    return answer;
  }//std::string encodeStateToUrl() const
    
  
  void ShieldingInfo::handleAppUrl( std::string query_str, MaterialDB *materialDb )
  {
    // "V=1&G=S&D1=1.2cm&N=Fe"
    // "V": version
    // "G": geometry
    // "F": for fitting; if not specified than false
    // "D1": "D2": Thickness, depth, etc
    // "FD1": "FD2": fit the cooresponding dimensions
    // "N": material name
    // "AN": atomic number
    // "FAN": fit atomic number - if not specified than false
    // "AD": areal density
    // "FAD": fit areal density - if not specified than false
    // ...
    
    const map<string,string> values = AppUtils::query_str_key_values( query_str );
    
    auto iter = values.find( "V" );
    if( (iter == end(values)) || ((iter->second != "1") && !SpecUtils::istarts_with(iter->second, "1.")) )
      throw runtime_error( "ShieldingInfo URI not compatible version" );
    
    iter = values.find( "F" );
    m_forFitting = (iter == end(values)) ? false : (iter->second == "1");
    
    const bool generic = values.count( "AN" );
    m_isGenericMaterial = generic;
    
    if( generic )
    {
      m_material = nullptr;
      m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
      
      // Get all the values, and validate, before setting valus
      iter = values.find( "AN" );
      assert( iter != end(values) ); //If we are here, we know `values` cintains 'AN'
      const string an_str = iter->second;
      
      if( (iter = values.find( "AD" )) == end(values) )
        throw runtime_error( "ShieldingSelect missing AD in URI." );
      
      const string ad_str = iter->second;
      
      // Validate AN and AD are numbers, within acceptable range
      try
      {
        float an, ad;
        if( !SpecUtils::parse_float( an_str.c_str(), an_str.size(), an) )
          throw runtime_error( "AN '" + an_str + "' not float" );
        if( !SpecUtils::parse_float( ad_str.c_str(), ad_str.size(), ad) )
          throw runtime_error( "AD '" + ad_str + "' not float" );
        
        if( (an < MassAttenuation::sm_min_xs_atomic_number)
           || (an > MassAttenuation::sm_max_xs_atomic_number) )
          throw runtime_error( "AN '" + an_str + "' is out of range" );
        
        if( (ad < 0.0f)
           || (ad > static_cast<float>(GammaInteractionCalc::sm_max_areal_density_g_cm2)) )
          throw runtime_error( "AD '" + ad_str + "' is out of range" );
        
        m_dimensions[0] = an;
        m_dimensions[1] = ad;
        m_dimensions[2] = 0.0f;
      }catch( std::exception &e )
      {
        throw runtime_error( "ShieldingSelect invalid AN or AD in URI: " + string(e.what()) );
      }
      
      iter = values.find( "FAN" );
      const bool fitAN = (!m_forFitting || (iter == end(values))) ? false : (iter->second == "1");
      
      iter = values.find( "FAD" );
      const bool fitAD = (!m_forFitting || (iter == end(values))) ? false : (iter->second == "1");
      
      m_fitDimensions[0] = fitAN;
      m_fitDimensions[1] = fitAD;
      m_fitDimensions[2] = false;
    }else
    {
      if( (iter = values.find( "G" )) == end(values) )
        throw runtime_error( "ShieldingSelect missing geometry in URI." );
      
      const string &geom_str = iter->second;
      m_geometry = GammaInteractionCalc::GeometryType::Spherical;
      
      if( SpecUtils::iequals_ascii(geom_str, "S") )
        m_geometry = GammaInteractionCalc::GeometryType::Spherical;
      else if( SpecUtils::iequals_ascii(geom_str, "CE") )
        m_geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
      else if( SpecUtils::iequals_ascii(geom_str, "CS") )
        m_geometry = GammaInteractionCalc::GeometryType::CylinderSideOn;
      else if( SpecUtils::iequals_ascii(geom_str, "R") )
        m_geometry = GammaInteractionCalc::GeometryType::Rectangular;
      else
        throw runtime_error( "ShieldingSelect invalid geometry '" + geom_str + "' in URI" );
      
      string d1str, d2str, d3str;
      bool fitD1 = false, fitD2 = false, fitD3 = false;
      
      auto getDim = [&]( const string &key, double &val, bool &fit ) {
        if( (iter = values.find(key)) == end(values) )
          throw runtime_error( "ShieldingSelect missing dimension '" + key + "' in URI." );

        const string &valstr = iter->second;
        if( !valstr.empty() )
        {
          try
          {
            val = static_cast<float>( PhysicalUnits::stringToDistance(valstr) );
            if( val < 0.0 )
              throw runtime_error( "Distance must be positive." );
          }catch( std::exception & )
          {
            throw runtime_error( "ShieldingSelect invalid dimension '" + valstr + "' for '" + key + "'" );
          }
        }//if( !valstr.empty() )
        
        iter = values.find( "F" + key );
        fit = (!m_forFitting || (iter == end(values))) ? false : (iter->second == "1");
      };//getDim
      
      
      for( size_t i = 0; i < 3; ++i )
      {
        m_dimensions[i] = 0.0f;
        m_fitDimensions[i] = false;
      }
      
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Rectangular:
           getDim("D3", m_dimensions[2], m_fitDimensions[2]);
          // drop-through intentional
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
           getDim("D2", m_dimensions[1], m_fitDimensions[1]);
          // drop-through intentional
        case GammaInteractionCalc::GeometryType::Spherical:
          getDim("D1", m_dimensions[0], m_fitDimensions[0]);
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( geometry )
      
      if( (iter = values.find( "N" )) == end(values) )
        throw runtime_error( "ShieldingSelect missing material name in URI." );
      
      string materialstr = iter->second;
      SpecUtils::ireplace_all(materialstr, "%23", "#" );
      SpecUtils::ireplace_all(materialstr, "%26", "&" );
      
      if( materialstr.empty() )
      {
        m_material = nullptr;
      }else
      {
        if( !materialDb )
          throw runtime_error( "ShieldingInfo::handleAppUrl: invalid MaterialDB passed in" );
        
        const Material * mat = materialDb->material(materialstr);
        
        if( !mat )
        {
          try
          {
            const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
            mat = materialDb->parseChemicalFormula( materialstr, db );
          }catch( std::exception & )
          {
          }
        }//if( !mat )
        
        if( !mat )
          throw runtime_error( "Invalid material name '" + materialstr + "'" );
        m_material = make_shared<Material>( *mat );
      }//if( !materialstr.empty() ) / else
      
      
      // Self-atten source stuff
      //bool m_fitMassFrac;
      //std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractions;
      
      // Trace-source stuff
      //std::vector<TraceSourceInfo> m_traceSources;
    }//if( generic ) / else
  }//void ShieldingInfo::handleAppUrl( std::string query_str )
    
  
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  void ShieldingInfo::equalEnough( const ShieldingInfo &lhs, const ShieldingInfo &rhs )
  {
    if( lhs.m_geometry != rhs.m_geometry )
      throw runtime_error( "ShieldingInfo LHS Geometry != RHS Geometry" );
    
    if( lhs.m_isGenericMaterial != rhs.m_isGenericMaterial )
      throw runtime_error( "ShieldingInfo LHS IsGenericMaterial != RHS IsGenericMaterial" );
    
    if( lhs.m_forFitting != rhs.m_forFitting )
      throw runtime_error( "ShieldingInfo LHS ForFitting != RHS ForFitting" );

    if( (!lhs.m_material) != (!rhs.m_material) )
      throw runtime_error( "ShieldingInfo LHS material validity != RHS material validity" );
    
    if( lhs.m_material && (lhs.m_material->name != rhs.m_material->name) )
      throw runtime_error( "ShieldingInfo LHS material name != RHS material name" );
    
    
    
    int ndim = 0;
    switch( lhs.m_geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        ndim = 1;
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        ndim = 2;
        break;
      
      case GammaInteractionCalc::GeometryType::Rectangular:
        ndim = 3;
        break;
        
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        if( lhs.m_isGenericMaterial )
          ndim = 2;
        break;
    }//switch( lhs.m_geometry )
    
    for( int i = 0; i < ndim; ++i )
    {
      {
        const float m = std::max( fabs(lhs.m_dimensions[i]), fabs(rhs.m_dimensions[i]) );
        const float d = fabs(lhs.m_dimensions[i] - rhs.m_dimensions[i]);
        if( (m > 1.0E-9) && (fabs(d) > 1.0E-6*m) )
          throw runtime_error( "ShieldingInfo LHS dimension " + std::to_string(i) + " differs from RHS." );
      }
      
      if( lhs.m_fitDimensions[i] != rhs.m_fitDimensions[i] )
        throw runtime_error( "ShieldingInfo LHS FitDimension " + std::to_string(i) + " differs from RHS." );
      
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
      if( lhs.m_truthDimensions[i].has_value() != rhs.m_truthDimensions[i].has_value() )
        throw runtime_error( "ShieldingInfo LHS TruthDimension validity != RHS TruthDimension validity" );
      
      if( lhs.m_truthDimensions[i].has_value() )
      {
        const float lval = *lhs.m_truthDimensions[i];
        const float rval = *rhs.m_truthDimensions[i];
        const float mval = std::max( fabs(lval), fabs(rval) );
        const float dval = fabs(lval - rval);
        
        if( (mval > 1.0E-9) && (fabs(dval) > 1.0E-6*mval) )
          throw runtime_error( "ShieldingInfo LHS TruthValue != RHS TruthValue" );
      }//if( lhs.m_truthDimensions[i].has_value() )
      
      
      if( lhs.m_truthDimensionsTolerances[i].has_value() != rhs.m_truthDimensionsTolerances[i].has_value() )
        throw runtime_error( "ShieldingInfo LHS TruthDimensionTolerance validity != RHS TruthDimensionTolerance validity" );
      
      if( lhs.m_truthDimensionsTolerances[i].has_value() )
      {
        const float lval = *lhs.m_truthDimensionsTolerances[i];
        const float rval = *rhs.m_truthDimensionsTolerances[i];
        const float mval = std::max( fabs(lval), fabs(rval) );
        const float dval = fabs(lval - rval);
        
        if( (mval > 1.0E-9) && (fabs(dval) > 1.0E-6*mval) )
          throw runtime_error( "ShieldingInfo LHS TruthTolerance != RHS TruthTolerance" );
      }//if( lhs.m_truthDimensions[i].has_value() )
#endif
    }//for( int i = 0; i < ndim; ++i )
    
    
    // Self-atten source stuff
    if( lhs.m_fitMassFrac != rhs.m_fitMassFrac )
      throw runtime_error( "ShieldingInfo LHS FitMassFraction != RHS FitMassFraction" );
    
    if( lhs.m_nuclideFractions.size() != rhs.m_nuclideFractions.size() )
      throw runtime_error( "ShieldingInfo LHS NumSourceElements != RHS NumSourceElements" );
    
    for( const auto &nucFrac : lhs.m_nuclideFractions )
    {
      const auto rfracpos = rhs.m_nuclideFractions.find(nucFrac.first);
      if( rfracpos == end(rhs.m_nuclideFractions) )
        throw runtime_error( "ShieldingInfo LHS has nuclide "
                            + (nucFrac.first ? nucFrac.first->symbol : string(""))
                            + " as source, but RHS doesnt" );
      
      const float lval = nucFrac.second;
      const float rval = rfracpos->second;
      const float mval = std::max( fabs(lval), fabs(rval) );
      const float dval = fabs(lval - rval);
      
      if( (mval > 1.0E-9) && (fabs(dval) > 1.0E-6*mval) )
        throw runtime_error( "ShieldingInfo LHS MassFraction for "
                            + (nucFrac.first ? nucFrac.first->symbol : string(""))
                            + " != RHS MassFraction ("
                            + std::to_string(lval) + " vs " + std::to_string(rval) + ")" );
    }//for( const auto &nucFrac : lmfrac )
    
    // Trace-source stuff
    if( lhs.m_traceSources.size() != rhs.m_traceSources.size() )
      throw runtime_error( "ShieldingInfo LHS NumTraceSources != RHS NumTraceSources" );
    
    for( size_t i = 0; i < lhs.m_traceSources.size(); ++i )
      TraceSourceInfo::equalEnough( lhs.m_traceSources[i], rhs.m_traceSources[i] );
  }//void equalEnough( const ShieldingInfo &lhs, const ShieldingInfo &rhs )
#endif //PERFORM_DEVELOPER_CHECKS
  
  
  
}//namespace ShieldingSourceFitCalc
