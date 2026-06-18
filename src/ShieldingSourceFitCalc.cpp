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
#include <optional>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include <Wt/WServer.h>
#include <Wt/WApplication.h>

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnUserParameters.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameterState.h"

#include "ceres/ceres.h"

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
#include "InterSpec/GammaInteractionCalc_imp.hpp"


using namespace std;

namespace ShieldingSourceFitCalc
{
  const int ShieldingInfo::sm_xmlSerializationMajorVersion = 0;

  /** Change log:
   - 20230705: no version changed, but converted from being member of ShieldingSelect to ShieldingInfo
   - 20241209: Depreciated <FitMassFraction> element, and made this a per-nuclide quantity by adding a <Fit> element under
     the <Nuclide> node.  Also added a "FitOtherFraction" attribute
   */
  const int ShieldingInfo::sm_xmlSerializationMinorVersion = 2;

  
  /** Change log:
   - 20230712, no version change: split this out from ShieldingSourceDisplay::sm_xmlSerializationMajorVersion
   */
  const int SourceFitDef::sm_xmlSerializationMajorVersion = 0;
  
  /** Change log:
   - 20211031, Version 1: added "SourceType" node under the <Nuclide> to say whether its a point, intrinsic, or trace source.
                    added "Geometry" node under the "ShieldingSourceFit" node
   - 20230712, no version change: split this out from ShieldingSourceDisplay::sm_xmlSerializationMinorVersion
   */
  const int SourceFitDef::sm_xmlSerializationMinorVersion = 1;

  
  
SourceFitDef::SourceFitDef()
    : nuclide( nullptr ), activity(0.0), fitActivity(false),
      age(0.0), fitAge(false), ageDefiningNuc( nullptr ),
      sourceType( ShieldingSourceFitCalc::ModelSourceType::Point )
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
      , truthActivity(), truthActivityTolerance(),
      truthAge(), truthAgeTolerance()
  #endif
{
}

#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
bool close_enough( const double l, const double r )
{
  const double diff = fabs(l - r);
  return (diff < 1.0E-12) || (diff < 1.0E-6*std::max(fabs(l), fabs(r)));
};
  
bool close_enough( const boost::optional<double> &l, const boost::optional<double> &r )
{
  if( l.has_value() != r.has_value() )
    return false;
  
  if( !l.has_value() )
    return true;
  
  const double diff = fabs((*l) - (*r));
  return (diff < 1.0E-12) || (diff < 1.0E-6*std::max(fabs(*l), fabs(*r)));
};

bool close_enough( const std::optional<double> &l, const std::optional<double> &r )
{
  if( l.has_value() != r.has_value() )
    return false;
  
  if( !l.has_value() )
    return true;
  
  const double diff = fabs((*l) - (*r));
  return (diff < 1.0E-12) || (diff < 1.0E-6*std::max(fabs(*l), fabs(*r)));
};
  
  
void SourceFitDef::equalEnough( const SourceFitDef &lhs, const SourceFitDef &rhs )
{
  if( lhs.nuclide != rhs.nuclide )
    throw runtime_error( "SourceFitDef LHS nuclide != RHS nuclide" );
  
  if( !close_enough(lhs.activity, rhs.activity) )
  {
    const double diff = fabs(lhs.activity - rhs.activity);
    const double max_val = std::max(fabs(lhs.activity), fabs(rhs.activity));
    const double abs_tol = 1.0E-12;
    const double rel_tol = 1.0E-6 * max_val;
    char buffer[512];
    snprintf( buffer, sizeof(buffer),
              "SourceFitDef LHS activity (%.17g) != RHS activity (%.17g): "
              "diff=%.17g, abs_tol=%.17g, rel_tol=%.17g, max_val=%.17g, rel_diff=%.6f%%",
              lhs.activity, rhs.activity, diff, abs_tol, rel_tol, max_val,
              100.0 * diff / max_val );
    throw runtime_error( buffer );
  }
  
  if( lhs.activityUncertainty )
    assert( lhs.activityUncertainty.value() > 0.0 );
  if( rhs.activityUncertainty )
    assert( rhs.activityUncertainty.value() > 0.0 );
  
  if( lhs.fitActivity != rhs.fitActivity )
    throw runtime_error( "SourceFitDef LHS fitActivity != RHS fitActivity" );
  
  if( !close_enough(lhs.age, rhs.age) )
    throw runtime_error( "SourceFitDef LHS age != RHS age" );
  
  if( lhs.fitAge != rhs.fitAge )
    throw runtime_error( "SourceFitDef LHS fitAge != RHS fitAge" );
  
  if( lhs.ageDefiningNuc != rhs.ageDefiningNuc )
    throw runtime_error( "SourceFitDef LHS ageDefiningNuc != RHS ageDefiningNuc" );
  
  if( lhs.sourceType != rhs.sourceType )
    throw runtime_error( "SourceFitDef LHS sourceType != RHS sourceType" );
  
  if( !close_enough( lhs.activityUncertainty, rhs.activityUncertainty ) )
  {
    const bool lhs_has = lhs.activityUncertainty.has_value();
    const bool rhs_has = rhs.activityUncertainty.has_value();
    if( lhs_has != rhs_has )
    {
      throw runtime_error( "SourceFitDef LHS activityUncertainty has_value (" + std::to_string(lhs_has) + ") "
                          "!= RHS activityUncertainty has_value (" + std::to_string(rhs_has) + ")" );
    }
    if( lhs_has && rhs_has )
    {
      const double diff = fabs( lhs.activityUncertainty.value() - rhs.activityUncertainty.value() );
      const double max_val = std::max(fabs(lhs.activityUncertainty.value()), fabs(rhs.activityUncertainty.value()));
      const double abs_tol = 1.0E-12;
      const double rel_tol = 1.0E-6 * max_val;
      char buffer[512];
      snprintf( buffer, sizeof(buffer),
                "SourceFitDef LHS activityUncertainty (%.17g) != RHS activityUncertainty (%.17g): "
                "diff=%.17g, abs_tol=%.17g, rel_tol=%.17g, max_val=%.17g",
                lhs.activityUncertainty.value(), rhs.activityUncertainty.value(),
                diff, abs_tol, rel_tol, max_val );
      throw runtime_error( buffer );
    }
    throw runtime_error( "SourceFitDef LHS activityUncertainty != RHS activityUncertainty (both unset)" );
  }
  
  if( lhs.ageUncertainty )
    assert( lhs.ageUncertainty.value() > 0.0 );
  if( rhs.ageUncertainty )
    assert( rhs.ageUncertainty.value() > 0.0 );
  
  if( !close_enough( lhs.ageUncertainty, rhs.ageUncertainty ) )
    throw runtime_error( "SourceFitDef LHS ageUncertainty != RHS ageUncertainty" );
  
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  if( !close_enough(lhs.truthActivity, rhs.truthActivity) )
    throw runtime_error( "SourceFitDef LHS truthActivity != RHS truthActivity" );
  if( !close_enough(lhs.truthActivityTolerance, rhs.truthActivityTolerance) )
    throw runtime_error( "SourceFitDef LHS truthActivityTolerance != RHS truthActivityTolerance" );
  if( !close_enough(lhs.truthAge, rhs.truthAge) )
    throw runtime_error( "SourceFitDef LHS truthAge != RHS truthAge" );
  if( !close_enough(lhs.truthAgeTolerance, rhs.truthAgeTolerance) )
    throw runtime_error( "SourceFitDef LHS truthAgeTolerance != RHS truthAgeTolerance" );
#endif //#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
}//void equalEnough( const SourceFitDef &lhs, const SourceFitDef &rhs )

#endif

  
TraceSourceInfo::TraceSourceInfo()
  : m_type( GammaInteractionCalc::TraceActivityType::NumTraceActivityType ),
    m_fitActivity( false ),
    m_nuclide( nullptr ),
    m_activity( 0.0 ),
    m_relaxationDistance( 0.0f )
{
}
   

void SourceFitDef::deSerialize( const ::rapidxml::xml_node<char> *src_node )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  
  const rapidxml::xml_node<> *name_node = src_node->first_node( "Name", 4 );
  const rapidxml::xml_node<> *activity_node = src_node->first_node( "Activity", 8 );
  const rapidxml::xml_node<> *self_atten_node = src_node->first_node( "ShieldingDeterminedActivity", 27 ); //XML version 1.0
  const rapidxml::xml_node<> *src_type_node = src_node->first_node( "SourceType", 10 ); //XML version 1.1
  
  const char *self_atten_value = self_atten_node ? self_atten_node->value() : nullptr;
  const char *src_type_value = src_type_node ? src_type_node->value() : nullptr;
  
  if( !name_node || !name_node->value() || !activity_node
     || (!self_atten_value && !src_type_value)  )
    throw runtime_error( "Missing necessary element for sources XML" );
  
  const rapidxml::xml_node<> *activity_value_node = activity_node->first_node( "Value", 5 );
  
  if( !activity_value_node )
    throw runtime_error( "No activity value node" );
  const rapidxml::xml_attribute<char> *fit_activity_attr = activity_value_node->first_attribute( "Fit", 3 );
  
  const rapidxml::xml_node<> *age_node = src_node->first_node( "Age", 3 );
  if( !age_node )
    throw runtime_error( "Missing necessary age element for sources XML" );
  
  const rapidxml::xml_node<> *age_value_node = age_node->first_node( "Value", 5 );
  if( !age_value_node )
    throw runtime_error( "Missing <Value> child of <Age> in sources XML" );

  const rapidxml::xml_attribute<char> *fit_age_attr = age_value_node->first_attribute( "Fit", 3 );
  const rapidxml::xml_attribute<char> *age_defining_attr = age_value_node->first_attribute( "AgeDefiningNuclide", 18 );
  if( !age_defining_attr )
    age_defining_attr = age_value_node->first_attribute( "AgeMaster", 9 ); //sm_xmlSerializationVersion


  if( !activity_value_node->value()
      || !age_value_node->value()
      || !fit_activity_attr || !fit_activity_attr->value()
      || !fit_age_attr || !fit_age_attr->value() )
    throw runtime_error( "Missing/invalid node for sources XML" );
  
  SourceFitDef &row = *this;
  row.nuclide = db->nuclide( name_node->value() );
  if( !row.nuclide )
    throw runtime_error( "Invalid nuclide for sources XML" );
  if( !(stringstream(activity_value_node->value()) >> row.activity) )
    throw runtime_error( "Failed to read activity" );
  
  if( !(stringstream(age_value_node->value()) >> row.age) )
    throw runtime_error( "Failed to read age" );
  if( !(stringstream(fit_activity_attr->value()) >> row.fitActivity) )
    throw runtime_error( "Failed to read fit_act" );
  if( !(stringstream(fit_age_attr->value()) >> row.fitAge) )
    throw runtime_error( "Failed to read fit_age" );
  
  if( self_atten_value )
  {
    // Depreciated as of XML version 0.1
    bool selfAtten = false;
    if( !(stringstream(self_atten_value) >> selfAtten) )
      throw runtime_error( "Failed to read shieldingIsSource" );
    row.sourceType = selfAtten ? ShieldingSourceFitCalc::ModelSourceType::Intrinsic : ShieldingSourceFitCalc::ModelSourceType::Point;
  }//if( self_atten_value )
  
  if( src_type_value )
  {
    // XML version >= 0.1
    if( SpecUtils::iequals_ascii(src_type_value, "Point") )
      row.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
    else if( SpecUtils::iequals_ascii(src_type_value, "Intrinsic") )
      row.sourceType = ShieldingSourceFitCalc::ModelSourceType::Intrinsic;
    else if( SpecUtils::iequals_ascii(src_type_value, "Trace") )
      row.sourceType = ShieldingSourceFitCalc::ModelSourceType::Trace;
    else
      throw runtime_error( "Invalid value of 'SourceType'" );
  }//if( src_type_value )

  if( !age_defining_attr || !age_defining_attr->value() )
    row.ageDefiningNuc = nullptr;
  else
    row.ageDefiningNuc = db->nuclide( age_defining_attr->value() );
  
  
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  auto getTruth = [src_node]( const char *truthName, const bool isActivity,
                              boost::optional<double> &value ){
    value.reset();
    
    auto node = src_node->first_node( truthName );
    if( !node )
      return;
    
    try
    {
      if( isActivity )
        value = PhysicalUnits::stringToActivity( node->value() );
      else
        value = PhysicalUnits::stringToTimeDuration( node->value() );
      
      //cout << "Set '" << truthName << "' to value " << *value << " from '" << node->value() <<  "' while deserializing" << endl;
    }catch(...)
    {
      cerr << "Failed to read back in " << truthName << " from " << node->value() << endl;
    }
  };//getTruth(...)
  
  getTruth( "TruthActivity", true, row.truthActivity );
  getTruth( "TruthActivityTolerance", true, row.truthActivityTolerance );
  getTruth( "TruthAge", false, row.truthAge );
  getTruth( "TruthAgeTolerance", false, row.truthAgeTolerance );
#endif
  
  row.activityUncertainty.reset();
  if( const rapidxml::xml_node<> *activity_uncert_node = activity_node->first_node( "Uncertainty", 11 ) )
  {
    double uncert = 0.0;
    if( activity_uncert_node->value() && (stringstream(activity_uncert_node->value()) >> uncert) && (uncert > 0.0) )
      row.activityUncertainty = uncert;
  }
  
  row.ageUncertainty.reset();
  if( const rapidxml::xml_node<> *age_uncert_node = age_node->first_node( "Uncertainty", 11 ) )
  {
    double uncert = 0.0;
    if( age_uncert_node->value() && (stringstream(age_uncert_node->value()) >> uncert) && (uncert > 0.0) )
      row.ageUncertainty = uncert;
  }
}//void SourceFitDef::deSerialize( const ::rapidxml::xml_node<char> *parent_node )
  
  
::rapidxml::xml_node<char> *SourceFitDef::serialize( rapidxml::xml_node<char> *isotope_nodes ) const
{
  if( !isotope_nodes )
    throw runtime_error( "SourceFitDef::serialize: invalid parent" );
  
  rapidxml::xml_document<char> *doc = isotope_nodes->document();
  if( !doc )
    throw runtime_error( "SourceFitDef::serialize: couldnt get document" );
  
  
  rapidxml::xml_node<char> *nuclide_node = doc->allocate_node( rapidxml::node_element, "Nuclide" );
  isotope_nodes->append_node( nuclide_node );
  
  const char *value = doc->allocate_string( nuclide->symbol.c_str() );
  rapidxml::xml_node<char> *name_node = doc->allocate_node( rapidxml::node_element, "Name", value );
  nuclide_node->append_node( name_node );
  
  value = "";
  const char *srctype = "";
  switch( sourceType )
  {
    case ShieldingSourceFitCalc::ModelSourceType::Point:
      srctype = "Point";
      value = "0";
      break;
      
    case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
      srctype = "Intrinsic";
      value = "1";
      break;
      
    case ShieldingSourceFitCalc::ModelSourceType::Trace:
      srctype = "Trace";
      value = "0";
      break;
  }//switch( sourceType )
  
  // <ShieldingDeterminedActivity> is depreciated as of XML version 0.1, but leaving in to be
  //   backward compatible with version 0.0 (trace sources will be treated as point sources)
  rapidxml::xml_node<char> *determined_note = doc->allocate_node( rapidxml::node_element, "ShieldingDeterminedActivity", value );
  nuclide_node->append_node( determined_note );
  
  // I dont think we actually need to note the source type here, because the ShieldingSelects
  //  already have this info, but when we deSerialize, we'll note
  //  SourceFitModel::SourceFitDef::sourceType - which is maybe not the best because it creates
  //  a condition for the XML to become inconsistent with itself
  value = srctype;
  determined_note = doc->allocate_node( rapidxml::node_element, "SourceType", value );
  nuclide_node->append_node( determined_note );
  
  rapidxml::xml_node<char> *activity_node = doc->allocate_node( rapidxml::node_element, "Activity" );
  nuclide_node->append_node( activity_node );
  
  value = doc->allocate_string( std::to_string(activity).c_str() );
  rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, "Value", value );
  activity_node->append_node( node );
  
  value = fitActivity ? "1" : "0";
  rapidxml::xml_attribute<char> *attr = doc->allocate_attribute( "Fit", value );
  node->append_attribute( attr );
  
  
  rapidxml::xml_node<char> *age_node = doc->allocate_node( rapidxml::node_element, "Age" );
  nuclide_node->append_node( age_node );
  
  value = doc->allocate_string( std::to_string(age).c_str() );
  node = doc->allocate_node( rapidxml::node_element, "Value", value );
  age_node->append_node( node );
  
  value = fitAge ? "1" : "0";
  attr = doc->allocate_attribute( "Fit", value );
  node->append_attribute( attr );
  
  if( ageDefiningNuc && (ageDefiningNuc != nuclide) )
  {
    value = ageDefiningNuc->symbol.c_str();
    
    if( SourceFitDef::sm_xmlSerializationMajorVersion == 0 )
    {
      //Depreciating tag 20201201, will remove when XML serialization is updated
      attr = doc->allocate_attribute( "AgeMaster", value );
      node->append_attribute( attr );
    }
    
    attr = doc->allocate_attribute( "AgeDefiningNuclide", value );
    node->append_attribute( attr );
  }//if( ageDefiningNuc && (ageDefiningNuc != nuclide) )
  
  
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  auto addTruth = [doc,nuclide_node]( const char *truthName, const bool isActivity,
                                     const boost::optional<double> &value ){
    if( !value )
      return;
    string strval;
    const bool useCuries = true;
    if( isActivity )
      strval = PhysicalUnits::printToBestActivityUnits( *value, 6, useCuries );
    else
      strval = PhysicalUnits::printToBestTimeUnits( *value, 6 );
    
    const char *txtvalue = doc->allocate_string( strval.c_str() );
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, truthName, txtvalue );
    nuclide_node->append_node( node );
  };//addTruth(...)
  
  addTruth( "TruthActivity", true, truthActivity );
  addTruth( "TruthActivityTolerance", true, truthActivityTolerance );
  addTruth( "TruthAge", false, truthAge );
  addTruth( "TruthAgeTolerance", false, truthAgeTolerance );
#endif
  
  const bool ageIsFittable = !PeakDef::ageFitNotAllowed( nuclide );
  const char *deprecated_comment = doc->allocate_string( "AgeIsFittable is deprecated; retained for backward compatibility" );
  nuclide_node->append_node( doc->allocate_node( rapidxml::node_comment, nullptr, deprecated_comment ) );

  const char *fit_value = doc->allocate_string( std::to_string(ageIsFittable ? 1 : 0).c_str() );
  rapidxml::xml_node<char> *age_fit_node = doc->allocate_node( rapidxml::node_element, "AgeIsFittable", fit_value );
  nuclide_node->append_node( age_fit_node );

  if( activityUncertainty && (*activityUncertainty > 0.0) )
  {
    const double uncert = *activityUncertainty;
    assert( uncert > 0.0 );
    rapidxml::xml_node<char> *activity_node = nuclide_node->first_node("Activity");
    if( !activity_node )
      throw runtime_error( "SourceFitDef::serialize: missing Activity node" );
    const char *uncert_str = doc->allocate_string( std::to_string(uncert).c_str() );
    rapidxml::xml_node<char> *uncert_node = doc->allocate_node( rapidxml::node_element, "Uncertainty", uncert_str );
    activity_node->append_node( uncert_node );
  }

  if( ageUncertainty && (*ageUncertainty > 0.0) )
  {
    const double uncert = *ageUncertainty;
    assert( uncert > 0.0 );
    rapidxml::xml_node<char> *age_node = nuclide_node->first_node("Age");
    if( !age_node )
      throw runtime_error( "SourceFitDef::serialize: missing Age node" );
    const char *uncert_str = doc->allocate_string( std::to_string(uncert).c_str() );
    rapidxml::xml_node<char> *uncert_node = doc->allocate_node( rapidxml::node_element, "Uncertainty", uncert_str );
    age_node->append_node( uncert_node );
  }

  return nuclide_node;
}//::rapidxml::xml_node<char> *SourceFitDef::serialize( rapidxml::xml_node<char> *parent_node )
  
  
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
  const double max_act = max(fabs(lhs.m_activity), fabs(rhs.m_activity));
  const double rel_tolerance = 1.0E-6 * max_act;
  const double abs_tolerance = 1.0E-12;
  if( (act_diff > abs_tolerance) && (act_diff > rel_tolerance) )
  {
    char buffer[512];
    snprintf( buffer, sizeof(buffer),
              "TraceSourceInfo LHS activity (%.17g) != RHS activity (%.17g): "
              "diff=%.17g, abs_tol=%.17g, rel_tol=%.17g, max_act=%.17g",
              lhs.m_activity, rhs.m_activity, act_diff, abs_tolerance, rel_tolerance, max_act );
    throw runtime_error( buffer );
  }
  
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
    m_nuclideFractions_{},
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
      
  #if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
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
      
      
      if( m_forFitting && !m_nuclideFractions_.empty() )
      {
        // For backward compatibility with version 0.1, we'
        bool vestigial_FitMassFraction = false;
        for( const auto &el_nucs : m_nuclideFractions_ )
        {
          for( const std::tuple<const SandiaDecay::Nuclide *,double,bool> &nuc : el_nucs.second )
            vestigial_FitMassFraction |= std::get<2>(nuc);
        }
        
        if( vestigial_FitMassFraction )
        {
          name = "FitMassFraction";
          value = "1";
          mass_frac_node = doc->allocate_node( rapidxml::node_element, name, value );
          material_node->append_node( mass_frac_node );
          
          auto attr = doc->allocate_attribute( "remark", "compatibility with v1.0.12 and before." );
          mass_frac_node->append_attribute( attr );
        }
      }//if( m_forFitting )
      
      for( const auto &el_nucs : m_nuclideFractions_ )
      {
        const SandiaDecay::Element *el = el_nucs.first;
        assert( el );
        
        for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_to_frac : el_nucs.second )
        {
          const SandiaDecay::Nuclide * const nuc = std::get<0>(nuc_to_frac);
          const float frac = std::get<1>(nuc_to_frac);
          const bool fit = std::get<2>(nuc_to_frac);
          
          iso_node = doc->allocate_node( rapidxml::node_element, "Nuclide" );
          material_node->append_node( iso_node );
          
          // If the "other" non-source fraction, we will rely on recognizing its a element, and not
          //  nuclide when we decode it.
          value = doc->allocate_string( nuc ? nuc->symbol.c_str() : el->symbol.c_str() );
          node = doc->allocate_node( rapidxml::node_element, "Name", value );
          iso_node->append_node( node );
          
          const string fracstr = SpecUtils::printCompact(frac, 12);
          value = doc->allocate_string( fracstr.c_str(), fracstr.size() + 1 );
          node = doc->allocate_node( rapidxml::node_element, "MassFrac", value );
          iso_node->append_node( node );
          
          value = fit ? "1" : "0";
          node = doc->allocate_node( rapidxml::node_element, "Fit", value );
          iso_node->append_node( node );
        }//for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_to_frac : el_nucs.second )
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
      
      if( !m_truthFitMassFractions.empty() )
      {
        char buffer[32] = { '\0' };
        
        name = "TruthMassFractions";
        rapidxml::xml_node<char> *mass_frac_node = doc->allocate_node( rapidxml::node_element, name );
        material_node->append_node( mass_frac_node );
        
        for( const auto &nuc_val_tol : m_truthFitMassFractions )
        {
          const SandiaDecay::Element * const el = nuc_val_tol.first;
          const map<const SandiaDecay::Nuclide *,pair<double,double>> &nucs_vals = nuc_val_tol.second;
          
          assert( el );
          if( !el )
            continue;
          
          for( const auto &nuc_val : nucs_vals )
          {
            const SandiaDecay::Nuclide * const nuc = nuc_val.first;
            const pair<double,double> &value_tolerance = nuc_val.second;
            
            name = "SelfAttenSrc";
            rapidxml::xml_node<char> *src_node = doc->allocate_node( rapidxml::node_element, name );
            mass_frac_node->append_node( src_node );
            
            name = "Nuclide";
            value = doc->allocate_string( nuc ? nuc->symbol.c_str() : el->symbol.c_str() );
            node = doc->allocate_node( rapidxml::node_element, name, value );
            src_node->append_node( node );
            
            name = "Value";
            snprintf( buffer, sizeof(buffer), "%.9g", value_tolerance.first );
            value = doc->allocate_string( buffer );
            node = doc->allocate_node( rapidxml::node_element, name, value );
            src_node->append_node( node );
            
            name = "Tolerance";
            snprintf( buffer, sizeof(buffer), "%.9g", value_tolerance.second );
            value = doc->allocate_string( buffer );
            node = doc->allocate_node( rapidxml::node_element, name, value );
            src_node->append_node( node );
          }//for( const auto &nuc_val : nucs_vals )
        }//for( const auto &nuc_val_tol : m_truthFitMassFractions )
      }//if( !m_truthFitMassFractions.empty() )
  #endif
    }//if( m_isGenericMaterial ) / else
    
    return base_node;
  }//void serialize( rapidxml::xml_node<char> *parent_node ) const
  
  
  void ShieldingInfo::deSerialize( const rapidxml::xml_node<char> *shield_node )
  {
    using GammaInteractionCalc::GeometryType;
    
    rapidxml::xml_attribute<char> *attr;
    rapidxml::xml_node<char> *node, *geom_node, *generic_node, *material_node;
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
   
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
      m_nuclideFractions_.clear();
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
      
  #if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
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
      if( !MaterialDB::instance() )
        throw runtime_error( "De-serializing shielding requires valid Material DB" );
        
      m_isGenericMaterial = false;
      m_nuclideFractions_.clear();
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
        const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
        std::shared_ptr<const Material> mat;

        try
        {
          mat = matdb->material( material_name );
        }catch( std::exception & )
        {
        }

        if( !mat )
        {
          // Maybe the user specified a chemical formula, like "U0.99Np0.01"
          try
          {
            const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
            mat = MaterialDB::materialFromChemicalFormula( material_name, db );
          }catch( std::exception & )
          {
          }
        }//if( !mat )

        if( !mat )
          throw runtime_error( "Invalid shielding material: '" + material_name + "'" );
        m_material = mat;
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

#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
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
      
      m_truthFitMassFractions.clear();
      const rapidxml::xml_node<> *mass_frac_node = XML_FIRST_NODE(material_node, "TruthMassFractions");
      XML_FOREACH_CHILD( src_node, mass_frac_node, "SelfAttenSrc" )
      {
        const string nuc_str = SpecUtils::xml_value_str( XML_FIRST_NODE(src_node, "Nuclide") );
        const string val_str = SpecUtils::xml_value_str( XML_FIRST_NODE(src_node, "Value") );
        const string tol_str = SpecUtils::xml_value_str( XML_FIRST_NODE(src_node, "Tolerance") );
        
        assert( !nuc_str.empty() && !val_str.empty() && !tol_str.empty() );
        if( nuc_str.empty() || val_str.empty() || tol_str.empty() )
          continue;
        
        const SandiaDecay::Nuclide * const nuc = db->nuclide( nuc_str );
        const SandiaDecay::Element * el = nullptr;
        
        if( nuc )
          el = db->element( nuc->atomicNumber );
        else
          el = db->element( nuc_str );
        
        assert( el );
        if( !el )
          continue;
        
        double value, tolerance;
        if( !(stringstream(val_str) >> value) || !(stringstream(tol_str) >> tolerance)
           || (value < 0.0) || (value > 1.0) || (tolerance < 0.0) || (tolerance > 1.0) )
        {
          assert( 0 );
          continue;
        }
        
        m_truthFitMassFractions[el][nuc] = make_pair( value, tolerance );
      }//XML_FOREACH_CHILD( src_node, mass_frac_node, "SelfAttenSrc" )
#endif //#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
      
      bool vestigual_FitMassFraction = false;
      if( m_forFitting )
      {
        const rapidxml::xml_node<> *fitmassfrac_node = XML_FIRST_NODE(material_node, "FitMassFraction");
        if( fitmassfrac_node && fitmassfrac_node->value() )
        {
          const string valstr = SpecUtils::xml_value_str( fitmassfrac_node );
          if( !(stringstream(valstr) >> vestigual_FitMassFraction) )
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
        const rapidxml::xml_node<> * const fit_node = XML_FIRST_NODE( iso_node, "Fit" );
        
        const string namestr = SpecUtils::xml_value_str( name_node );
        const string fractionstr = SpecUtils::xml_value_str( frac_node );
        if( namestr.empty() || fractionstr.empty() )
          throw runtime_error( "Missing invalid name/mass frac node form iso" );
        
        const SandiaDecay::Nuclide * const nuc = db->nuclide( namestr );
        const SandiaDecay::Element * const el = nuc ? db->element(nuc->atomicNumber) : db->element(namestr);
        if( !el )
          throw runtime_error( namestr + " is not a valid isotope or element" );
        
        double fraction;
        if( !(stringstream(fractionstr) >> fraction) )
          throw runtime_error( "Invalid mass fraction: '" + fractionstr + "'" );
        
        for( const auto &v : m_nuclideFractions_[el] )
        {
          if( std::get<0>(v) == nuc )
            throw runtime_error( "Nuclide '" + namestr + "' specified multiple times." );
        }
        
        bool fitFraction = vestigual_FitMassFraction;
        if( fit_node )
        {
          const string valstr = SpecUtils::xml_value_str( fit_node );
          if( !(stringstream(valstr) >> fitFraction) )
            throw runtime_error( "Invalid FitMassFraction value '" + valstr + "' for " + namestr );
        }//if( fit_node )
        
        if( !m_material )
          throw runtime_error( "Mass-fraction specified, without a valid shielding material" );
        
        bool valid_element = false;
        for( size_t i = 0; !valid_element && (i < m_material->elements.size()); ++i )
          valid_element = (m_material->elements[i].first == el);
        
        // I dont recall if m_material->nuclides are a subset of elements that are in
        //  m_material->elements so we'll check, just to be sure, for the moment
        for( size_t i = 0; !valid_element && (i < m_material->nuclides.size()); ++i )
          valid_element = (m_material->nuclides[i].first == nuc);
        
        if( !valid_element )
          throw runtime_error( "A mass-fraction for '" + nuc->symbol
                              + "' was specified, but material '" + m_material->name
                              + "' doesnt have that element." );
        
        m_nuclideFractions_[el].emplace_back( nuc, fraction, fitFraction );
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
    // "FD1": "FD2": fit the corresponding dimensions
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
      const double ad_g_cm2 = m_dimensions[1] * PhysicalUnits::cm2 / PhysicalUnits::gram;
      answer += "&AD=" + SpecUtils::printCompact(ad_g_cm2, 7);
      answer += "&AN=" + SpecUtils::printCompact(m_dimensions[0], 7);
      
      if( m_forFitting && m_fitDimensions[0] )
        answer += "&FAN=1";
      if( m_forFitting && m_fitDimensions[1] )
        answer += "&FAD=1";
    }else
    {
      std::string material = m_material ? m_material->name : string();
      SpecUtils::ireplace_all(material, "#", "%23" );
      SpecUtils::ireplace_all(material, "&", "%26" );
      const string::size_type open_pos = material.find('(');
      if( open_pos != string::npos )
      {
        const string::size_type close_pos = material.find(')', open_pos);
        if( close_pos != string::npos )
          material.erase(open_pos, close_pos - open_pos + 1);
      }
      SpecUtils::trim( material );
      
      answer += "&N=" + material;
      
      using PhysicalUnits::cm;
        
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          answer += "&G=S&D1=" + SpecUtils::printCompact(m_dimensions[0]/cm, 7) + "cm";
          if( m_forFitting && m_fitDimensions[0] )
            answer += "&FD1=1";
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          if( m_geometry == GammaInteractionCalc::GeometryType::CylinderEndOn )
            answer += "&G=CE";
          else
            answer += "&G=CS";
          answer += "&D1=" + SpecUtils::printCompact(m_dimensions[0]/cm, 7) + "cm";
          answer += "&D2=" + SpecUtils::printCompact(m_dimensions[1]/cm, 7) + "cm";
          if( m_forFitting && m_fitDimensions[0] )
            answer += "&FD1=1";
          if( m_forFitting && m_fitDimensions[1] )
            answer += "&FD2=1";
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          answer += "&G=R";
          answer += "&D1=" + SpecUtils::printCompact(m_dimensions[0]/cm, 7) + "cm";
          answer += "&D2=" + SpecUtils::printCompact(m_dimensions[1]/cm, 7) + "cm";
          answer += "&D3=" + SpecUtils::printCompact(m_dimensions[2]/cm, 7) + "cm";
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
    
  
  void ShieldingInfo::handleAppUrl( std::string query_str )
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
        m_dimensions[1] = ad * (PhysicalUnits::gram / PhysicalUnits::cm2);
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
        const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();
        if( !matdb )
          throw runtime_error( "ShieldingInfo::handleAppUrl: MaterialDB not available" );

        std::shared_ptr<const Material> mat = matdb->material( materialstr );

        if( !mat )
        {
          try
          {
            const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
            mat = MaterialDB::materialFromChemicalFormula( materialstr, db );
          }catch( std::exception & )
          {
          }
        }//if( !mat )

        if( !mat )
          throw runtime_error( "Invalid material name '" + materialstr + "'" );
        m_material = mat;
      }//if( !materialstr.empty() ) / else
      
      
      // Self-atten source stuff
      //std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractions;
      
      // Trace-source stuff
      //std::vector<TraceSourceInfo> m_traceSources;
    }//if( generic ) / else
  }//void ShieldingInfo::handleAppUrl( std::string query_str )
    
  
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
  void ShieldingInfo::equalEnough( const ShieldingInfo &lhs, const ShieldingInfo &rhs )
  {
    if( lhs.m_geometry != rhs.m_geometry )
      throw runtime_error( "ShieldingInfo LHS Geometry != RHS Geometry ("
                          + string(GammaInteractionCalc::to_str(lhs.m_geometry))
                          + " != " + string(GammaInteractionCalc::to_str(rhs.m_geometry)) + ")" ) ;
    
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
        {
          throw runtime_error( "ShieldingInfo for dimension " + std::to_string(i) + " differs:"
                               + " LHS=" + std::to_string(lhs.m_dimensions[i])
                               +  ", RHS=" + std::to_string(rhs.m_dimensions[i]) );
        }
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
    if( lhs.m_nuclideFractions_.size() != rhs.m_nuclideFractions_.size() )
      throw runtime_error( "ShieldingInfo LHS NumSourceElements (" + std::to_string(lhs.m_nuclideFractions_.size()) + ") != RHS NumSourceElements (" + std::to_string(rhs.m_nuclideFractions_.size()) + ")" );
    
    for( const auto &elToNucFrac : lhs.m_nuclideFractions_ )
    {
      const SandiaDecay::Element * const el = elToNucFrac.first;
      assert( el );
      if( !rhs.m_nuclideFractions_.count(el) )
        throw runtime_error( "ShieldingInfo LHS has element "
                            + (el ? el->symbol : string("nullptr"))
                            + " as source, but RHS doesnt" );
      if( !el )
        throw logic_error( "Got nullptr SandiaDecay::Element in ShieldingInfo::equalEnough" );
      
      const auto lhs_el_pos = lhs.m_nuclideFractions_.find(el);
      const auto rhs_el_pos = rhs.m_nuclideFractions_.find(el);
      
      if( lhs_el_pos == end(lhs.m_nuclideFractions_) )
        throw logic_error( "Failed to find SandiaDecay::Element in lhs" );
      
      if( rhs_el_pos == end(rhs.m_nuclideFractions_) )
        throw logic_error( "RHS m_nuclideFractions_ does not have " + (el ? el->symbol : "nullptr") + " element." );
      
      const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> &lhs_nucs = lhs_el_pos->second;
      const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> &rhs_nucs = rhs_el_pos->second;
      
      if( lhs_nucs.size() != rhs_nucs.size() )
        throw runtime_error( "ShieldingInfo LHS number of nucs (" + std::to_string(lhs_nucs.size()) 
                            + ") != RHS number of nucs (" + std::to_string(rhs_nucs.size()) 
                            + ") for " + el->symbol );
      
      for( const tuple<const SandiaDecay::Nuclide *,double,bool> &lhs_nuc : lhs_nucs )
      {
        const SandiaDecay::Nuclide * const nuc = std::get<0>( lhs_nuc );
        const double lval = std::get<1>( lhs_nuc );
        const double lfit = std::get<2>( lhs_nuc );
        
        const tuple<const SandiaDecay::Nuclide *,double,bool> *rhs_nuc = nullptr;
        for( const tuple<const SandiaDecay::Nuclide *,double,bool> &rhs_nuc_try : rhs_nucs )
        {
          if( nuc == std::get<0>(rhs_nuc_try) )
          {
            rhs_nuc = &rhs_nuc_try;
            break;
          }
        }
        
        if( !rhs_nuc )
          throw runtime_error( "ShieldingInfo LHS has nuclide "
                              + (nuc ? nuc->symbol : string("nullptr"))
                              + " as source, but RHS doesnt" );
        
        const double rval = std::get<1>( *rhs_nuc );
        const double rfit = std::get<2>( *rhs_nuc );
        
        const float mval = std::max( fabs(lval), fabs(rval) );
        const float dval = fabs(lval - rval);
        
        if( (mval > 1.0E-9) && (fabs(dval) > 1.0E-6*mval) )
          throw runtime_error( "ShieldingInfo LHS MassFraction for "
                              + (nuc ? nuc->symbol : string("other"))
                              + " != RHS MassFraction ("
                              + std::to_string(lval) + " vs " + std::to_string(rval) + ")" );
        
        if( lfit != rfit )
          throw runtime_error( "ShieldingInfo LHS Fit MassFraction for "
                              + (nuc ? nuc->symbol : string("nullptr"))
                              + " != RHS Fit MassFraction ("
                              + std::to_string(lfit) + " vs " + std::to_string(rfit) + ")" );
      }//for( const tuple<const SandiaDecay::Nuclide *,double,bool> &rhs_nuc : lhs_nucs )
    }//for( const auto &nucFrac : lmfrac )
    
    // Trace-source stuff
    if( lhs.m_traceSources.size() != rhs.m_traceSources.size() )
      throw runtime_error( "ShieldingInfo LHS NumTraceSources != RHS NumTraceSources" );
    
    for( size_t i = 0; i < lhs.m_traceSources.size(); ++i )
      TraceSourceInfo::equalEnough( lhs.m_traceSources[i], rhs.m_traceSources[i] );
  }//void equalEnough( const ShieldingInfo &lhs, const ShieldingInfo &rhs )
#endif //PERFORM_DEVELOPER_CHECKS
  
  
FitShieldingInfo::FitShieldingInfo()
  : ShieldingInfo(),
  m_dimensionUncerts{ -1.0, -1.0, -1.0 },
  m_nuclideFractionUncerts{},
  m_traceSourceActivityUncerts{}
{
}
  
  
void ShieldingSourceFitOptions::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  // For backwards compatibility, we wont create a <ShieldingSourceFitOptions> node,
  //  but just put values directly under `parent_node`.
  if( !parent_node )
    throw runtime_error( "ShieldingSourceFitOptions::serialize: invalid parent" );
  
  rapidxml::xml_document<char> *doc = parent_node->document();
  if( !doc )
    throw runtime_error( "ShieldingSourceFitOptions::serialize: couldnt get document" );
  
  const char *name = "MultipleIsotopesPerPeak";
  const char *value = multiple_nucs_contribute_to_peaks ? "1" : "0";
  rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  name = "AttenuateForAir";
  value = attenuate_for_air ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  name = "DecayCorrect";
  value = account_for_decay_during_meas ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  name = "MultithreadComputation";
  value = multithread_self_atten ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  
  name = "BackgroundPeakSubtraction";
  value = background_peak_subtract ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  name = "SameAgeIsotopes";
  value = same_age_isotopes ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  name = "ComputeEffectiveShielding";
  value = compute_effective_shielding ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
  
  name = "PhotopeakClusterSigma";
  char buffer[64] = { '\0' };
  snprintf( buffer, sizeof(buffer), "%.9g", photopeak_cluster_sigma );
  value = doc->allocate_string( buffer );
  node = doc->allocate_node( rapidxml::node_element, name, value );
  parent_node->append_node( node );
}//void serialize( rapidxml::xml_node<char> *parent_node ) const
  
  
void ShieldingSourceFitOptions::deSerialize( const rapidxml::xml_node<char> *parent_node )
{
  if( !parent_node )
    throw runtime_error( "ShieldingSourceFitOptions::deSerialize: invalid parent" );
  
  // Reset values
  *this = ShieldingSourceFitOptions();
  
  auto boolval = []( const rapidxml::xml_node<char> *n ) -> bool {
    const string val = SpecUtils::xml_value_str( n );
    if( val.empty() )
      throw runtime_error( "ShieldingSourceFitOptions node " + SpecUtils::xml_name_str(n) + " missing value" );
    return (val == "1") || SpecUtils::iequals_ascii(val, "true")  || SpecUtils::iequals_ascii(val, "yes");
  };
  
  const rapidxml::xml_node<char> *node = XML_FIRST_NODE( parent_node, "MultipleIsotopesPerPeak" );
  if( !node )
    throw runtime_error( "ShieldingSourceFitOptions missing MultipleIsotopesPerPeak node" );
  
  multiple_nucs_contribute_to_peaks = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "AttenuateForAir" );
  if( node )
    attenuate_for_air = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "DecayCorrect" );
  if( node )
    account_for_decay_during_meas = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "MultithreadComputation" );
  if( node )
    multithread_self_atten = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "BackgroundPeakSubtraction" );
  if( node )
    background_peak_subtract = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "SameAgeIsotopes" );
  if( node )
    same_age_isotopes = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "ComputeEffectiveShielding" );
  if( node )
    compute_effective_shielding = boolval( node );
  
  node = XML_FIRST_NODE( parent_node, "PhotopeakClusterSigma" );
  if( node )
  {
    const string clusterstr = SpecUtils::xml_value_str(node);
    if( !(stringstream(clusterstr) >> photopeak_cluster_sigma) )
      throw runtime_error( "ShieldingSourceFitOptions invalid cluster sigma: '" + clusterstr + "'" );
  }//if( node )
}//void deSerialize( const rapidxml::xml_node<char> *parent_node )
  
  
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
void ShieldingSourceFitOptions::equalEnough( const ShieldingSourceFitOptions &lhs, const ShieldingSourceFitOptions &rhs )
{
  if( lhs.multiple_nucs_contribute_to_peaks != rhs.multiple_nucs_contribute_to_peaks )
    throw runtime_error( "ShieldingSourceFitOptions LHS multiple_nucs_contribute_to_peaks != RHS multiple_nucs_contribute_to_peaks" );
  
  if( lhs.attenuate_for_air != rhs.attenuate_for_air )
    throw runtime_error( "ShieldingSourceFitOptions LHS attenuate_for_air != RHS attenuate_for_air" );
  
  if( lhs.account_for_decay_during_meas != rhs.account_for_decay_during_meas )
    throw runtime_error( "ShieldingSourceFitOptions LHS account_for_decay_during_meas != RHS account_for_decay_during_meas" );
  
  if( lhs.multithread_self_atten != rhs.multithread_self_atten )
    throw runtime_error( "ShieldingSourceFitOptions LHS multithread_self_atten != RHS multithread_self_atten" );
  
  if( fabs(lhs.photopeak_cluster_sigma - rhs.photopeak_cluster_sigma) > 1.0E-8
     && fabs(lhs.photopeak_cluster_sigma - rhs.photopeak_cluster_sigma) > 1.0E-6*std::max(lhs.photopeak_cluster_sigma,rhs.photopeak_cluster_sigma) )
    throw runtime_error( "ShieldingSourceFitOptions LHS photopeak_cluster_sigma != RHS photopeak_cluster_sigma" );
  
  if( lhs.background_peak_subtract != rhs.background_peak_subtract )
    throw runtime_error( "ShieldingSourceFitOptions LHS background_peak_subtract != RHS background_peak_subtract" );
  
  if( lhs.same_age_isotopes != rhs.same_age_isotopes )
    throw runtime_error( "ShieldingSourceFitOptions LHS same_age_isotopes != RHS same_age_isotopes" );
  
  if( lhs.compute_effective_shielding != rhs.compute_effective_shielding )
    throw runtime_error( "ShieldingSourceFitOptions LHS compute_effective_shielding != RHS compute_effective_shielding" );
}//void equalEnough( const ShieldingSourceFitOptions &lhs, const ShieldingSourceFitOptions &rhs )
#endif
  
  
ModelFitProgress::ModelFitProgress()
  : m_mutex{},
  chi2( std::numeric_limits<double>::max() ),
  elapsedTime( 0.0 ),
  numFcnCalls( 0 ),
  parameters{}
{    
}

  
/** Detects non-fatal warning conditions for a completed (Final) fit and appends human-readable
 messages to `results->warnings`.  Centralizing this means every consumer - the desktop fit-message
 area, the phone fit bar's "issue" pill, and the text/HTML/CSV reports - picks up each warning
 automatically; add future warning conditions here rather than at call sites.

 English strings (no WString::tr), matching the `errormsgs` convention: this runs on the fit worker
 thread / in batch, where the Wt message bundles aren't necessarily available.  Avoid '<'/'>'/'&' in
 the text so it renders correctly in the HTML report and the XHTML message area without escaping.
 */
static void check_for_fit_warnings( ShieldingSourceFitCalc::ModelFitResults &results )
{
  // x-ray peaks: the activity/shielding model doesn't account for x-ray production, so a fit that
  //  relies on x-ray peaks may be unreliable.  (Previously a transient toast at peak-select time.)
  for( const PeakDef &peak : results.foreground_peaks )
  {
    const SandiaDecay::RadParticle * const radpart = peak.decayParticle();
    if( radpart && (radpart->type == SandiaDecay::XrayParticle) )
    {
      results.warnings.push_back( "One or more peaks used in the fit are x-rays. X-ray production is"
        " not modeled by the activity/shielding fit, so results that rely on these peaks may be"
        " unreliable." );
      break;
    }
  }//for( const PeakDef &peak : results.foreground_peaks )

  // Questionable fit: the average per-peak deviation shown on the fit chart is
  //  dev = sqrt(chi2 / num_peaks) - it divides by the NUMBER OF PEAKS, not the degrees of freedom
  //  (see `dev` in ShieldingSourceFitPlot.js).  A good fit is around 1 sigma, so flag dev above 3.
  //  Compute and format it the same way so the warning's number matches the chart's "<dev>".
  const size_t num_peaks = results.foreground_peaks.size();
  if( (num_peaks > 0) && (results.chi2 > 0.0) )
  {
    const double dev = std::sqrt( results.chi2 / static_cast<double>(num_peaks) );
    if( dev > 3.0 )
    {
      char devbuf[32] = { '\0' };
      snprintf( devbuf, sizeof(devbuf), "%.2f", dev );  // match the chart's dev.toFixed(2)
      results.warnings.push_back( "The fit is questionable: the average peak deviation (dev) is "
        + string(devbuf) + " sigma, well above the roughly 1 sigma expected for a good fit - the"
        " model may not describe the data well." );
    }
  }//if( num_peaks > 0 && chi2 > 0 )
}//void check_for_fit_warnings( ModelFitResults & )


/** Fills the source/shielding interpretation of the fit results (fit_src_info,
 final_shieldings, and the calculation-detail logs) from the final parameter
 values and their uncertainties.

 Shared by the Minuit2 and Ceres fit drivers; the caller must hold results->m_mutex.
 */
static void fill_fit_results( std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                              const std::vector<double> &params,
                              const std::vector<double> &errors,
                              const unsigned int ndof,
                              std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results )
{
  {// Begin set fit source info
    results->fit_src_info.clear();
      
    const vector<ShieldingSourceFitCalc::SourceFitDef> &initialSources = chi2Fcn->initialSourceDefinitions();
      
    // Finally we'll set activities and ages
    const size_t nnucs = chi2Fcn->numNuclides();
    results->fit_src_info.resize( nnucs );
      
    assert( initialSources.size() == nnucs );
    if( initialSources.size() != nnucs )
      throw runtime_error( "Somehow different number of initial source and nuclides." );
      
      
    //Go through and set the ages and activities fit for
    for( size_t nucn = 0; nucn < nnucs; ++nucn )
    {
      const SandiaDecay::Nuclide *nuc = chi2Fcn->nuclide( nucn );
        
      size_t initial_row = nnucs;
      for( size_t row = 0; (initial_row == nnucs) && (row < nnucs); ++row )
      {
        if( initialSources[row].nuclide == nuc )
          initial_row = row;
      }//for( size_t row = 0; row < nnucs; ++row )
        
      assert( initial_row != nnucs );
      if( initial_row == nnucs )
        throw runtime_error( "Unable to finish initial source for nuclide "
                            + (nuc ? nuc->symbol : string("null")) );
        
      const ShieldingSourceFitCalc::SourceFitDef &initialdef = initialSources[initial_row];
      ShieldingSourceFitCalc::SourceFitDef &row = results->fit_src_info[initial_row];
        
        
      const double age = chi2Fcn->age( nuc, params );
      const double total_activity = chi2Fcn->totalActivity( nuc, params );
        
      row.nuclide = nuc;
      row.fitAge = initialdef.fitAge;
      row.fitActivity = initialdef.fitActivity;
        
      row.age = age;
      row.activity = total_activity;
      row.ageDefiningNuc = initialdef.ageDefiningNuc;
      row.sourceType = initialdef.sourceType;
        
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
      row.truthActivity = initialdef.truthActivity;
      row.truthActivityTolerance = initialdef.truthActivityTolerance;
      row.truthAge = initialdef.truthAge;
      row.truthAgeTolerance = initialdef.truthAgeTolerance;
#endif
        
      row.activityUncertainty.reset();
      if( row.fitActivity || chi2Fcn->isTraceSource(nuc) || chi2Fcn->isSelfAttenSource(nuc) )
      {
        try
        {
          const double activityUncert = chi2Fcn->totalActivityUncertainty( nuc, params, errors );
          if( activityUncert > FLT_EPSILON )
          {
            assert( activityUncert > 0.0 );
            row.activityUncertainty = activityUncert;
          }
        }catch( std::exception & )
        {
          //We dont seem to get here - but JUC
          assert( 0 );
        }
      }
        
      row.ageUncertainty.reset();
      if( row.fitAge )
      {
        const double ageUncert = chi2Fcn->age( nuc, errors );
        if( ageUncert > FLT_EPSILON )
        {
          assert( ageUncert > 0.0 );
          row.ageUncertainty = ageUncert;
        }
      }else
      {
        assert( (max(row.age, initialdef.age) < 1.0E-6)
               || (fabs(row.age - initialdef.age) < 1.0E-5*max(row.age,initialdef.age)) );
      }
        
    }//for( int ison = 0; ison < niso; ++ison )
  }// End set fit source info
    
    
    
  results->final_shieldings.clear();
  assert( results->initial_shieldings.size() == chi2Fcn->numMaterials() );
  if( results->initial_shieldings.size() != chi2Fcn->numMaterials() )
    throw std::logic_error( "Pre/Post fit number of shieldings not equal" );
    
  for( size_t shielding_index = 0; shielding_index < chi2Fcn->numMaterials(); ++shielding_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &initial_shield = results->initial_shieldings[shielding_index];
      
    ShieldingSourceFitCalc::FitShieldingInfo shield;
    shield.m_forFitting = true;
    shield.m_geometry = chi2Fcn->geometry();
    shield.m_isGenericMaterial = chi2Fcn->isGenericMaterial(shielding_index);
    assert( shield.m_isGenericMaterial == initial_shield.m_isGenericMaterial );
    shield.m_material = initial_shield.m_material;
    assert( shield.m_material == chi2Fcn->material(shielding_index) );
      
#if( INCLUDE_ANALYSIS_TEST_SUITE || PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
    shield.m_truthFitMassFractions = initial_shield.m_truthFitMassFractions;
      
    for( size_t i = 0; i < 3; ++i )
    {
      shield.m_truthDimensions[i] = initial_shield.m_truthDimensions[i];
      shield.m_truthDimensionsTolerances[i] = initial_shield.m_truthDimensionsTolerances[i];
    }
#endif
      
    for( size_t i = 0; i < 3; ++i )
    {
      shield.m_dimensions[i] = 0.0;
      shield.m_fitDimensions[i] = false;
    }
      
    const int shield_start_par = static_cast<int>(2*chi2Fcn->numNuclides() + 3*shielding_index);
      
      
    if( shield.m_isGenericMaterial )
    {
      const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
      const double an = chi2Fcn->atomicNumber( shielding_index, params );
      const double ad = chi2Fcn->arealDensity( shielding_index, params ) / adUnits;
        
      shield.m_dimensions[0] = an;
      shield.m_fitDimensions[0] = initial_shield.m_fitDimensions[0];
      if( shield.m_fitDimensions[0] )
        shield.m_dimensionUncerts[0] = chi2Fcn->atomicNumber( shielding_index, errors );
        
      shield.m_dimensions[1] = ad;
      shield.m_fitDimensions[1] = initial_shield.m_fitDimensions[1];
      if( shield.m_fitDimensions[1] )
        shield.m_dimensionUncerts[1] = chi2Fcn->arealDensity( shielding_index, errors ) / adUnits;
        
      // There looks to be a bug in Minuit that IsFixed() doesnt work
      //assert( shield.m_fitDimensions[0] != fitParams.Parameter(shield_start_par).IsFixed() );
      //assert( shield.m_fitDimensions[1] != fitParams.Parameter(shield_start_par + 1).IsFixed() );
    }else
    {
      const map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,double,bool>>>
      selfAttenSrcInfo = chi2Fcn->selfAttenSrcInfo( shielding_index, params, errors );
        
      const vector<const SandiaDecay::Nuclide *> trace_nucs = chi2Fcn->traceNuclidesForMaterial( shielding_index );
      
      
      for( const auto &el_nucs : selfAttenSrcInfo )
      {
        const SandiaDecay::Element * const el = el_nucs.first;
        assert( el );
          
        double post_fit_frac = 0.0;
        for( const tuple<const SandiaDecay::Nuclide *,double,double,bool> &nuc_info : el_nucs.second )
        {
          const SandiaDecay::Nuclide * const nuc = get<0>(nuc_info);
          const double &mass_frac = get<1>(nuc_info);
          const double &mass_frac_uncert = get<2>(nuc_info);
          const bool fit_frac = get<3>(nuc_info);
            
          assert( fit_frac || (mass_frac_uncert <= 0.0) );
            
          shield.m_nuclideFractions_[el].emplace_back( nuc, mass_frac, fit_frac );
            
          if( fit_frac && (mass_frac_uncert > 0.0) )
            shield.m_nuclideFractionUncerts[el][nuc] = mass_frac_uncert;
            
          if( fit_frac )
            post_fit_frac += mass_frac;
        }//for( loop over nuclides )
          
        double pre_fit_frac = 0.0;
        assert( initial_shield.m_nuclideFractions_.count(el) );
          
        for( const auto &el_infos : initial_shield.m_nuclideFractions_ )
        {
          if( el_infos.first != el )
            continue;
            
          for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_infos.second )
          {
            if( get<2>(nuc_info) )
              pre_fit_frac += get<1>(nuc_info);
          }
        }//for( loop over initial shielding self-atten elements )
          
        const double frac_diff = fabs(pre_fit_frac - post_fit_frac);
        assert( (frac_diff < 1.0E-12) || (frac_diff < 1.0E-5*std::max(pre_fit_frac,post_fit_frac)) );
          
        if( (frac_diff > 1.0E-12) && ((frac_diff/std::max(pre_fit_frac,post_fit_frac)) > 1.0E-5) ) //limits chosen arbitrarily
          throw logic_error( "Mass fraction for self-atten src element " + el->symbol
                            + " should be " + to_string(pre_fit_frac) + " but calculation yielded "
                            + to_string(post_fit_frac) );
      }//for( const auto &el_nucs : selfAttenSrcInfo )
        
      for( const SandiaDecay::Nuclide * const nuc : trace_nucs )
      {
        ShieldingSourceFitCalc::TraceSourceInfo trace;
        trace.m_nuclide = nuc;
        trace.m_type = chi2Fcn->traceSourceActivityType( nuc );
        const int ind = static_cast<int>( chi2Fcn->nuclideIndex( nuc ) );
          
        // There looks to be a bug in Minuit that IsFixed() doesnt work
        //trace.m_fitActivity = !fitParams.Parameter(2*ind).IsFixed();
        bool foundTrace = false;
        for( size_t i = 0; !foundTrace && (i < initial_shield.m_traceSources.size()); ++i )
        {
          foundTrace = (initial_shield.m_traceSources[i].m_nuclide == nuc);
          if( foundTrace )
            trace.m_fitActivity = initial_shield.m_traceSources[i].m_fitActivity;
        }//
          
        trace.m_activity = chi2Fcn->activity( nuc, params );
          
        if( trace.m_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
          trace.m_relaxationDistance = chi2Fcn->relaxationLength( nuc );
          
        shield.m_traceSources.push_back( trace );
          
        if( trace.m_fitActivity )
          shield.m_traceSourceActivityUncerts[nuc] = chi2Fcn->activityUncertainty( nuc, params, errors );
      }//for( const SandiaDecay::Nuclide * const nuc : trace_nucs )
        
      switch( shield.m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          shield.m_dimensions[0] = chi2Fcn->sphericalThickness( shielding_index, params );
          // There looks to be a bug in Minuit that IsFixed() doesnt work
          //shield.m_fitDimensions[0] = !fitParams.Parameter(shield_start_par).IsFixed();
          //assert( shield.m_fitDimensions[0] == initial_shield.m_fitDimensions[0] );
          shield.m_fitDimensions[0] = initial_shield.m_fitDimensions[0];
            
          if( shield.m_fitDimensions[0] )
            shield.m_dimensionUncerts[0] = chi2Fcn->sphericalThickness( shielding_index, errors );
          break;
            
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          shield.m_dimensions[0] = chi2Fcn->cylindricalRadiusThickness( shielding_index, params );
          shield.m_dimensions[1] = chi2Fcn->cylindricalLengthThickness( shielding_index, params );
          // There looks to be a bug in Minuit that IsFixed() doesnt work
          //shield.m_fitDimensions[0] = !fitParams.Parameter(shield_start_par).IsFixed();
          //shield.m_fitDimensions[1] = !fitParams.Parameter(shield_start_par + 1 ).IsFixed();
          //assert( shield.m_fitDimensions[0] == initial_shield.m_fitDimensions[0] );
          //assert( shield.m_fitDimensions[1] == initial_shield.m_fitDimensions[1] );
          shield.m_fitDimensions[0] = initial_shield.m_fitDimensions[0];
          shield.m_fitDimensions[1] = initial_shield.m_fitDimensions[1];
            
          if( shield.m_fitDimensions[0] )
            shield.m_dimensionUncerts[0] = chi2Fcn->cylindricalRadiusThickness( shielding_index, errors );
          if( shield.m_fitDimensions[1] )
            shield.m_dimensionUncerts[1] = chi2Fcn->cylindricalLengthThickness( shielding_index, errors );
          break;
            
        case GammaInteractionCalc::GeometryType::Rectangular:
          shield.m_dimensions[0] = chi2Fcn->rectangularWidthThickness( shielding_index, params );
          shield.m_dimensions[1] = chi2Fcn->rectangularHeightThickness( shielding_index, params );
          shield.m_dimensions[2] = chi2Fcn->rectangularDepthThickness( shielding_index, params );
          // There looks to be a bug in Minuit that IsFixed() doesnt work
          //shield.m_fitDimensions[0] = !fitParams.Parameter(shield_start_par ).IsFixed();
          //shield.m_fitDimensions[1] = !fitParams.Parameter(shield_start_par + 1 ).IsFixed();
          //shield.m_fitDimensions[2] = !fitParams.Parameter(shield_start_par + 2 ).IsFixed();
          //assert( shield.m_fitDimensions[0] == initial_shield.m_fitDimensions[0] );
          //assert( shield.m_fitDimensions[1] == initial_shield.m_fitDimensions[1] );
          //assert( shield.m_fitDimensions[2] == initial_shield.m_fitDimensions[2] );
          shield.m_fitDimensions[0] = initial_shield.m_fitDimensions[0];
          shield.m_fitDimensions[1] = initial_shield.m_fitDimensions[1];
          shield.m_fitDimensions[2] = initial_shield.m_fitDimensions[2];
            
          if( shield.m_fitDimensions[0] )
            shield.m_dimensionUncerts[0] = chi2Fcn->rectangularWidthThickness( shielding_index, errors );
          if( shield.m_fitDimensions[1] )
            shield.m_dimensionUncerts[1] = chi2Fcn->rectangularHeightThickness( shielding_index, errors );
          if( shield.m_fitDimensions[2] )
            shield.m_dimensionUncerts[2] = chi2Fcn->rectangularDepthThickness( shielding_index, errors );
          break;
            
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( shield.m_geometry )
    }//if( shield.m_isGenericMaterial ) / else
      
    results->final_shieldings.push_back( shield );
  }//for( int i = 0; i < nshieldings; ++i )
    
    
    
  {// Begin logging detailed info, we'll later use to template reports
    GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache mixcache;
    auto peak_calc_details = make_unique<vector<GammaInteractionCalc::PeakDetail>>();
    const auto peak_comparisons = chi2Fcn->energy_chi_contributions( params, errors, mixcache,
                                                                    &(results->peak_calc_log),
                                                                    peak_calc_details.get() );
      
    results->peak_calc_details = std::move(peak_calc_details);
    results->peak_comparisons.reset( new vector<GammaInteractionCalc::PeakResultPlotInfo>(peak_comparisons) );
      
      
    if( !results->peak_calc_log.empty() )
    {
      char buffer[64];
      snprintf( buffer, sizeof(buffer), "There %s %i parameter%s fit for",
               (ndof>1 ? "were" : "was"), int(ndof), (ndof>1 ? "s" : "") );
      results->peak_calc_log.push_back( "&nbsp;" );
      results->peak_calc_log.push_back( buffer );
    }//if( !m_calcLog.empty() )
      
    try
    {
      auto shielding_details = make_unique<vector<GammaInteractionCalc::ShieldingDetails>>();
      chi2Fcn->log_shield_info( params, errors, results->fit_src_info, *shielding_details );
        
      // Now fill in a little info we need the results for
      assert( results->initial_shieldings.size() == shielding_details->size() );
      for( size_t i = 0; i < results->initial_shieldings.size(); ++i )
      {
        const ShieldingSourceFitCalc::ShieldingInfo &initial_shield = results->initial_shieldings[i];
        if( i >= shielding_details->size() )
          continue; //wont happen, but jic
        if( i >= results->final_shieldings.size() )
          continue; //wont happen, but jic
          
        const ShieldingSourceFitCalc::FitShieldingInfo &final_shield = results->final_shieldings[i];
          
        GammaInteractionCalc::ShieldingDetails &calc_detail = shielding_details->at( i );
        for( size_t j = 0; j < 3; ++j )
        {
          calc_detail.m_fit_dimension[j] = initial_shield.m_fitDimensions[j];
          if( calc_detail.m_fit_dimension[j] )
            calc_detail.m_dimension_uncert[j] = final_shield.m_dimensionUncerts[j];
        }//
      }//for( size_t i = 0; i < results->initial_shieldings.size(); ++i )
        
        
      results->shield_calc_details = std::move(shielding_details);
    }catch( std::exception &e )
    {
      results->errormsgs.push_back( e.what() );
    }
      
    try
    {
      auto shielding_details = make_unique<vector<GammaInteractionCalc::SourceDetails>>();
      chi2Fcn->log_source_info( params, errors, results->fit_src_info, *shielding_details );
      results->source_calc_details = std::move(shielding_details);
    }catch( std::exception &e )
    {
      results->errormsgs.push_back( e.what() );
    }
      
    if( chi2Fcn->options().compute_effective_shielding )
    {
      try
      {
        GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache eff_mixcache;
        results->effective_shielding = make_unique<vector<GammaInteractionCalc::EffectiveShieldingInfo>>(
                                chi2Fcn->computeEffectiveShielding( params, eff_mixcache ) );
      }catch( std::exception &e )
      {
        results->errormsgs.push_back( "Failed to compute effective shielding: " + string(e.what()) );
      }
    }//if( options().compute_effective_shielding )
    
    // Check if background subtraction is enabled, but no peaks actually background subtracted
    assert( results->peak_calc_details );
    if( results->peak_calc_details && chi2Fcn->options().background_peak_subtract )
    {
      size_t num_back_sub_peaks = 0;
      for( const GammaInteractionCalc::PeakDetail &p : *results->peak_calc_details )
        num_back_sub_peaks += (p.backgroundCounts > 0.0f);
      if( num_back_sub_peaks == 0 )
        results->errormsgs.push_back( "Background peak subtraction requested, but no background"
                                     " peak overlapped a foreground peak.");
    }//if( background peak subtraction selected )
  }// end logging detailed info, we'll later use to template reports

  // Surface non-fatal warnings (poor average deviation, x-ray peaks, ...) on the results so the
  //  GUI, phone fit bar, and reports can all show them.  Only reached on the Final/success path.
  check_for_fit_warnings( *results );
}//fill_fit_results(...)


void fit_model_minuit2( const std::string wtsession,
                         std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                         std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams,
                         std::shared_ptr<ShieldingSourceFitCalc::ModelFitProgress> progress,
                         std::function<void()> progress_fcn,
                         std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results,
                         std::function<void()> finished_fcn )
{
  //The self attenuating probing questions are not tested.
  assert( results );
  
  results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
  results->distance = chi2Fcn->distance();
  results->geometry = chi2Fcn->geometry();
  results->foreground_peaks = chi2Fcn->peaks();
  results->background_peaks = chi2Fcn->backgroundPeaks();
  results->background_normalization_factor = chi2Fcn->backgroundNormalizationFactor();

  results->options = chi2Fcn->options();
  results->initial_shieldings = chi2Fcn->initialShieldings();

  
  shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn::GuiProgressUpdateInfo> gui_progress_info;
  
  if( progress_fcn ) //if we are wanted to post updates to the GUI periodically.
  {
    auto updatefcn = [progress_fcn,wtsession,progress]( size_t ncalls, double elapsed_time,
                                                       double chi2, std::vector<double> pars){
      if( progress )
      {
        std::lock_guard<std::mutex> scoped_lock( progress->m_mutex );
        progress->chi2 = chi2;
        progress->elapsedTime = elapsed_time;
        progress->parameters = pars;
        progress->numFcnCalls = ncalls;
      }
      
      Wt::WApplication *app = Wt::WApplication::instance();
      if( wtsession.empty() || (app && (app->sessionId() == wtsession)) )
        progress_fcn();
      else
        Wt::WServer::instance()->post( wtsession, progress_fcn );
    };//define progressUpdatInfo->m_guiUpdater
    
    //Create a object that will be shared by the Chi2 function, and its callback
    //  which will in turn update the variable shared with the function that
    //  gets posted to the GUI thread.  Its a little convoluted, but I think
    //  this lets us better make a consistent handoff of information (e.g., we
    //  can be sure all the member variables of ModelFitProgress are not-changed
    //  by the time the GUI update function gets executed in the Wt event loop).
    gui_progress_info = std::make_shared<GammaInteractionCalc::ShieldingSourceChi2Fcn::GuiProgressUpdateInfo>( sm_model_update_frequency_ms, updatefcn );
    
    chi2Fcn->setGuiProgressUpdater( gui_progress_info );
  }//if( progress_fcn )
  
  // We will update the GUI with results, unless the status code is
  //  #GammaInteractionCalc::ShieldingSourceChi2Fcn::CalcStatus::CanceledNoUpdate
  bool update_gui = true;
  
  try
  {
    if( !chi2Fcn )
      throw runtime_error( "Programming logic error - Chi2Function pointer is null." );
    
    //I think inputPrams.VariableParameters() == num_fit_params
    if( inputPrams->VariableParameters() > chi2Fcn->peaks().size() )
    {
      
      const string msg = "You are asking to fit "
                         + std::to_string( inputPrams->VariableParameters() )
                         + " parameters, however there are only "
                         + std::to_string( chi2Fcn->peaks().size() )
                         + " peaks, which leads to this being an under-constrained problem";
      throw runtime_error( msg );
    }//if( num_fit_params > peaks.size() )
    
    if( inputPrams->VariableParameters() < 1 )
      throw runtime_error( "No parameters are selected for fitting." );
    
    chi2Fcn->fittingIsStarting( sm_max_model_fit_time_ms );
    
    ROOT::Minuit2::MnUserParameterState inputParamState( *inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
    ROOT::Minuit2::MnMinimize fitter( *chi2Fcn, inputParamState, strategy );
    
    //cout << "Parameters are: {";
    //for( const auto &par : inputPrams->Parameters() )
    //  cout << par.Name() << ", ";
    //cout << endl;
    
    const double tolerance = 2.0*inputPrams->VariableParameters();
    const unsigned int maxFcnCall = 50000;  //default minuit2: 200 + 100 * npar + 5 * npar**2
    
    ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
    
    
    //Try two more times to get a valid fit... a stupid hack
    for( int i = 0; !minimum.IsValid() && i < 2; ++i )
      minimum = fitter( maxFcnCall, tolerance );
    
    ROOT::Minuit2::MnUserParameters fitParams = minimum.UserParameters();
    
    //For some reason the fit to atomic number of generic shielding is horrible!
    //  I'm probably screwing something up in the setup for Minuit, but as a
    //  work around, lets manually coarsely scan AN, and then use the best found
    //  AN to then perform another fit to fine tune it.
    //Not a lot of validation went into this, but it does seem to work better.
    //  For the life of me, I cant get Minuit to fit for AN and AD simultaneously well
    vector<string> fit_generic_an;
    for( auto &par : inputPrams->Parameters() )
    {
      const string &name = par.GetName();
      if( (name.find("Generic_")!=string::npos)
         && (name.find("_AN")!=string::npos)
         && (name.find("_FIXED")==string::npos)
         && !par.IsFixed() )  //Looks to be bug in Minuit, so we need the above "_FIXED"
      {
        fit_generic_an.push_back( name );
      }
    }//for( auto &par : inputPrams->Parameters() )
    
    
    // If we have a self-attenuating or trace source, computation time becomes pretty large, so we
    //  wont do the detailed AN scan in this case
    //  TODO: check if we are fitting anything related to self-attenuating materials dimensions, and
    //        if not, dont reject doing the AN scan.
    if( !fit_generic_an.empty() )
    {
      for( size_t i = 0; i < chi2Fcn->numNuclides(); ++i )
      {
        const SandiaDecay::Nuclide * const nuc = chi2Fcn->nuclide( static_cast<int>(i) );
        if( chi2Fcn->isVolumetricSource(nuc) )
          fit_generic_an.clear();
      }
    }//if( check if we also have a self-attenuating source )
    
    
    if( !fit_generic_an.empty() )  //should add in a check that we arnt doing self attenuating sources.
    {
      //Re-perform the fit scanning in atomic number to find best.
      // TODO: Fix all other AN for the scan, and update inputPrams, so subsequent generic materials will use new AN
      for( const auto &parname : fit_generic_an )
      {
        std::mutex best_an_mutex;
        const double orig_fit_an = minimum.UserParameters().Value(parname);
        const double orig_an_error = minimum.UserParameters().Error(parname);
        const double orig_chi2 = minimum.Fval();
        
        int best_an = static_cast<int>( std::round(orig_fit_an) );
        double best_chi2 = orig_chi2;
        ROOT::Minuit2::FunctionMinimum best_min = minimum;
        
        vector<double> an_chi2s( MassAttenuation::sm_max_xs_atomic_number + 1, std::numeric_limits<double>::max() );
        
        auto calc_for_an = [&an_chi2s,&best_an_mutex,&best_an,&best_chi2,&best_min,inputPrams,&parname,chi2Fcn]( const int an ){
          ROOT::Minuit2::MnUserParameters testpar;
          
          for( const auto &p : inputPrams->Parameters() )
          {
            const string name = p.Name();
            if( name == parname )
            {
              testpar.Add( name, static_cast<double>(an) );
            }else if( p.IsConst() || p.IsFixed() || (name.find("_FIXED") != string::npos) )
            {
              testpar.Add( name, p.Value() );
            }else
            {
              testpar.Add( name, p.Value(), p.Error() );
              if( p.HasLowerLimit() )
                testpar.SetLowerLimit( name, p.LowerLimit() );
              if( p.HasUpperLimit() )
                testpar.SetUpperLimit( name, p.UpperLimit() );
            }
          }//for( const auto &p : inputPrams->Parameters() )
          
          try
          {
            ROOT::Minuit2::MnUserParameterState anInputParam( testpar );
            ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
            ROOT::Minuit2::MnMinimize anfitter( *chi2Fcn, anInputParam, strategy );
            
            const double tolerance = 2.0*inputPrams->VariableParameters();
            const unsigned int maxFcnCall = 50000;  //default minuit2: 200 + 100 * npar + 5 * npar**2
            
            ROOT::Minuit2::FunctionMinimum anminimum = anfitter( maxFcnCall, tolerance );
            for( int i = 0; !anminimum.IsValid() && i < 2; ++i )
              anminimum = anfitter( maxFcnCall, tolerance );
            
            {// begin lock on best_an_mutex
              std::lock_guard<std::mutex> lock( best_an_mutex );
              
              an_chi2s[an] = anminimum.Fval();
              if( an_chi2s[an] < best_chi2 )
              {
                best_an = an;
                best_chi2 = an_chi2s[an];
                best_min = anminimum;
              }
            }// end lock on best_an_mutex
          }catch( std::exception &e )
          {
            cerr << "Unexpected exception scanning in AN for shielding/src fit!" << endl;
          }//try / catch
        };//calc_for_an lambda
        
        
        const bool origMultithread = chi2Fcn->options().multithread_self_atten;
        chi2Fcn->setSelfAttMultiThread( false ); //shouldnt affect anything, but JIC
        
        
        SpecUtilsAsync::ThreadPool pool;
        
        // Note 20250106: prior to this, there was a likely bug where the attenuation coefficient,
        //                mu, was being truncated to an integer, and possibly causing some of
        //                the issue with fitting the correct atomic number, and possibly
        //                contributing to the horribleness below.  However, a few quick checks after
        //                fixing this bug shows finding the correct AN is still difficult because
        //                there are local-minima in the AN-AD plane, and in fact a the scans of AN
        //                should include a few different starting AD's, to avoid false-minima still.
        const int initial_an_skip = 5;
        for( int an = MassAttenuation::sm_min_xs_atomic_number;
            an < MassAttenuation::sm_max_xs_atomic_number;
            an += initial_an_skip )
        {
          pool.post( [&calc_for_an,an](){ calc_for_an(an); } );
        }//for( loop over AN )
        
        pool.join();
        
        int detail_scan_start, detail_scan_end, best_coarse_an;
        {// begin lock on best_an_mutex
          std::lock_guard<std::mutex> lock( best_an_mutex );
          
          const auto min_coarse_chi2_iter = std::min_element( begin(an_chi2s), end(an_chi2s) );
          best_coarse_an = static_cast<int>( min_coarse_chi2_iter - begin(an_chi2s) );
          
          //Now scan best_coarse_an +- (initial_an_skip - 1),
          // (clamp to MassAttenuation::sm_min_xs_atomic_number, MassAttenuation::sm_max_xs_atomic_number)
          detail_scan_start = static_cast<int>( best_coarse_an - (initial_an_skip - 1) );
          detail_scan_end = static_cast<int>( best_coarse_an + initial_an_skip );
          detail_scan_start = std::max( detail_scan_start, MassAttenuation::sm_min_xs_atomic_number );
          detail_scan_end = std::min( detail_scan_end, MassAttenuation::sm_max_xs_atomic_number );
          
          assert( detail_scan_start >= 1 );
          assert( detail_scan_end <= MassAttenuation::sm_max_xs_atomic_number );
          
          if( (best_coarse_an < MassAttenuation::sm_min_xs_atomic_number)
             || (best_coarse_an > MassAttenuation::sm_max_xs_atomic_number) )
          {
            detail_scan_start = detail_scan_end;
          }
        }// end lock on best_an_mutex
        
        
        for( int an = detail_scan_start; an <= detail_scan_end; an += 1 )
        {
          if( an != best_coarse_an )
            pool.post( [&calc_for_an,an](){ calc_for_an(an); } );
        }//for( loop over AN )
        
        pool.join();
        
        {// begin lock on best_an_mutex
          std::lock_guard<std::mutex> lock( best_an_mutex );
          
          best_an = std::min( best_an, MassAttenuation::sm_max_xs_atomic_number );
          best_an = std::max( best_an, MassAttenuation::sm_min_xs_atomic_number );
          
          minimum = best_min;
          fitParams = minimum.UserParameters();
          
          // TODO: explore updating inputPrams incase we are fitting multiple generic materials
          //inputPrams->SetValue( parname, best_an );
          
          // We'll approximate errors on AN.  Extraordinarily rough right now.
          // TODO: better estimate the actual lower_and upper AN; for now overestimating error,
          //       sometimes horribly.
          const double up = chi2Fcn->Up();
          double lower_an = best_an, upper_an = best_an;
          for( ; lower_an > MassAttenuation::sm_min_xs_atomic_number; --lower_an )
          {
            if( (an_chi2s[lower_an] != std::numeric_limits<double>::max())
               && ((an_chi2s[lower_an] - best_chi2) >= up) )
              break;
          }
          
          for( ; upper_an < MassAttenuation::sm_max_xs_atomic_number; ++upper_an )
          {
            if( (an_chi2s[upper_an] != std::numeric_limits<double>::max())
               && ((an_chi2s[upper_an] - best_chi2) >= up) )
              break;
          }
          
          float error = 0.5*(upper_an - lower_an);
          
          // If we are near the limits of atomic number, be a little more conservative since our
          //  naive estimation may not have been able to fully go to the atomic number.
          if( ((best_an + error) >= MassAttenuation::sm_max_xs_atomic_number)
             || ((best_an - error) <= MassAttenuation::sm_min_xs_atomic_number) )
          {
            error = std::max( error, static_cast<float>(upper_an - best_an) );
            error = std::max( error, static_cast<float>(best_an - lower_an) );
          }
          
          
          error = std::max( 1.0f, error );
          
          fitParams.SetError( parname, error );
          
          //cout << "Originally fit AN for " << parname
          //<< " was " << orig_fit_an << " +- " << orig_an_error << " with chi2=" << orig_chi2
          //<< ", but final AN is " << best_an << " +- " << error << " with chi2=" << best_chi2 << endl;
        }// end lock on best_an_mutex
        
        chi2Fcn->setSelfAttMultiThread( origMultithread ); //shouldnt affect anything, but JIC
      }//for( const auto &parname : fit_generic_an )
      
    }//if( !fit_generic_an.empty() )
    
    
    std::lock_guard<std::mutex> lock( results->m_mutex );
    
    if( !minimum.IsValid() )
    {
      string msg = "Fit status is not valid:";
      if( minimum.HasMadePosDefCovar() )
        msg += "<br />&nbsp;&nbsp;- Covariance matrix forced positive-definite";
      if( !minimum.HasAccurateCovar() )
        msg += "<br />&nbsp;&nbsp;- Does not have accurate covariance matrix";
      if( minimum.HasReachedCallLimit() )
        msg += "<br />&nbsp;&nbsp;- Optimization reached call limit.";
      if( !minimum.HasValidCovariance() )
        msg += "<br />&nbsp;&nbsp;- Did not have valid covariance,";
      if( !minimum.HasValidParameters() )
        msg += "<br />&nbsp;&nbsp;- Invalid fit parameters.";
      if( minimum.IsAboveMaxEdm() )
        msg += "<br />&nbsp;&nbsp;- The estimated distance to minimum too large.";
      
      results->errormsgs.push_back( msg );
    }//if( !minimum.IsValid() )
    
    const vector<double> params = fitParams.Params();
    const vector<double> errors = fitParams.Errors();
    const unsigned int num_fit_pars = inputPrams->VariableParameters();
    unsigned int ndof = (num_fit_pars > chi2Fcn->peaks().size()) ? static_cast<unsigned int>(0)
                                                                 : static_cast<unsigned int>(chi2Fcn->peaks().size() - num_fit_pars);
    
    results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final;
    results->paramValues = params;
    results->paramErrors = errors;
    results->edm = minimum.Edm();
    results->num_fcn_calls = minimum.NFcn();
    results->numDOF = ndof;
    results->chi2 = minimum.Fval();  //chi2Fcn->DoEval( results->paramValues );
    
    
    fill_fit_results( chi2Fcn, params, errors, ndof, results );
  }catch( GammaInteractionCalc::ShieldingSourceChi2Fcn::CancelException &e )
  {
    const size_t nFunctionCallsSoFar = gui_progress_info->numFunctionCallsSoFar();
    const double bestChi2SoFar = gui_progress_info->bestChi2SoFar();
    const vector<double> bestParsSoFar = gui_progress_info->bestParametersSoFar();
    
    std::lock_guard<std::mutex> lock( results->m_mutex );
    
    switch( e.m_code )
    {
      case GammaInteractionCalc::ShieldingSourceChi2Fcn::CalcStatus::UserCanceled:
        results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::UserCancelled;
        results->errormsgs.push_back( "User canceled fit." );
        break;
        
      case GammaInteractionCalc::ShieldingSourceChi2Fcn::CalcStatus::Timeout:
        results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::TimedOut;
        results->errormsgs.push_back( "Fit took to long." );
        
        if( gui_progress_info )
        {
          results->edm = -1.0;
          results->num_fcn_calls = static_cast<int>( nFunctionCallsSoFar );
          results->chi2 = bestChi2SoFar;
          results->paramValues = bestParsSoFar;
          results->paramErrors.clear();
        }//if( gui_progress_info )
        break;
        
      case GammaInteractionCalc::ShieldingSourceChi2Fcn::CalcStatus::CanceledNoUpdate:
        update_gui = false;
        break;
        
      default:
        results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
        results->errormsgs.push_back( "Fit not performed: " + string(e.what()) );
        break;
    }//switch( e.m_code )
    
  }catch( exception &e )
  {
    std::lock_guard<std::mutex> lock( results->m_mutex );
    results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
    results->errormsgs.push_back( "Fit not performed: " + string(e.what()) );
  }// try / catch
  
  chi2Fcn->fittingIsFinished();
 
  if( finished_fcn && update_gui )
  {
    Wt::WApplication *app = Wt::WApplication::instance();
    if( wtsession.empty() || (app && (app->sessionId() == wtsession)) )
      finished_fcn();
    else
      Wt::WServer::instance()->post( wtsession, finished_fcn );
  }//if( finished_fcn )
}//void fit_model_minuit2( std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams )


namespace
{
  /** A Minuit2-independent description of one fit parameter.

   The fit interface hands the Ceres driver a ROOT::Minuit2::MnUserParameters (it is the
   projects current lingua franca for the activity/shielding fit setup), but everything
   inside the driver works off these instead, so that Minuit2 can eventually be removed
   from the project without touching the Ceres machinery -
   #to_fit_parameter_defs is the single conversion point.
   */
  struct FitParameterDef
  {
    std::string name;
    double value = 0.0;

    /** Initial step size (what Minuit calls the parameter "error"); zero if not set. */
    double step = 0.0;

    /** True for constant/fixed parameters (including the "_FIXED"-named ones). */
    bool is_const = true;

    bool has_lower = false, has_upper = false;
    double lower = 0.0, upper = 0.0;
  };//struct FitParameterDef


  /** The single place the Minuit2 parameter description is converted for the Ceres driver. */
  std::vector<FitParameterDef> to_fit_parameter_defs( const ROOT::Minuit2::MnUserParameters &input )
  {
    std::vector<FitParameterDef> defs;

    for( const ROOT::Minuit2::MinuitParameter &p : input.Parameters() )
    {
      FitParameterDef def;
      def.name = p.GetName();
      def.value = p.Value();
      def.step = p.Error();
      def.is_const = p.IsConst() || p.IsFixed() || (def.name.find("_FIXED") != std::string::npos);
      def.has_lower = p.HasLowerLimit();
      def.has_upper = p.HasUpperLimit();
      def.lower = def.has_lower ? p.LowerLimit() : 0.0;
      def.upper = def.has_upper ? p.UpperLimit() : 0.0;
      defs.push_back( def );
    }//for( loop over Minuit parameters )

    return defs;
  }//to_fit_parameter_defs(...)


  /** The Minuit2-style transform between the bounded "external" parameters the physics
   uses, and the unbounded "internal" parameters the optimizer sees.

   Ceres' projected trust region has a known failure mode for this problem: a parameter
   that rides into its bound far from the optimum gets pinned there (the projected steps
   go to ~zero norm, and parameter-tolerance "convergence" is declared at a terrible
   chi2).  Minuit2 does not suffer this because it optimizes unbounded internal
   parameters, with x = c + s*sin(u) for two-sided limits - so we do the same; the
   transform is smooth, so automatic differentiation just chains through it.

   Note: the alternative used by RelActCalcAuto - keeping plain Ceres bounds but
   loosening them slightly (1e-5 of the range) past the physical limits - was evaluated
   here too (20260612).  With exact bounds, the u_np_1kg analyst case reliably pins its
   mass-fraction parameters at a bound far from the optimum (chi2 ~1e12 vs the true
   ~159); the loosening does rescue it (the slightly-unphysical region past the bound
   restores a usable gradient), but took 55 residual evaluations vs 24 for this
   transform on that problem, and can return slightly-unphysical values (e.g., mass
   fraction 1+1e-5) that would need clamping on readout.  So the transform is used.
   TODO.md has a note to evaluate this transform approach for RelActCalcAuto.
   */
  struct BoundTransform
  {
    // NOTE: from_par() currently only produces Unbounded (plain Ceres bounds) or TwoSided
    //  (sin-transform) - every fit parameter is added with both limits or neither.  The
    //  LowerOnly/UpperOnly arms are retained for a future one-sided-limit parameter.
    enum class Type : int { Unbounded, TwoSided, LowerOnly, UpperOnly };

    Type type = Type::Unbounded;
    double lo = 0.0, up = 0.0;

    /** When type==Unbounded but the parameter has limits, the limits are
     enforced directly by Ceres bound constraints instead of a transform.
     */
    bool use_plain_bounds = false;

    static BoundTransform from_par( const FitParameterDef &p )
    {
      BoundTransform tf;

      if( !p.has_lower && !p.has_upper )
        return tf;

      tf.lo = p.lower;
      tf.up = p.upper;

      // Only "compact" parameters - whose full range is within a couple orders of
      //  magnitude of their step size (mass fractions in [0,1], atomic number in
      //  [1,98], areal density in [0,400 g/cm2]) - get the sin() transform; those are
      //  the ones the projected trust region pins at a bound far from the optimum.
      //  Parameters with sanity-check limits spanning many orders of magnitude
      //  (thicknesses up to 1 km, ages up to centuries) sit so close to their lower
      //  limit, relative to the range, that the transform is horribly conditioned
      //  there - they keep direct coordinates with plain Ceres bounds.
      const double range = (p.has_lower && p.has_upper) ? (tf.up - tf.lo) : -1.0;
      const bool compact = (range > 0.0) && (p.step > 0.0) && (range <= 100.0*p.step);

      if( compact )
        tf.type = Type::TwoSided;
      else
        tf.use_plain_bounds = true;

      return tf;
    }//from_par(...)

    /** External (physical, bounded) value from internal (unbounded) value. */
    template<typename T>
    T to_external( const T &u ) const
    {
      using namespace std;
      using namespace ceres;

      switch( type )
      {
        case Type::Unbounded: return u;
        case Type::TwoSided:  return 0.5*(lo + up) + 0.5*(up - lo)*sin(u);
        case Type::LowerOnly: return lo - 1.0 + sqrt(u*u + 1.0);
        case Type::UpperOnly: return up + 1.0 - sqrt(u*u + 1.0);
      }
      assert( 0 );
      return u;
    }//to_external(...)

    /** Internal value whose to_external() is (approximately) x; for two-sided limits,
     a starting value exactly at a limit is pulled slightly inside, since the
     transform has zero derivative exactly at the limits.
     */
    double to_internal( const double x ) const
    {
      switch( type )
      {
        case Type::Unbounded:
          return x;

        case Type::TwoSided:
        {
          // Pull a value exactly at (or beyond) a limit just inside it, since the
          //  transform derivative is zero exactly at the limits.
          const double ratio = (2.0*x - lo - up) / (up - lo);
          const double max_ratio = 1.0 - 1.0E-10;
          return std::asin( std::max( -max_ratio, std::min( max_ratio, ratio ) ) );
        }

        case Type::LowerOnly:
        {
          const double shifted = std::max( x - lo + 1.0, 1.0 + 1.0E-7 );
          return std::sqrt( shifted*shifted - 1.0 );
        }

        case Type::UpperOnly:
        {
          const double shifted = std::max( up - x + 1.0, 1.0 + 1.0E-7 );
          return std::sqrt( shifted*shifted - 1.0 );
        }
      }
      assert( 0 );
      return x;
    }//to_internal(...)

    /** dx/du - for propagating internal-parameter uncertainties to external ones. */
    double derivative( const double u ) const
    {
      switch( type )
      {
        case Type::Unbounded: return 1.0;
        case Type::TwoSided:  return 0.5*(up - lo)*std::cos(u);
        case Type::LowerOnly: return u / std::sqrt(u*u + 1.0);
        case Type::UpperOnly: return -u / std::sqrt(u*u + 1.0);
      }
      assert( 0 );
      return 1.0;
    }//derivative(...)
  };//struct BoundTransform


  /** The Ceres cost functor for the activity/shielding fit: residual_i =
   (observed_i - expected_i)/uncert_i for each peak included in the fit, with the
   expected counts computed by the templated
   #GammaInteractionCalc::ShieldingSourceChi2Fcn::expected_peak_counts_imp - so with
   T = ceres::Jet<>, the Jacobian is computed by automatic differentiation (including
   through the volumetric-source integration).

   Only the *varying* parameters are exposed to Ceres (constant/fixed ones are baked
   in from their initial values), so no autodiff lanes are wasted on constants.
   */
  struct ShieldSrcCeresCostFcn
  {
    /** How many parameters get their derivatives computed in a single pass; typical
     fits vary ~3-20 parameters, so one or two passes covers everything.
     */
    static const int sm_auto_diff_stride = 16;

    std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> m_fcn;

    /** Full Minuit-layout parameter vector, supplying the constant parameters. */
    std::vector<double> m_initial_params;

    /** Index into the full parameter vector, for each parameter Ceres is varying. */
    std::vector<size_t> m_variable_indices;

    /** The internal<->external transform for each varied parameter. */
    std::vector<BoundTransform> m_transforms;

    /** Background-subtracted observed counts (and uncert) for each peak included in
     the fit - same peak set and ordering as expected_peak_counts_imp returns.
     */
    std::vector<double> m_observed;
    std::vector<double> m_observed_uncert;

    /** Cache of decay mixtures; Ceres evaluates a single residual block serially, so
     no synchronization is needed (each AN-scan solve has its own functor).
     */
    mutable GammaInteractionCalc::ShieldingSourceChi2Fcn::NucMixtureCache m_mix_cache;

    template<typename T>
    bool operator()( T const *const *parameters, T *residuals ) const
    {
      using GammaInteractionCalc::ShieldingSourceChi2Fcn;

      if( m_fcn->currentCancelStatus() != ShieldingSourceChi2Fcn::CalcStatus::NotCanceled )
        return false;

      try
      {
        std::vector<T> full( m_initial_params.size() );
        for( size_t i = 0; i < m_initial_params.size(); ++i )
          full[i] = T( m_initial_params[i] );
        for( size_t i = 0; i < m_variable_indices.size(); ++i )
          full[m_variable_indices[i]] = m_transforms[i].to_external( parameters[0][i] );

        const std::vector<T> expected = m_fcn->expected_peak_counts_imp( full, m_mix_cache );

        assert( expected.size() == m_observed.size() );
        if( expected.size() != m_observed.size() )
          return false;

        for( size_t i = 0; i < expected.size(); ++i )
          residuals[i] = (m_observed[i] - expected[i]) / m_observed_uncert[i];

        // Report progress (and track best-so-far parameters) on the value-only passes
        if constexpr ( std::is_same_v<T,double> )
        {
          double chi2 = 0.0;
          for( size_t i = 0; i < expected.size(); ++i )
            chi2 += residuals[i]*residuals[i];

          m_fcn->reportCompletedEval( chi2, full );
        }

        return true;
      }catch( std::exception &e )
      {
        if( getenv("INTERSPEC_DEBUG_CERES_FIT") )
          std::cerr << "ShieldSrcCeresCostFcn eval failed: " << e.what() << std::endl;
        return false;
      }
    }//operator()(...)
  };//struct ShieldSrcCeresCostFcn


  /** Aborts the Ceres minimization between iterations if the fit was cancelled or
   timed out (the cost function also polls, for faster response).
   */
  class CheckCancelCallback : public ceres::IterationCallback
  {
  public:
    const GammaInteractionCalc::ShieldingSourceChi2Fcn * const m_fcn;

    explicit CheckCancelCallback( const GammaInteractionCalc::ShieldingSourceChi2Fcn *fcn )
      : m_fcn( fcn )
    {
    }

    virtual ceres::CallbackReturnType operator()( const ceres::IterationSummary & ) override
    {
      typedef GammaInteractionCalc::ShieldingSourceChi2Fcn::CalcStatus CalcStatus;
      return (m_fcn->currentCancelStatus() == CalcStatus::NotCanceled)
              ? ceres::SOLVER_CONTINUE : ceres::SOLVER_ABORT;
    }
  };//class CheckCancelCallback


  /** Performs a single Ceres solve of the activity/shielding problem.

   @param chi2Fcn The chi2 function (only expected_peak_counts_imp and cancel status used).
   @param par_defs Source of the parameter bounds and step sizes.
   @param start_full Full Minuit-layout starting parameter values.
   @param var_indices Which entries of start_full to vary.
   @param observed,observed_uncert Per-included-peak observed counts.
   @param[out] final_full start_full with the varied entries updated to the solution.
   @param[out] var_errors If non-null: 1-sigma uncertainties for the varied parameters,
               from the Ceres covariance (zeros if covariance fails).
   @param[out] cov_errmsg If non-null and the covariance fails, a user-displayable note.
   @param[in,out] num_evals If non-null, incremented by the number of residual evaluations.
   @returns the chi2 (=2*final_cost) at the solution.
   */
  double run_ceres_solve( std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                          const std::vector<FitParameterDef> &par_defs,
                          const std::vector<double> &start_full,
                          const std::vector<size_t> &var_indices,
                          const std::vector<double> &observed,
                          const std::vector<double> &observed_uncert,
                          std::vector<double> &final_full,
                          std::vector<double> *var_errors,
                          std::string *cov_errmsg,
                          size_t *num_evals,
                          std::string *convergence_msg = nullptr )
  {
    const size_t num_vars = var_indices.size();
    assert( num_vars > 0 );

    auto functor = std::make_unique<ShieldSrcCeresCostFcn>();
    functor->m_fcn = chi2Fcn;
    functor->m_initial_params = start_full;
    functor->m_variable_indices = var_indices;
    functor->m_observed = observed;
    functor->m_observed_uncert = observed_uncert;

    // The optimizer works on unbounded internal parameters (Minuit2-style transforms);
    //  the functor maps them to the bounded physical values - see #BoundTransform.
    std::vector<BoundTransform> transforms( num_vars );
    for( size_t i = 0; i < num_vars; ++i )
      transforms[i] = BoundTransform::from_par( par_defs[var_indices[i]] );

    functor->m_transforms = transforms;

    auto cost_function = std::make_unique<ceres::DynamicAutoDiffCostFunction<ShieldSrcCeresCostFcn,
                            ShieldSrcCeresCostFcn::sm_auto_diff_stride>>( functor.release() );
    cost_function->AddParameterBlock( static_cast<int>(num_vars) );
    cost_function->SetNumResiduals( static_cast<int>(observed.size()) );

    std::vector<double> u( num_vars );
    for( size_t i = 0; i < num_vars; ++i )
      u[i] = transforms[i].to_internal( start_full[var_indices[i]] );

    ceres::Problem problem;
    problem.AddResidualBlock( cost_function.release(), nullptr, u.data() );  //Problem takes ownership

    // Parameters that didnt get a transform keep their limits as direct bound constraints
    for( size_t i = 0; i < num_vars; ++i )
    {
      if( !transforms[i].use_plain_bounds )
        continue;

      const FitParameterDef &p = par_defs[var_indices[i]];
      if( p.has_lower )
        problem.SetParameterLowerBound( u.data(), static_cast<int>(i), p.lower );
      if( p.has_upper )
        problem.SetParameterUpperBound( u.data(), static_cast<int>(i), p.upper );
    }//for( loop over varied parameters )

    CheckCancelCallback cancel_callback( chi2Fcn.get() );

    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    options.use_nonmonotonic_steps = true;
    options.max_num_iterations = 1000;
    options.function_tolerance = 1.0E-9;
    options.logging_type = ceres::SILENT;
    options.minimizer_progress_to_stdout = false;
    options.callbacks.push_back( &cancel_callback );

    // Even on the transformed (unbounded) parameters, a solution can stall in the
    //  asymptotically-flat region near a limit (where the transform derivative goes to
    //  zero), so if any parameter ends up essentially at a limit, nudge it inside and
    //  re-solve, keeping the best result.
    ceres::Solver::Summary summary;
    double best_chi2 = std::numeric_limits<double>::max();
    std::vector<double> best_u = u;
    ceres::TerminationType best_termination = ceres::FAILURE;

    const int max_attempts = 4;
    for( int attempt = 0; attempt < max_attempts; ++attempt )
    {
      ceres::Solve( options, &problem, &summary );

      if( num_evals )
        *num_evals += static_cast<size_t>( summary.num_residual_evaluations );

      // For development debugging of fit issues: set INTERSPEC_DEBUG_CERES_FIT=1
      if( getenv("INTERSPEC_DEBUG_CERES_FIT") )
      {
        std::cout << summary.FullReport() << std::endl;
        std::cout << "Attempt " << attempt << " solution x={";
        for( size_t i = 0; i < num_vars; ++i )
          std::cout << transforms[i].to_external( u[i] ) << ", ";
        std::cout << "}" << std::endl;
      }//if( INTERSPEC_DEBUG_CERES_FIT )

      const double this_chi2 = 2.0*summary.final_cost;
      const bool improved = (this_chi2 < best_chi2);
      if( improved )
      {
        best_chi2 = this_chi2;
        best_u = u;
        best_termination = summary.termination_type;
      }

      if( chi2Fcn->currentCancelStatus() != GammaInteractionCalc::ShieldingSourceChi2Fcn::CalcStatus::NotCanceled )
        break;

      // A re-solve that didnt improve on the previous attempt wont get better by nudging again
      if( (attempt > 0) && !improved )
        break;

      bool nudged = false;
      for( size_t i = 0; i < num_vars; ++i )
      {
        const BoundTransform &tf = transforms[i];
        const FitParameterDef &p = par_defs[var_indices[i]];

        // "Pinned at a limit" is judged in external (physical) units, on the scale of
        //  the parameters Minuit step size - a thickness of 22 mm with limits
        //  [0, 1 km] is fine, but a mass fraction of 0.9999999 with limits [0,1] and
        //  a step size of ~0.1 is pinned.
        const double step_scale = (p.step > 0.0) ? p.step
                  : ((tf.type == BoundTransform::Type::TwoSided) ? 0.01*(tf.up - tf.lo) : 0.0);
        if( step_scale <= 0.0 )
          continue;

        const double x_ext = tf.to_external( u[i] );

        if( p.has_lower && ((x_ext - tf.lo) < 1.0E-3*step_scale) )
        {
          double x_new = tf.lo + step_scale;
          if( tf.type == BoundTransform::Type::TwoSided )
            x_new = std::min( x_new, 0.5*(tf.lo + tf.up) );
          u[i] = tf.to_internal( x_new );
          nudged = true;
        }else if( p.has_upper && ((tf.up - x_ext) < 1.0E-3*step_scale) )
        {
          double x_new = tf.up - step_scale;
          if( tf.type == BoundTransform::Type::TwoSided )
            x_new = std::max( x_new, 0.5*(tf.lo + tf.up) );
          u[i] = tf.to_internal( x_new );
          nudged = true;
        }
      }//for( loop over varied parameters )

      if( !nudged )
        break;
    }//for( int attempt = 0; attempt < max_attempts; ++attempt )

    // If no solve attempt got a valid cost, the problem is unevaluatable at/near the
    //  starting parameters - report it rather than returning garbage.
    if( best_chi2 == std::numeric_limits<double>::max() )
      throw std::runtime_error( "The optimizer could not evaluate the model"
                                " (check starting values and limits)." );

    // Surface non-convergence (e.g. hit the iteration cap, or a numeric failure) so the
    //  user gets the same "fit may not be reliable" signal the Minuit path provides.  A
    //  user-canceled solve (USER_FAILURE / USER_SUCCESS) is handled separately and not warned.
    if( convergence_msg
        && (best_termination != ceres::CONVERGENCE)
        && (best_termination != ceres::USER_SUCCESS)
        && (best_termination != ceres::USER_FAILURE) )
    {
      *convergence_msg = "The fit did not fully converge (solver: "
            + std::string( ceres::TerminationTypeToString(best_termination) )
            + ") - results may be unreliable or sit at a parameter limit.";
    }

    u = best_u;

    final_full = start_full;
    for( size_t i = 0; i < num_vars; ++i )
      final_full[var_indices[i]] = transforms[i].to_external( u[i] );

    if( var_errors )
    {
      var_errors->assign( num_vars, 0.0 );

      try
      {
        ceres::Covariance::Options cov_opts;
        cov_opts.algorithm_type = ceres::DENSE_SVD;
        cov_opts.null_space_rank = -1;

        ceres::Covariance covariance( cov_opts );
        std::vector<std::pair<const double *, const double *>> cov_blocks;
        cov_blocks.emplace_back( u.data(), u.data() );

        if( covariance.Compute(cov_blocks, &problem) )
        {
          std::vector<double> cov_mat( num_vars*num_vars, 0.0 );
          covariance.GetCovarianceBlock( u.data(), u.data(), cov_mat.data() );

          // Propagate the internal-parameter uncertainties to the external (physical)
          //  parameters: sigma_x = |dx/du| * sigma_u  (same as Minuit2 does).
          for( size_t i = 0; i < num_vars; ++i )
          {
            const double sigma_u = std::sqrt( std::max( 0.0, cov_mat[i*num_vars + i] ) );
            (*var_errors)[i] = std::fabs( transforms[i].derivative( u[i] ) ) * sigma_u;
          }
        }else if( cov_errmsg )
        {
          *cov_errmsg = "Failed to compute parameter uncertainties (a parameter may be"
                        " at a limit, or the problem nearly degenerate).";
        }
      }catch( std::exception &e )
      {
        if( cov_errmsg )
          *cov_errmsg = "Failed to compute parameter uncertainties: " + std::string(e.what());
      }
    }//if( var_errors )

    return best_chi2;
  }//run_ceres_solve(...)
}//namespace


void fit_model_ceres( const std::string wtsession,
                      std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                      std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams,
                      std::shared_ptr<ShieldingSourceFitCalc::ModelFitProgress> progress,
                      std::function<void()> progress_fcn,
                      std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results,
                      std::function<void()> finished_fcn )
{
  using GammaInteractionCalc::ShieldingSourceChi2Fcn;
  typedef ShieldingSourceChi2Fcn::CalcStatus CalcStatus;

  assert( results );

  results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
  results->distance = chi2Fcn->distance();
  results->geometry = chi2Fcn->geometry();
  results->foreground_peaks = chi2Fcn->peaks();
  results->background_peaks = chi2Fcn->backgroundPeaks();
  results->background_normalization_factor = chi2Fcn->backgroundNormalizationFactor();

  results->options = chi2Fcn->options();
  results->initial_shieldings = chi2Fcn->initialShieldings();

  shared_ptr<ShieldingSourceChi2Fcn::GuiProgressUpdateInfo> gui_progress_info;

  if( progress_fcn ) //if we are wanted to post updates to the GUI periodically.
  {
    auto updatefcn = [progress_fcn,wtsession,progress]( size_t ncalls, double elapsed_time,
                                                       double chi2, std::vector<double> pars){
      if( progress )
      {
        std::lock_guard<std::mutex> scoped_lock( progress->m_mutex );
        progress->chi2 = chi2;
        progress->elapsedTime = elapsed_time;
        progress->parameters = pars;
        progress->numFcnCalls = ncalls;
      }

      Wt::WApplication *app = Wt::WApplication::instance();
      if( wtsession.empty() || (app && (app->sessionId() == wtsession)) )
        progress_fcn();
      else
        Wt::WServer::instance()->post( wtsession, progress_fcn );
    };//updatefcn

    gui_progress_info = std::make_shared<ShieldingSourceChi2Fcn::GuiProgressUpdateInfo>(
                                                    sm_model_update_frequency_ms, updatefcn );
    chi2Fcn->setGuiProgressUpdater( gui_progress_info );
  }//if( progress_fcn )

  // We will update the GUI with results, unless the status code is CanceledNoUpdate
  bool update_gui = true;

  try
  {
    if( !chi2Fcn )
      throw runtime_error( "Programming logic error - Chi2Function pointer is null." );

    // Convert to the Minuit2-independent parameter description the driver works off of,
    //  and decompose into constant values and varied indices.  All further logic uses these,
    //  so to_fit_parameter_defs() is the only Minuit2 dependency in this driver.
    const vector<FitParameterDef> par_defs = to_fit_parameter_defs( *inputPrams );

    vector<double> initial_full( par_defs.size(), 0.0 );
    vector<size_t> variable_indices;
    for( size_t i = 0; i < par_defs.size(); ++i )
    {
      initial_full[i] = par_defs[i].value;
      if( !par_defs[i].is_const )
        variable_indices.push_back( i );
    }//for( loop over parameters )

    const size_t num_var_params = variable_indices.size();
    if( num_var_params > chi2Fcn->peaks().size() )
    {
      const string msg = "You are asking to fit "
                         + std::to_string( num_var_params )
                         + " parameters, however there are only "
                         + std::to_string( chi2Fcn->peaks().size() )
                         + " peaks, which leads to this being an under-constrained problem";
      throw runtime_error( msg );
    }//if( num_fit_params > peaks.size() )

    if( num_var_params < 1 )
      throw runtime_error( "No parameters are selected for fitting." );

    chi2Fcn->fittingIsStarting( sm_max_model_fit_time_ms );

    // The background-subtracted observed counts for each peak in the fit - same peak
    //  inclusion and background-subtraction as #ShieldingSourceChi2Fcn::expected_observed_chis
    vector<double> observed, observed_uncert;
    bool warned_zero_uncert = false;
    const vector<PeakDef> &fit_peaks = chi2Fcn->peaks();
    const vector<PeakDef> &back_peaks = chi2Fcn->backgroundPeaks();
    for( const PeakDef &peak : fit_peaks )
    {
      if( !peak.decayParticle() && (peak.sourceGammaType() != PeakDef::AnnihilationGamma) )
        continue;

      double observed_counts = peak.peakArea();
      double observed_uncertainty = peak.peakAreaUncert();
      double backCounts = 0.0, backUncert2 = 0.0;
      const double nsigmaNear = 1.0;

      for( const PeakDef &backPeak : back_peaks )
      {
        const double sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();
        if( fabs(backPeak.mean() - peak.mean()) < (nsigmaNear*sigma) )
        {
          backCounts += backPeak.peakArea();
          backUncert2 += backPeak.peakAreaUncert()*backPeak.peakAreaUncert();
        }
      }//for( const PeakDef &backPeak : back_peaks )

      if( backCounts > 0.0 )
      {
        observed_counts -= backCounts;
        observed_uncertainty = sqrt( observed_uncertainty*observed_uncertainty + backUncert2 );
      }

      // A non-positive area uncertainty would make the residual (obs-exp)/uncert blow up to
      //  inf/NaN and the solver would fail opaquely.  Real fit peaks always have a positive
      //  Poisson area uncertainty; clamp to a Poisson floor so the fit can proceed, and note it.
      if( !(observed_uncertainty > 0.0) )
      {
        observed_uncertainty = std::sqrt( std::max( std::fabs(observed_counts), 1.0 ) );
        if( results && !warned_zero_uncert )
        {
          warned_zero_uncert = true;
          results->errormsgs.push_back( "A peak had a non-positive area uncertainty; using a"
                                        " Poisson estimate so the fit could proceed." );
        }
      }//if( non-positive area uncertainty )

      observed.push_back( observed_counts );
      observed_uncert.push_back( observed_uncertainty );
    }//for( const PeakDef &peak : fit_peaks )

    if( observed.empty() )
      throw runtime_error( "No peaks suitable for fitting." );

    // Figure out the generic-material atomic numbers being fit - if there are any (and
    //  no volumetric sources), the chi2 is multi-modal in AN, so they get a scan.
    vector<size_t> fit_generic_an_indices;
    for( size_t i = 0; i < par_defs.size(); ++i )
    {
      const FitParameterDef &p = par_defs[i];
      if( !p.is_const
          && (p.name.find("Generic_") != string::npos)
          && (p.name.find("_AN") != string::npos) )
      {
        fit_generic_an_indices.push_back( i );
      }
    }//for( loop over parameters )

    // With a self-attenuating or trace source, computation time becomes pretty large,
    //  so skip the detailed AN scan in that case (same policy as the Minuit2 driver).
    if( !fit_generic_an_indices.empty() )
    {
      for( size_t i = 0; i < chi2Fcn->numNuclides(); ++i )
      {
        if( chi2Fcn->isVolumetricSource( chi2Fcn->nuclide(static_cast<int>(i)) ) )
          fit_generic_an_indices.clear();
      }
    }//if( check if we also have a volumetric source )

    size_t num_evals = 0;
    string cov_errmsg, conv_msg;
    vector<double> fit_full, var_errors;

    // Per-parameter error overrides from the AN scans (full-vector index -> error)
    std::map<size_t,double> an_scan_errors;

    // The varied set used for the final (covariance-bearing) solve; scanned ANs are
    //  excluded since their value comes from the scan and their error from the
    //  chi2-vs-AN profile.
    vector<size_t> final_var_indices = variable_indices;

    if( fit_generic_an_indices.empty() )
    {
      const double chi2 = run_ceres_solve( chi2Fcn, par_defs, initial_full, variable_indices,
                                           observed, observed_uncert, fit_full, &var_errors,
                                           &cov_errmsg, &num_evals, &conv_msg );
      results->chi2 = chi2;
    }else
    {
      // First a fit with everything (including AN) free, as the starting reference
      vector<double> free_fit_full;
      double best_chi2 = run_ceres_solve( chi2Fcn, par_defs, initial_full, variable_indices,
                                          observed, observed_uncert, free_fit_full, nullptr,
                                          nullptr, &num_evals );
      vector<double> best_full = free_fit_full;

      for( const size_t an_index : fit_generic_an_indices )
      {
        std::mutex best_an_mutex;
        vector<double> an_chi2s( MassAttenuation::sm_max_xs_atomic_number + 1,
                                 std::numeric_limits<double>::max() );

        // The varied parameters for the scan solves: everything except this AN
        vector<size_t> scan_var_indices;
        for( const size_t index : variable_indices )
        {
          if( index != an_index )
            scan_var_indices.push_back( index );
        }

        if( scan_var_indices.empty() )
          continue;  //fitting only the AN - the scan itself is the fit, but this shouldnt really happen

        std::atomic<size_t> scan_num_evals( 0 );

        auto calc_for_an = [&]( const int an ){
          try
          {
            vector<double> start_full = initial_full;
            start_full[an_index] = static_cast<double>( an );

            vector<double> an_fit_full;
            size_t this_num_evals = 0;
            const double chi2 = run_ceres_solve( chi2Fcn, par_defs, start_full,
                                                 scan_var_indices, observed, observed_uncert,
                                                 an_fit_full, nullptr, nullptr, &this_num_evals );
            scan_num_evals += this_num_evals;

            std::lock_guard<std::mutex> lock( best_an_mutex );
            an_chi2s[an] = chi2;
            if( chi2 < best_chi2 )
            {
              best_chi2 = chi2;
              best_full = an_fit_full;
            }
          }catch( std::exception & )
          {
            cerr << "Unexpected exception scanning in AN for shielding/src fit!" << endl;
          }
        };//calc_for_an lambda

        const bool origMultithread = chi2Fcn->options().multithread_self_atten;
        chi2Fcn->setSelfAttMultiThread( false ); //shouldnt affect anything, but JIC

        SpecUtilsAsync::ThreadPool pool;

        const int initial_an_skip = 5;
        for( int an = MassAttenuation::sm_min_xs_atomic_number;
             an < MassAttenuation::sm_max_xs_atomic_number;
             an += initial_an_skip )
        {
          pool.post( [&calc_for_an,an](){ calc_for_an(an); } );
        }
        pool.join();

        int detail_scan_start, detail_scan_end, best_coarse_an;
        {// begin lock on best_an_mutex
          std::lock_guard<std::mutex> lock( best_an_mutex );

          const auto min_coarse_chi2_iter = std::min_element( begin(an_chi2s), end(an_chi2s) );
          best_coarse_an = static_cast<int>( min_coarse_chi2_iter - begin(an_chi2s) );

          detail_scan_start = best_coarse_an - (initial_an_skip - 1);
          detail_scan_end = best_coarse_an + initial_an_skip;
          detail_scan_start = std::max( detail_scan_start, MassAttenuation::sm_min_xs_atomic_number );
          detail_scan_end = std::min( detail_scan_end, MassAttenuation::sm_max_xs_atomic_number );
        }// end lock on best_an_mutex

        for( int an = detail_scan_start; an <= detail_scan_end; an += 1 )
        {
          if( an != best_coarse_an )
            pool.post( [&calc_for_an,an](){ calc_for_an(an); } );
        }
        pool.join();

        chi2Fcn->setSelfAttMultiThread( origMultithread );

        num_evals += scan_num_evals.load();

        {// begin lock on best_an_mutex
          std::lock_guard<std::mutex> lock( best_an_mutex );

          // Estimate the AN uncertainty by walking the chi2-vs-AN profile out by Up()=1
          //  (extraordinarily rough - same approach as the Minuit2 driver)
          const double up = chi2Fcn->Up();
          const int best_an = static_cast<int>( std::round( best_full[an_index] ) );

          int lower_an = best_an, upper_an = best_an;
          for( ; lower_an > MassAttenuation::sm_min_xs_atomic_number; --lower_an )
          {
            if( (an_chi2s[lower_an] != std::numeric_limits<double>::max())
                && ((an_chi2s[lower_an] - best_chi2) >= up) )
              break;
          }

          for( ; upper_an < MassAttenuation::sm_max_xs_atomic_number; ++upper_an )
          {
            if( (an_chi2s[upper_an] != std::numeric_limits<double>::max())
                && ((an_chi2s[upper_an] - best_chi2) >= up) )
              break;
          }

          double error = 0.5*(upper_an - lower_an);

          if( ((best_an + error) >= MassAttenuation::sm_max_xs_atomic_number)
              || ((best_an - error) <= MassAttenuation::sm_min_xs_atomic_number) )
          {
            error = std::max( error, static_cast<double>(upper_an - best_an) );
            error = std::max( error, static_cast<double>(best_an - lower_an) );
          }

          an_scan_errors[an_index] = std::max( 1.0, error );
        }// end lock on best_an_mutex

        // The scanned AN keeps its scan value; drop it from the final varied set
        final_var_indices.erase( std::remove( begin(final_var_indices), end(final_var_indices), an_index ),
                                 end(final_var_indices) );
      }//for( const size_t an_index : fit_generic_an_indices )

      // Final solve from the best point found (cheap - already at/near the minimum),
      //  to get the covariance for the non-scanned parameters
      if( final_var_indices.empty() )
      {
        fit_full = best_full;
        var_errors.clear();
        results->chi2 = best_chi2;
      }else
      {
        const double chi2 = run_ceres_solve( chi2Fcn, par_defs, best_full, final_var_indices,
                                             observed, observed_uncert, fit_full, &var_errors,
                                             &cov_errmsg, &num_evals, &conv_msg );
        // Report the cost of the parameters we are actually returning (fit_full); the final
        //  solve restarts from best_full so this is normally == best_chi2, but using it keeps
        //  results->chi2 consistent with the reported parameters and per-peak residuals.
        results->chi2 = chi2;
      }
    }//if( fit_generic_an_indices.empty() ) / else

    // If the fit was cancelled or timed out, the solves returned early - reuse the
    //  CancelException handling for reporting.
    const CalcStatus cancel_status = chi2Fcn->currentCancelStatus();
    if( cancel_status != CalcStatus::NotCanceled )
      throw ShieldingSourceChi2Fcn::CancelException( cancel_status );

    // Assemble the full-layout error vector: covariance errors for the varied
    //  parameters, scan-profile errors for scanned ANs, zero for constants.
    vector<double> errors( initial_full.size(), 0.0 );
    if( var_errors.size() == final_var_indices.size() )
    {
      for( size_t i = 0; i < final_var_indices.size(); ++i )
        errors[final_var_indices[i]] = var_errors[i];
    }
    for( const auto &index_error : an_scan_errors )
      errors[index_error.first] = index_error.second;

    std::lock_guard<std::mutex> lock( results->m_mutex );

    if( !cov_errmsg.empty() )
      results->errormsgs.push_back( cov_errmsg );
    if( !conv_msg.empty() )
      results->errormsgs.push_back( conv_msg );

    const unsigned int num_fit_pars = static_cast<unsigned int>( num_var_params );
    const unsigned int ndof = (num_fit_pars > chi2Fcn->peaks().size())
                  ? static_cast<unsigned int>(0)
                  : static_cast<unsigned int>(chi2Fcn->peaks().size() - num_fit_pars);

    results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final;
    results->paramValues = fit_full;
    results->paramErrors = errors;
    results->edm = -1.0;  //no Minuit-style estimated-distance-to-minimum from Ceres
    results->num_fcn_calls = static_cast<int>( num_evals );
    results->numDOF = ndof;

    fill_fit_results( chi2Fcn, fit_full, errors, ndof, results );
  }catch( GammaInteractionCalc::ShieldingSourceChi2Fcn::CancelException &e )
  {
    std::lock_guard<std::mutex> lock( results->m_mutex );

    switch( e.m_code )
    {
      case CalcStatus::UserCanceled:
        results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::UserCancelled;
        results->errormsgs.push_back( "User canceled fit." );
        break;

      case CalcStatus::Timeout:
        results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::TimedOut;
        results->errormsgs.push_back( "Fit took to long." );

        if( gui_progress_info )
        {
          results->edm = -1.0;
          results->num_fcn_calls = static_cast<int>( gui_progress_info->numFunctionCallsSoFar() );
          results->chi2 = gui_progress_info->bestChi2SoFar();
          results->paramValues = gui_progress_info->bestParametersSoFar();
          results->paramErrors.clear();
        }//if( gui_progress_info )
        break;

      case CalcStatus::CanceledNoUpdate:
        update_gui = false;
        break;

      default:
        results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
        results->errormsgs.push_back( "Fit not performed: " + string(e.what()) );
        break;
    }//switch( e.m_code )
  }catch( exception &e )
  {
    std::lock_guard<std::mutex> lock( results->m_mutex );
    results->successful = ShieldingSourceFitCalc::ModelFitResults::FitStatus::InvalidOther;
    results->errormsgs.push_back( "Fit not performed: " + string(e.what()) );
  }// try / catch

  chi2Fcn->fittingIsFinished();

  if( finished_fcn && update_gui )
  {
    Wt::WApplication *app = Wt::WApplication::instance();
    if( wtsession.empty() || (app && (app->sessionId() == wtsession)) )
      finished_fcn();
    else
      Wt::WServer::instance()->post( wtsession, finished_fcn );
  }//if( finished_fcn )
}//void fit_model_ceres(...)


void fit_model( const std::string wtsession,
                std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> chi2Fcn,
                std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams,
                std::shared_ptr<ShieldingSourceFitCalc::ModelFitProgress> progress,
                std::function<void()> progress_fcn,
                std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> results,
                std::function<void()> finished_fcn )
{
#if( USE_CERES_FOR_ACTIVITY_FIT )
  fit_model_ceres( wtsession, chi2Fcn, inputPrams, progress, progress_fcn, results, finished_fcn );
#else
  fit_model_minuit2( wtsession, chi2Fcn, inputPrams, progress, progress_fcn, results, finished_fcn );
#endif
}//void fit_model(...)

}//namespace ShieldingSourceFitCalc
