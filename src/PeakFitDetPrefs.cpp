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
#include <cstdio>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

const int PeakFitDetPrefs::sm_xmlSerializationMajorVersion = 0;
const int PeakFitDetPrefs::sm_xmlSerializationMinorVersion = 0;


PeakFitDetPrefs::PeakFitDetPrefs()
  : m_det_type( PeakFitUtils::CoarseResolutionType::Unknown )
  , m_peak_skew_type( PeakDef::NoSkew )
  , m_lower_energy_skew{}
  , m_upper_energy_skew{}
  , m_roi_independent_skew( false )
  , m_fwhm_method( FwhmMethod::Normal )
  , m_source( LoadingSource::Default )
{
}


bool PeakFitDetPrefs::operator==( const PeakFitDetPrefs &rhs ) const
{
  if( m_det_type != rhs.m_det_type )
    return false;
  if( m_peak_skew_type != rhs.m_peak_skew_type )
    return false;
  if( m_source != rhs.m_source )
    return false;
  if( m_roi_independent_skew != rhs.m_roi_independent_skew )
    return false;
  if( m_fwhm_method != rhs.m_fwhm_method )
    return false;

  for( int i = 0; i < 4; ++i )
  {
    if( m_lower_energy_skew[i] != rhs.m_lower_energy_skew[i] )
      return false;
    if( m_upper_energy_skew[i] != rhs.m_upper_energy_skew[i] )
      return false;
  }

  return true;
}//operator==


bool PeakFitDetPrefs::operator!=( const PeakFitDetPrefs &rhs ) const
{
  return !( *this == rhs );
}//operator!=


const char *PeakFitDetPrefs::to_str( const PeakFitUtils::CoarseResolutionType type )
{
  switch( type )
  {
    case PeakFitUtils::CoarseResolutionType::Low:     return "Low";
    case PeakFitUtils::CoarseResolutionType::LaBr:    return "LaBr";
    case PeakFitUtils::CoarseResolutionType::CZT:     return "CZT";
    case PeakFitUtils::CoarseResolutionType::High:    return "High";
    case PeakFitUtils::CoarseResolutionType::Unknown:  return "Unknown";
  }//switch( type )

  assert( 0 );
  return "Unknown";
}//to_str( CoarseResolutionType )


PeakFitUtils::CoarseResolutionType PeakFitDetPrefs::coarse_res_from_str( const string &str )
{
  if( SpecUtils::iequals_ascii( str, "Low" ) || SpecUtils::iequals_ascii( str, "NaI" ) )
    return PeakFitUtils::CoarseResolutionType::Low;
  if( SpecUtils::iequals_ascii( str, "LaBr" ) || SpecUtils::iequals_ascii( str, "LaBr3" ) )
    return PeakFitUtils::CoarseResolutionType::LaBr;
  if( SpecUtils::iequals_ascii( str, "CZT" ) )
    return PeakFitUtils::CoarseResolutionType::CZT;
  if( SpecUtils::iequals_ascii( str, "High" ) || SpecUtils::iequals_ascii( str, "HPGe" ) )
    return PeakFitUtils::CoarseResolutionType::High;
  if( SpecUtils::iequals_ascii( str, "Unknown" ) )
    return PeakFitUtils::CoarseResolutionType::Unknown;

  throw runtime_error( "PeakFitDetPrefs::coarse_res_from_str: unrecognized value '" + str + "'" );
}//coarse_res_from_str


const char *PeakFitDetPrefs::to_str( const FwhmMethod method )
{
  switch( method )
  {
    case FwhmMethod::Normal:        return "Normal";
    case FwhmMethod::DetFwhm:       return "DetFwhm";
    case FwhmMethod::DetPlusRefine: return "DetPlusRefine";
  }//switch( method )

  assert( 0 );
  return "Normal";
}//to_str( FwhmMethod )


PeakFitDetPrefs::FwhmMethod PeakFitDetPrefs::fwhm_method_from_str( const string &str )
{
  if( SpecUtils::iequals_ascii( str, "Normal" ) )
    return FwhmMethod::Normal;
  if( SpecUtils::iequals_ascii( str, "DetFwhm" ) )
    return FwhmMethod::DetFwhm;
  if( SpecUtils::iequals_ascii( str, "DetPlusRefine" ) )
    return FwhmMethod::DetPlusRefine;

  throw runtime_error( "PeakFitDetPrefs::fwhm_method_from_str: unrecognized value '" + str + "'" );
}//fwhm_method_from_str


const char *PeakFitDetPrefs::to_str( const LoadingSource src )
{
  switch( src )
  {
    case LoadingSource::Default:                  return "Default";
    case LoadingSource::UserInputInGui:           return "UserInputInGui";
    case LoadingSource::FromDetectorPeakResponse: return "FromDetectorPeakResponse";
    case LoadingSource::DefaultForDetectorType:   return "DefaultForDetectorType";
    case LoadingSource::FromSpectralData:         return "FromSpectralData";
  }//switch( src )

  assert( 0 );
  return "Default";
}//to_str( LoadingSource )


PeakFitDetPrefs::LoadingSource PeakFitDetPrefs::loading_source_from_str( const string &str )
{
  if( SpecUtils::iequals_ascii( str, "Default" ) )
    return LoadingSource::Default;
  if( SpecUtils::iequals_ascii( str, "UserInputInGui" ) )
    return LoadingSource::UserInputInGui;
  if( SpecUtils::iequals_ascii( str, "FromDetectorPeakResponse" ) )
    return LoadingSource::FromDetectorPeakResponse;
  if( SpecUtils::iequals_ascii( str, "DefaultForDetectorType" ) )
    return LoadingSource::DefaultForDetectorType;
  if( SpecUtils::iequals_ascii( str, "FromSpectralData" ) )
    return LoadingSource::FromSpectralData;

  throw runtime_error( "PeakFitDetPrefs::loading_source_from_str: unrecognized value '" + str + "'" );
}//loading_source_from_str


void PeakFitDetPrefs::toXml( ::rapidxml::xml_node<char> *parent_node,
                              ::rapidxml::xml_document<char> *doc ) const
{
  assert( parent_node && doc );
  if( !parent_node || !doc )
    throw runtime_error( "PeakFitDetPrefs::toXml: null parent or doc" );

  rapidxml::xml_node<char> *base_node
    = doc->allocate_node( rapidxml::node_element, "PeakFitDetPrefs" );
  parent_node->append_node( base_node );

  const string versionstr = to_string( sm_xmlSerializationMajorVersion )
                            + "." + to_string( sm_xmlSerializationMinorVersion );
  rapidxml::xml_attribute<char> *ver_attr
    = doc->allocate_attribute( "version", doc->allocate_string( versionstr.c_str() ) );
  base_node->append_attribute( ver_attr );

  XmlUtils::append_string_node( base_node, "DetType", to_str( m_det_type ) );
  XmlUtils::append_string_node( base_node, "SkewType", PeakDef::to_string( m_peak_skew_type ) );

  static const char * const lower_names[] = { "LowerSkew0", "LowerSkew1", "LowerSkew2", "LowerSkew3" };
  static const char * const upper_names[] = { "UpperSkew0", "UpperSkew1", "UpperSkew2", "UpperSkew3" };

  for( int i = 0; i < 4; ++i )
  {
    if( m_lower_energy_skew[i].has_value() )
      XmlUtils::append_float_node( base_node, lower_names[i], m_lower_energy_skew[i].value() );
    if( m_upper_energy_skew[i].has_value() )
      XmlUtils::append_float_node( base_node, upper_names[i], m_upper_energy_skew[i].value() );
  }//for( int i = 0; i < 4; ++i )

  XmlUtils::append_bool_node( base_node, "RoiIndepSkew", m_roi_independent_skew );
  XmlUtils::append_string_node( base_node, "FwhmMethod", to_str( m_fwhm_method ) );
  XmlUtils::append_string_node( base_node, "Source", to_str( m_source ) );
}//toXml


void PeakFitDetPrefs::fromXml( const ::rapidxml::xml_node<char> *node )
{
  assert( node );
  if( !node )
    throw runtime_error( "PeakFitDetPrefs::fromXml: null node" );

  const rapidxml::xml_attribute<char> *ver_attr = node->first_attribute( "version" );
  if( !ver_attr || !ver_attr->value() )
    throw runtime_error( "PeakFitDetPrefs::fromXml: missing version attribute" );

  int majorVersion = -1;
  if( !(stringstream( ver_attr->value() ) >> majorVersion) )
    throw runtime_error( "PeakFitDetPrefs::fromXml: invalid version" );

  if( majorVersion != sm_xmlSerializationMajorVersion )
    throw runtime_error( "PeakFitDetPrefs::fromXml: unsupported major version "
                          + to_string( majorVersion ) );

  const string det_type_str = XmlUtils::get_string_node_value( node, "DetType" );
  m_det_type = coarse_res_from_str( det_type_str );

  const string skew_type_str = XmlUtils::get_string_node_value( node, "SkewType" );
  m_peak_skew_type = PeakDef::skew_from_string( skew_type_str );

  // Reset skew parameters
  for( int i = 0; i < 4; ++i )
  {
    m_lower_energy_skew[i].reset();
    m_upper_energy_skew[i].reset();
  }

  // Parse optional skew parameters
  static const char * const lower_names[] = { "LowerSkew0", "LowerSkew1", "LowerSkew2", "LowerSkew3" };
  static const char * const upper_names[] = { "UpperSkew0", "UpperSkew1", "UpperSkew2", "UpperSkew3" };

  for( int i = 0; i < 4; ++i )
  {
    // Use XML_FIRST_NODE directly rather than get_float_node_value (which throws if missing)
    const rapidxml::xml_node<char> *lower_node = node->first_node( lower_names[i] );
    if( lower_node && lower_node->value_size() )
    {
      double val = 0.0;
      if( SpecUtils::parse_double( lower_node->value(), lower_node->value_size(), val ) )
        m_lower_energy_skew[i] = val;
      else
        throw runtime_error( "PeakFitDetPrefs::fromXml: invalid value for " + string( lower_names[i] ) );
    }//if( lower_node )

    const rapidxml::xml_node<char> *upper_node = node->first_node( upper_names[i] );
    if( upper_node && upper_node->value_size() )
    {
      double val = 0.0;
      if( SpecUtils::parse_double( upper_node->value(), upper_node->value_size(), val ) )
        m_upper_energy_skew[i] = val;
      else
        throw runtime_error( "PeakFitDetPrefs::fromXml: invalid value for " + string( upper_names[i] ) );
    }//if( upper_node )
  }//for( int i = 0; i < 4; ++i )

  // ROI-independent skew (default true for backward compatibility)
  const rapidxml::xml_node<char> *roi_indep_node = node->first_node( "RoiIndepSkew" );
  if( roi_indep_node && roi_indep_node->value_size() )
    m_roi_independent_skew = (string( roi_indep_node->value(), roi_indep_node->value_size() ) == "1");
  else
    m_roi_independent_skew = true;

  // FWHM method (default Normal for backward compat with version 0.0)
  const rapidxml::xml_node<char> *fwhm_node = node->first_node( "FwhmMethod" );
  if( fwhm_node && fwhm_node->value_size() )
    m_fwhm_method = fwhm_method_from_str( string( fwhm_node->value(), fwhm_node->value_size() ) );
  else
    m_fwhm_method = FwhmMethod::Normal;

  const string source_str = XmlUtils::get_string_node_value( node, "Source" );
  m_source = loading_source_from_str( source_str );
}//fromXml


string PeakFitDetPrefs::toUrlQueryParts() const
{
  string result;

  // Always include det type and skew type (they're small and important)
  result += "DT=";
  result += to_str( m_det_type );

  result += "&SK=";
  result += PeakDef::to_label( m_peak_skew_type );

  // FWHM method - only include if non-default
  if( m_fwhm_method != FwhmMethod::Normal )
  {
    result += "&FM=";
    result += to_str( m_fwhm_method );
  }

  // Skew parameters - only include if set
  static const char * const lower_keys[] = { "LS0", "LS1", "LS2", "LS3" };
  static const char * const upper_keys[] = { "US0", "US1", "US2", "US3" };

  for( int i = 0; i < 4; ++i )
  {
    if( m_lower_energy_skew[i].has_value() )
    {
      char buf[64];
      snprintf( buf, sizeof( buf ), "&%s=%1.8e", lower_keys[i], m_lower_energy_skew[i].value() );
      result += buf;
    }

    if( m_upper_energy_skew[i].has_value() )
    {
      char buf[64];
      snprintf( buf, sizeof( buf ), "&%s=%1.8e", upper_keys[i], m_upper_energy_skew[i].value() );
      result += buf;
    }
  }//for( int i = 0; i < 4; ++i )

  return result;
}//toUrlQueryParts


void PeakFitDetPrefs::fromUrlQueryParts( const map<string,string> &parts )
{
  // Det type
  {
    const auto it = parts.find( "DT" );
    if( it == parts.end() )
      throw runtime_error( "PeakFitDetPrefs::fromUrlQueryParts: missing DT key" );
    m_det_type = coarse_res_from_str( it->second );
  }

  // Skew type
  {
    const auto it = parts.find( "SK" );
    if( it == parts.end() )
      throw runtime_error( "PeakFitDetPrefs::fromUrlQueryParts: missing SK key" );
    m_peak_skew_type = PeakDef::skew_from_string( it->second );
  }

  // Reset and parse optional skew parameters
  for( int i = 0; i < 4; ++i )
  {
    m_lower_energy_skew[i].reset();
    m_upper_energy_skew[i].reset();
  }

  static const char * const lower_keys[] = { "LS0", "LS1", "LS2", "LS3" };
  static const char * const upper_keys[] = { "US0", "US1", "US2", "US3" };

  for( int i = 0; i < 4; ++i )
  {
    {
      const auto it = parts.find( lower_keys[i] );
      if( it != parts.end() )
      {
        double val = 0.0;
        if( !SpecUtils::parse_double( it->second.c_str(), it->second.size(), val ) )
          throw runtime_error( "PeakFitDetPrefs::fromUrlQueryParts: invalid " + string( lower_keys[i] ) );
        m_lower_energy_skew[i] = val;
      }
    }

    {
      const auto it = parts.find( upper_keys[i] );
      if( it != parts.end() )
      {
        double val = 0.0;
        if( !SpecUtils::parse_double( it->second.c_str(), it->second.size(), val ) )
          throw runtime_error( "PeakFitDetPrefs::fromUrlQueryParts: invalid " + string( upper_keys[i] ) );
        m_upper_energy_skew[i] = val;
      }
    }
  }//for( int i = 0; i < 4; ++i )

  // FWHM method (optional, defaults to Normal)
  {
    const auto it = parts.find( "FM" );
    if( it != parts.end() )
      m_fwhm_method = fwhm_method_from_str( it->second );
    else
      m_fwhm_method = FwhmMethod::Normal;
  }

  // Source defaults to Default for URL-loaded prefs
  m_source = LoadingSource::Default;
}//fromUrlQueryParts


shared_ptr<const PeakFitDetPrefs>
PeakFitDetPrefs::defaultForDetectorType( const int specutils_det_type,
                                          const string &manufacturer,
                                          const string &model )
{
  using DT = SpecUtils::DetectorType;
  using CRT = PeakFitUtils::CoarseResolutionType;

  const DT det_type = static_cast<DT>( specutils_det_type );

  CRT coarse_type = CRT::Unknown;

  switch( det_type )
  {
    // HPGe detectors
    case DT::Fulcrum:
    case DT::Fulcrum40h:
    case DT::DetectiveUnknown:
    case DT::DetectiveEx:
    case DT::DetectiveEx100:
    case DT::DetectiveEx200:
    case DT::DetectiveX:
    case DT::MicroDetective:
    case DT::Falcon5000:
      coarse_type = CRT::High;
      break;

    // LaBr3 detectors
    case DT::IdentiFinderLaBr3:
    case DT::IdentiFinderR425LaBr:
    case DT::IdentiFinderR500LaBr:
    case DT::RadHunterLaBr3:
    case DT::OrtecRadEagleLaBr:
    case DT::Sam940LaBr3:
    case DT::RIIDEyeLaBr:
    case DT::RadSeekerLaBr:
    case DT::VerifinderLaBr:
      coarse_type = CRT::LaBr;
      break;

    // CZT detectors
    case DT::MicroRaider:
    case DT::KromekGR1:
      coarse_type = CRT::CZT;
      break;

    // CeBr3 - resolution between LaBr and NaI; use LaBr as closer match
    case DT::OrtecRadEagleCeBr2Inch:
    case DT::OrtecRadEagleCeBr3Inch:
    case DT::KromekD5:
      coarse_type = CRT::LaBr;
      break;

    // NaI detectors
    case DT::Exploranium:
    case DT::IdentiFinder:
    case DT::IdentiFinderNG:
    case DT::IdentiFinderTungsten:
    case DT::IdentiFinderR425NaI:
    case DT::IdentiFinderR500NaI:
    case DT::IdentiFinderUnknown:
    case DT::SAIC8:
    case DT::RadiaCodeCsI10:
    case DT::RadiaCodeCsI14:
    case DT::RadiaCodeGAGG10:
    case DT::Raysid:
    case DT::RadHunterNaI:
    case DT::Rsi701:
    case DT::Rsi705:
    case DT::AvidRsi:
    case DT::OrtecRadEagleNai:
    case DT::Sam940:
    case DT::Sam945:
    case DT::Sam950:
    case DT::Srpm210:
    case DT::RIIDEyeNaI:
    case DT::RadSeekerNaI:
    case DT::VerifinderNaI:
    case DT::KromekD3S:
      coarse_type = CRT::Low;
      break;

    case DT::Interceptor:
    case DT::Unknown:
      // Fall through to keyword matching
      break;
  }//switch( det_type )


  // If still unknown, try keyword matching on manufacturer/model
  if( coarse_type == CRT::Unknown )
  {
    const string mfr_model = manufacturer + " " + model;

    if( SpecUtils::icontains( mfr_model, "HPGe" )
       || SpecUtils::icontains( mfr_model, "Germanium" )
       || SpecUtils::icontains( mfr_model, "Detective" )
       || SpecUtils::icontains( mfr_model, "Falcon" )
       || SpecUtils::icontains( mfr_model, "Fulcrum" ) )
    {
      coarse_type = CRT::High;
    }
    else if( SpecUtils::icontains( mfr_model, "LaBr" )
            || SpecUtils::icontains( mfr_model, "Lanthanum" ) )
    {
      coarse_type = CRT::LaBr;
    }
    else if( SpecUtils::icontains( mfr_model, "CZT" )
            || SpecUtils::icontains( mfr_model, "Cadmium Zinc" )
            || SpecUtils::icontains( mfr_model, "CdZnTe" ) )
    {
      coarse_type = CRT::CZT;
    }
    else if( SpecUtils::icontains( mfr_model, "NaI" )
            || SpecUtils::icontains( mfr_model, "Sodium Iodide" )
            || SpecUtils::icontains( mfr_model, "CsI" )
            || SpecUtils::icontains( mfr_model, "Cesium Iodide" )
            || SpecUtils::icontains( mfr_model, "CeBr" )
            || SpecUtils::icontains( mfr_model, "GAGG" ) )
    {
      coarse_type = CRT::Low;
    }
  }//if( coarse_type == CRT::Unknown )

  if( coarse_type == CRT::Unknown )
    return nullptr;

  auto prefs = make_shared<PeakFitDetPrefs>();
  prefs->m_det_type = coarse_type;
  prefs->m_peak_skew_type = PeakDef::NoSkew;
  prefs->m_source = LoadingSource::DefaultForDetectorType;

  return prefs;
}//defaultForDetectorType


shared_ptr<const PeakFitDetPrefs>
PeakFitDetPrefs::guessFromSpectralData( const shared_ptr<const SpecUtils::Measurement> & )
{
  // Stub for future implementation
  return nullptr;
}//guessFromSpectralData


void apply_fwhm_method_to_peaks(
  vector<shared_ptr<PeakDef>> &peaks,
  const shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitDetPrefs &prefs,
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> &fit_options )
{
  if( prefs.m_fwhm_method == PeakFitDetPrefs::FwhmMethod::Normal )
    return;

  if( !drf || !drf->hasResolutionInfo() )
    return; // Silently fall back to Normal

  for( shared_ptr<PeakDef> &peak : peaks )
  {
    if( !peak )
      continue;

    const float energy = static_cast<float>( peak->mean() );
    if( prefs.m_fwhm_method == PeakFitDetPrefs::FwhmMethod::Normal )
    {
      const double drf_sigma = drf->peakResolutionSigma( energy );
      peak->setSigma( drf_sigma );

      if( prefs.m_fwhm_method == PeakFitDetPrefs::FwhmMethod::DetFwhm )
        peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
      // DetPlusRefine: leave fitFor(Sigma) as-is (true); SmallFwhmRefinementOnly
      // will be added to fit_options below.
    }
  }//for( peak : peaks )

  if( prefs.m_fwhm_method == PeakFitDetPrefs::FwhmMethod::DetPlusRefine )
    fit_options |= PeakFitLM::SmallFwhmRefinementOnly;
}//apply_fwhm_method_to_peaks
