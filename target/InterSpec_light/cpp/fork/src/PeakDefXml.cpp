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

/* InterSpec's peak XML format (the `<PeakContinuum>` and `<Peak>` elements it writes into the
 `<DHS:InterSpec>` node of N42-2012 files), ported from InterSpec's PeakDef.cpp.

 The only intentional difference is the peak source: InterSpec stores SandiaDecay/ReactionGamma
 pointers, and resolves them from these elements; here the same fields go to and from
 `PeakDef::Source`.  Only compiled when the LIGHT_PEAK_FILE_IO build option is on.
 */
#include "InterSpec_config.h"

#include <map>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakDists.h"

using namespace std;

namespace
{
  /** The version of InterSpec's peak XML this writes (InterSpec's `PeakDef::sm_xmlSerialization*`).
   For backward compatibility, peaks not using the Voigt/Bortel-family skews write version 0.2.
   */
  const int sm_peak_xml_major = 1;
  const int sm_peak_xml_minor = 1;
  const int sm_pre_voigt_peak_xml_major = 0;
  const int sm_pre_voigt_peak_xml_minor = 2;

  const float sm_annihilation_energy = 510.9989f;

  /** Replacement for rapidxml::internal::compare (which this rapidxml does not have). */
  bool compare( const char *a, const size_t alen, const char *b, const size_t blen, const bool case_sensitive )
  {
    if( alen != blen )
      return false;
    for( size_t i = 0; i < alen; ++i )
    {
      const bool same = case_sensitive ? (a[i] == b[i])
                          : (std::tolower( static_cast<unsigned char>(a[i]) )
                             == std::tolower( static_cast<unsigned char>(b[i]) ));
      if( !same )
        return false;
    }
    return true;
  }//compare(...)


  //clones 'source' into the document that 'result' is a part of.
  //  'result' is cleared and set lexically equal to 'source'.
  void clone_node_deep( const ::rapidxml::xml_node<char> *source, ::rapidxml::xml_node<char> *result )
  {
    using namespace ::rapidxml;

    xml_document<char> *doc = result->document();
    if( !doc )
      throw runtime_error( "clone_node_deep: insert result into document before calling" );

    result->remove_all_attributes();
    result->remove_all_nodes();
    result->type( source->type() );

    char *str = doc->allocate_string( source->name(), source->name_size() );
    result->name( str, source->name_size() );

    if( source->value() )
    {
      str = doc->allocate_string( source->value(), source->value_size() );
      result->value( str, source->value_size() );
    }

    for( xml_node<char> *child = source->first_node(); child; child = child->next_sibling() )
    {
      xml_node<char> *clone = doc->allocate_node( child->type() );
      result->append_node( clone );
      clone_node_deep( child, clone );
    }

    for( xml_attribute<char> *attr = source->first_attribute(); attr; attr = attr->next_attribute() )
    {
      const char *name = doc->allocate_string( attr->name(), attr->name_size() );
      const char *value = doc->allocate_string( attr->value(), attr->value_size() );
      xml_attribute<char> *clone = doc->allocate_attribute( name, value, attr->name_size(), attr->value_size() );
      result->append_attribute( clone );
    }
  }//void clone_node_deep(...)


  /** Reads a `"<value> <uncertainty>"` element, as written by `PeakDef::toXml`. */
  bool read_peak_value_node( const rapidxml::xml_node<char> * const peak_node,
                             const char * const name, const size_t name_len, double &value )
  {
    const rapidxml::xml_node<char> * const node = peak_node->first_node( name, name_len );
    double uncert = 0.0;
    return (node && node->value() && (sscanf( node->value(), "%lf %lf", &value, &uncert ) >= 1));
  }//read_peak_value_node(...)


  /** Sums over the peaks referencing continuum `cont_id`, for converting a pre-version-3
   BiLinearStepCDF (see `PeakContinuum::fromXml`).  `<Peak>` nodes are directly under `<Peaks>` for
   peak-serialization version 2, and under `<Peaks><PeakSet>` for version 1, so both are walked.
   */
  struct LegacyRoiPeakSums
  {
    double total_amp = 0.0;   //SUM_j( amp_j )
    double amp_cdf0 = 0.0;    //SUM_j( amp_j * CDF_j(anchor) )
  };


  LegacyRoiPeakSums legacy_roi_peak_sums( const rapidxml::xml_node<char> * const peaks_node,
                                          const int cont_id, const double anchor_energy )
  {
    using namespace rapidxml;

    LegacyRoiPeakSums sums;
    if( !peaks_node )
      return sums;

    const auto sum_peaks_under = [&sums,cont_id,anchor_energy]( const xml_node<char> * const parent ){
      for( const xml_node<char> *peak_node = parent->first_node("Peak",4);
           peak_node; peak_node = peak_node->next_sibling("Peak",4) )
      {
        const xml_attribute<char> *att = peak_node->first_attribute( "continuumID", 11 );
        int this_id = -1;
        if( !att || !att->value() || (sscanf(att->value(), "%i", &this_id) != 1) || (this_id != cont_id) )
          continue;

        double amp = 0.0, mean = 0.0, sigma = 0.0;
        if( !read_peak_value_node( peak_node, "Amplitude", 9, amp ) || !std::isfinite(amp) )
          continue;

        amp = (std::max)( amp, 0.0 );
        sums.total_amp += amp;

        if( !read_peak_value_node( peak_node, "Centroid", 8, mean )
           || !read_peak_value_node( peak_node, "Width", 5, sigma )
           || !std::isfinite(mean) || !std::isfinite(sigma) || (sigma <= 0.0) )
          continue;

        PeakDef::SkewType skew = PeakDef::SkewType::NoSkew;
        const xml_node<char> * const skew_node = peak_node->first_node( "Skew", 4 );
        if( skew_node && skew_node->value_size() )
        {
          try
          {
            skew = PeakDef::skew_from_string( std::string(skew_node->value(), skew_node->value_size()) );
          }catch( std::exception & )
          {
            skew = PeakDef::SkewType::NoSkew;   // e.g. the retired "LandauSkew"
          }
        }//if( skew_node && skew_node->value_size() )

        double skew_pars[PeakDef::CoefficientType::SkewPar5 - PeakDef::CoefficientType::SkewPar0 + 1] = { 0.0 };
        const size_t num_skew = PeakDef::num_skew_parameters( skew );
        bool have_skew_pars = true;
        for( size_t k = 0; (k < num_skew) && have_skew_pars; ++k )
        {
          const std::string name = "Skew" + std::to_string(k);
          have_skew_pars = read_peak_value_node( peak_node, name.c_str(), name.size(), skew_pars[k] );
        }

        if( !have_skew_pars )
          skew = PeakDef::SkewType::NoSkew;

        const double cdf0 = PeakDists::peak_cdf( anchor_energy, mean, sigma, skew, skew_pars );
        if( std::isfinite(cdf0) )
          sums.amp_cdf0 += amp * cdf0;
      }//for( loop over Peak nodes )
    };//sum_peaks_under

    sum_peaks_under( peaks_node );
    for( const xml_node<char> *set_node = peaks_node->first_node("PeakSet",7);
         set_node; set_node = set_node->next_sibling("PeakSet",7) )
      sum_peaks_under( set_node );

    return sums;
  }//LegacyRoiPeakSums legacy_roi_peak_sums(...)


  bool gamma_type_from_str( const char *val, const size_t len, PeakDef::SourceGammaType &type )
  {
    const pair<const char *,PeakDef::SourceGammaType> types[] = {
      { "NormalGamma", PeakDef::NormalGamma }, { "AnnihilationGamma", PeakDef::AnnihilationGamma },
      { "SingleEscapeGamma", PeakDef::SingleEscapeGamma }, { "DoubleEscapeGamma", PeakDef::DoubleEscapeGamma },
      { "XrayGamma", PeakDef::XrayGamma }
    };

    for( const auto &t : types )
    {
      if( compare( val, len, t.first, strlen(t.first), false ) )
      {
        type = t.second;
        return true;
      }
    }
    return false;
  }//gamma_type_from_str(...)
}//namespace


void PeakContinuum::toXml( rapidxml::xml_node<char> *parent, const int contId ) const
{
  using namespace rapidxml;

  xml_document<char> *doc = parent ? parent->document() : (xml_document<char> *)0;
  if( !doc )
    throw runtime_error( "PeakContinuum::toXml(...): invalid input" );

  char buffer[128];
  xml_node<char> *node = 0;
  xml_node<char> *cont_node = doc->allocate_node( node_element, "PeakContinuum" );

  static_assert( PeakContinuum::sm_xmlSerializationVersion == 3,
                "PeakContinuum::toXml needs to be updated for new serialization version." );

  // Write the lowest version that can hold this continuum type, so older InterSpec can read it
  int version = PeakContinuum::sm_xmlSerializationVersion;
  switch( m_type )
  {
    case NoOffset: case External: case Constant: case Linear: case Quadratic: case Cubic:
      version = 0;
      break;
    case FlatStep: case LinearStep: case BiLinearStep:
      version = 1;
      break;
    case FlatStepCDF: case LinearStepCDF:
      version = 2;
      break;
    case BiLinearStepCDF:
      version = 3;
      break;
  }//switch( m_type )

  snprintf( buffer, sizeof(buffer), "%i", version );
  const char *val = doc->allocate_string( buffer );
  xml_attribute<char> *att = doc->allocate_attribute( "version", val );
  cont_node->append_attribute( att );

  snprintf( buffer, sizeof(buffer), "%i", contId );
  val = doc->allocate_string( buffer );
  att = doc->allocate_attribute( "id", val );
  cont_node->append_attribute( att );

  parent->append_node( cont_node );

  const char *type = offset_type_str( m_type );
  node = doc->allocate_node( node_element, "Type", type );
  cont_node->append_node( node );

  snprintf( buffer, sizeof(buffer), "%1.8e", m_lowerEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "LowerEnergy", val );
  cont_node->append_node( node );

  snprintf( buffer, sizeof(buffer), "%1.8e", m_upperEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "UpperEnergy", val );
  cont_node->append_node( node );

  snprintf( buffer, sizeof(buffer), "%1.8e", m_referenceEnergy );
  val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "ReferenceEnergy", val );
  cont_node->append_node( node );

  if( (m_type != NoOffset) && (m_type != External) )
  {
    stringstream valsstrm, uncertstrm, fitstrm;
    for( size_t i = 0; i < m_values.size(); ++i )
    {
      const char *spacer = (i ? " " : "");
      snprintf( buffer, sizeof(buffer), "%1.8e", m_values[i] );
      valsstrm << spacer << buffer;

      snprintf( buffer, sizeof(buffer), "%1.8e", m_uncertainties[i] );
      uncertstrm << spacer << buffer;

      fitstrm << spacer << (m_fitForValue[i] ? '1': '0');
    }//for( size_t i = 0; i < m_values.size(); ++i )

    xml_node<char> *coeffs_node = doc->allocate_node( node_element, "Coefficients" );
    cont_node->append_node( coeffs_node );

    val = doc->allocate_string( valsstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Values", val );
    coeffs_node->append_node( node );

    val = doc->allocate_string( uncertstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Uncertainties", val );
    coeffs_node->append_node( node );

    val = doc->allocate_string( fitstrm.str().c_str() );
    node = doc->allocate_node( node_element, "Fittable", val );
    coeffs_node->append_node( node );
  }//if( m_type != NoOffset && m_type != External )

  if( m_externalContinuum )
  {
    stringstream contXml;
    m_externalContinuum->write_2006_N42_xml( contXml );

    // Parse the XML, and insert it into this document
    const string datastr = contXml.str();
    std::unique_ptr<char []> data( new char [datastr.size()+1] );
    memcpy( data.get(), datastr.c_str(), datastr.size()+1 );

    xml_document<char> contdoc;
    contdoc.parse<rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace>( data.get() );

    node = doc->allocate_node( node_element, "ExternalContinuum" );
    cont_node->append_node( node );

    xml_node<char> *spec_node = contdoc.first_node( "Measurement", 11 );
    if( !spec_node )
      throw runtime_error( "Didnt get expected Measurement node" );
    spec_node = spec_node->first_node( "Spectrum", 8 );
    if( !spec_node )
      throw runtime_error( "Didnt get expected Spectrum node" );

    xml_node<char> *new_spec_node = doc->allocate_node( node_element );
    node->append_node( new_spec_node );

    clone_node_deep( spec_node, new_spec_node );
  }//if( m_externalContinuum )
}//void PeakContinuum::toXml(...)


void PeakContinuum::fromXml( const rapidxml::xml_node<char> *cont_node, int &contId )
{
  using namespace rapidxml;

  if( !cont_node )
    throw runtime_error( "PeakContinuum::fromXml(...): invalid input" );

  if( !compare( cont_node->name(), cont_node->name_size(), "PeakContinuum", 13, false ) )
    throw std::logic_error( "PeakContinuum::fromXml(...): invalid input node name" );

  xml_attribute<char> *att = cont_node->first_attribute( "version", 7 );

  int version;
  if( !att || !att->value() || (sscanf(att->value(), "%i", &version) != 1) )
    throw runtime_error( "PeakContinuum invalid version" );

  static_assert( PeakContinuum::sm_xmlSerializationVersion == 3,
                "PeakContinuum::fromXml needs to be updated for new serialization version." );

  // Versions 1 and 2 only add continuum types; version 3 only changes what BiLinearStepCDF's
  //  parameters mean, which is converted below.
  if( (version < 0) || (version > PeakContinuum::sm_xmlSerializationVersion) )
    throw runtime_error( "Invalid PeakContinuum version: " + std::to_string(version) + ".  "
                    + "Only up to version " + to_string(PeakContinuum::sm_xmlSerializationVersion)
                    + " supported." );

  att = cont_node->first_attribute( "id", 2 );
  if( !att || !att->value() || (sscanf(att->value(), "%i", &contId) != 1) )
    throw runtime_error( "PeakContinuum invalid ID" );

  xml_node<char> *node = cont_node->first_node( "Type", 4 );
  if( !node || !node->value() )
    throw runtime_error( "PeakContinuum not Type node" );

  m_type = str_to_offset_type_str( node->value(), node->value_size() );

  float dummyval;
  node = cont_node->first_node( "LowerEnergy", 11 );
  if( node )
  {
    if( !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
      throw runtime_error( "Continuum didnt have valid LowerEnergy" );

    m_lowerEnergy = dummyval;

    node = cont_node->first_node( "UpperEnergy", 11 );
    if( !node || !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
      throw runtime_error( "Continuum didnt have UpperEnergy" );
    m_upperEnergy = dummyval;
  }else
  {
    m_lowerEnergy = m_upperEnergy = 0.0;
    node = cont_node->first_node( "UpperEnergy", 11 );
    if( node )
      throw runtime_error( "Continuum didnt have LowerEnergy, but did have UpperEnergy" );
  }//if( have <LowerEnergy> ) / else

  node = cont_node->first_node( "ReferenceEnergy", 15 );
  if( node )
  {
    if( !node->value() || (sscanf(node->value(),"%e",&dummyval) != 1) )
      throw runtime_error( "Continuum didnt have ReferenceEnergy" );
    m_referenceEnergy = dummyval;
  }else
  {
    if( m_lowerEnergy != m_upperEnergy )
      throw runtime_error( "Continuum didnt have ReferenceEnergy, but did have energy range defined" );
    m_referenceEnergy = 0.0;
  }//if( have <ReferenceEnergy> ) / else

  if( (m_type != NoOffset) && (m_type != External) )
  {
    xml_node<char> *coeffs_node = cont_node->first_node( "Coefficients", 12 );
    if( !coeffs_node )
      throw runtime_error( "Continuum didnt have Coefficients node" );

    std::vector<float> contents;
    node = coeffs_node->first_node( "Values", 6 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Values" );

    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_values.assign( begin(contents), end(contents) );

    node = coeffs_node->first_node( "Uncertainties", 13 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Uncertainties" );

    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_uncertainties.assign( begin(contents), end(contents) );

    node = coeffs_node->first_node( "Fittable", 8 );
    if( !node || !node->value() )
      throw runtime_error( "Continuum didnt have Coefficient Fittable" );

    SpecUtils::split_to_floats( node->value(), node->value_size(), contents );
    m_fitForValue.resize( contents.size() );
    for( size_t i = 0; i < contents.size(); ++i )
      m_fitForValue[i] = (contents[i] > 0.5f);

    if( (m_values.size() != m_uncertainties.size()) || (m_fitForValue.size() != m_values.size()) )
      throw runtime_error( "Continuum coefficients not consistent" );

    // Code elsewhere indexes the coefficients by type alone
    const size_t num_expected = PeakContinuum::num_parameters( m_type );
    if( m_values.size() != num_expected )
      throw runtime_error( "PeakContinuum::fromXml: continuum type '"
                          + string(offset_type_str(m_type)) + "' expects "
                          + std::to_string(num_expected) + " coefficients, but XML had "
                          + std::to_string(m_values.size()) + "." );
  }else
  {
    m_values.clear();
    m_uncertainties.clear();
    m_fitForValue.clear();
  }//if( m_type != NoOffset && m_type != External ) / else

  node = cont_node->first_node( "ExternalContinuum", 17 );
  if( node )
  {
    node = node->first_node( "Spectrum", 8 );
    if( !node )
      throw runtime_error( "Spectrum node expected under ExternalContinuum" );
    std::shared_ptr<SpecUtils::Measurement> meas = std::make_shared<SpecUtils::Measurement>();
    meas->set_info_from_2006_N42_spectrum_node( node );
    m_externalContinuum = meas;
  }//if( node )

  // Versions 2 and earlier stored BiLinearStepCDF differently; the ROI's peak sums needed to convert
  //  are in the <Peak> nodes alongside this one.
  if( (version < 3) && (m_type == BiLinearStepCDF) )
  {
    assert( m_values.size() == 4 );
    const LegacyRoiPeakSums sums = legacy_roi_peak_sums( cont_node->parent(), contId, m_lowerEnergy );
    convert_legacy_bilinear_step_cdf( m_values, m_uncertainties, sums.total_amp, sums.amp_cdf0 );
  }//if( a pre-version-3 BiLinearStepCDF )
}//void PeakContinuum::fromXml(...)


void PeakContinuum::convert_legacy_bilinear_step_cdf( vector<double> &values,
                                                      vector<double> &uncertainties,
                                                      const double total_amp,
                                                      const double amp_cdf0 )
{
  if( values.size() != 4 )
    return;

  double step0 = 0.0, step1 = 0.0;
  if( std::isfinite(total_amp) && (total_amp > 0.0) )
  {
    step0 = (values[2] - values[0]) / total_amp;
    step1 = (values[3] - values[1]) / total_amp;
  }

  if( !std::isfinite(step0) )
    step0 = 0.0;
  if( !std::isfinite(step1) )
    step1 = 0.0;

  // Version 2 blended with the CDF measured from -infinity, so it carried a constant offset of
  //  `step_k * SUM_j(amp_j*CDF_j(roi_lower))` that the ROI-anchored version-3 model does not;
  //  folding it into the polynomial keeps the continuum the same shape.
  if( std::isfinite(amp_cdf0) )
  {
    values[0] += step0 * amp_cdf0;
    values[1] += step1 * amp_cdf0;
  }

  values[2] = step0;
  values[3] = step1;

  // The old uncertainties were on the right-hand line, not on a step
  if( uncertainties.size() == 4 )
    uncertainties[2] = uncertainties[3] = 0.0;
}//void PeakContinuum::convert_legacy_bilinear_step_cdf(...)


rapidxml::xml_node<char> *PeakDef::toXml( rapidxml::xml_node<char> *parent,
                                          rapidxml::xml_node<char> *continuum_parent,
                                          std::map<std::shared_ptr<PeakContinuum>,int> &continuums ) const
{
  using namespace rapidxml;

  xml_document<char> *doc = parent ? parent->document() : (xml_document<char> *)0;
  if( !doc )
    throw runtime_error( "PeakDef::toXml(...): invalid input" );

  if( !m_continuum )
    throw logic_error( "PeakDef::toXml(...): continuum should be valid" );

  if( !continuums.count(m_continuum) )
  {
    const int index = static_cast<int>( continuums.size() + 1 );
    m_continuum->toXml( continuum_parent, index );
    continuums[m_continuum] = index;
  }//if( !continuums.count(m_continuum) )

  char buffer[128];
  const int contID = continuums[m_continuum];

  xml_node<char> *node = 0;
  xml_node<char> *peak_node = doc->allocate_node( node_element, "Peak" );

  // For backward compatibility, write version 0.2 for peaks without newer skew types
  const bool uses_new_skew = (m_skewType == SkewType::VoigtPlusBortel)
                             || (m_skewType == SkewType::GaussPlusBortel)
                             || (m_skewType == SkewType::DoubleBortel)
                             || (m_skewType == SkewType::GadrasGeneric)
                             || (m_skewType == SkewType::GadrasCZT);
  const int majorVersion = uses_new_skew ? sm_peak_xml_major : sm_pre_voigt_peak_xml_major;
  const int minorVersion = uses_new_skew ? sm_peak_xml_minor : sm_pre_voigt_peak_xml_minor;

  snprintf( buffer, sizeof(buffer), "%i.%i", majorVersion, minorVersion );
  const char *val = doc->allocate_string( buffer );
  xml_attribute<char> *att = doc->allocate_attribute( "version", val );
  peak_node->append_attribute( att );

  snprintf( buffer, sizeof(buffer), "%i", contID );
  val = doc->allocate_string( buffer );
  att = doc->allocate_attribute( "continuumID", val );
  peak_node->append_attribute( att );

  parent->append_node( peak_node );

  if( m_userLabel.size() )
  {
    val = doc->allocate_string( m_userLabel.c_str() );
    node = doc->allocate_node( node_element, "UserLabel", val );
    peak_node->append_node( node );
  }//if( m_userLabel.size() )

  val = PeakDef::to_str( m_type );
  node = doc->allocate_node( node_element, "Type", val );
  peak_node->append_node( node );

  switch( m_skewType )
  {
    case NumSkewType:
      assert( 0 );
      [[fallthrough]];
    case NoSkew:                 val = "NoSkew";                 break;
    case Bortel:                 val = "ExGauss";                break;
    case DoubleBortel:           val = "DoubleBortel";           break;
    case GaussPlusBortel:        val = "GaussPlusBortel";        break;
    case GaussExp:               val = "GaussExp";               break;
    case CrystalBall:            val = "CrystalBall";            break;
    case ExpGaussExp:            val = "ExpGaussExp";            break;
    case DoubleSidedCrystalBall: val = "DoubleSidedCrystalBall"; break;
    case VoigtPlusBortel:        val = "VoigtPlusBortel";        break;
    case GadrasGeneric:          val = "GadrasGeneric";          break;
    case GadrasCZT:              val = "GadrasCZT";              break;
  }//switch( m_skewType )

  node = doc->allocate_node( node_element, "Skew", val );
  peak_node->append_node( node );

  for( CoefficientType t = CoefficientType(0); t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const char *label = to_string( t );

    // InterSpec v1.0.11 and before always expect the "LandauAmplitude", "LandauMode", and
    //  "LandauSigma" elements, so they are written for peaks without skew.
    const char *extra_label = nullptr;
    switch( t )
    {
      case Mean: case Sigma: case GaussAmplitude:
      case Chi2DOF: case NumCoefficientTypes:
        break;

      case SkewPar0:
        if( m_skewType == NoSkew )
        {
          extra_label = "LandauAmplitude";
          label = nullptr;
        }
        break;

      case SkewPar1:
        if( m_skewType == NoSkew )
          extra_label = "LandauMode";
        if( (m_skewType != CrystalBall)
           && (m_skewType != DoubleSidedCrystalBall)
           && (m_skewType != ExpGaussExp)
           && (m_skewType != VoigtPlusBortel)
           && (m_skewType != GaussPlusBortel)
           && (m_skewType != DoubleBortel)
           && (m_skewType != GadrasGeneric)
           && (m_skewType != GadrasCZT) )
        {
          label = nullptr;
        }
        break;

      case SkewPar2:
        if( m_skewType == NoSkew )
          extra_label = "LandauSigma";
        if( (m_skewType != DoubleSidedCrystalBall)
           && (m_skewType != VoigtPlusBortel)
           && (m_skewType != DoubleBortel)
           && (m_skewType != GadrasGeneric)
           && (m_skewType != GadrasCZT) )
          label = nullptr;
        break;

      case SkewPar3:
        if( (m_skewType != DoubleSidedCrystalBall)
           && (m_skewType != GadrasGeneric)
           && (m_skewType != GadrasCZT) )
          label = nullptr;
        break;

      case SkewPar4:
      case SkewPar5:
        if( (m_skewType != GadrasGeneric) && (m_skewType != GadrasCZT) )
          label = nullptr;
        break;
    }//switch( t )

    auto add_node = [this, &buffer, doc, peak_node, t]( const char *label, const bool is_deprecated ){
      snprintf( buffer, sizeof(buffer), "%1.8e %1.8e", m_coefficients[t], m_uncertainties[t] );
      const char *val = doc->allocate_string( buffer );
      xml_node<char> *node = doc->allocate_node( node_element, label, val );
      xml_attribute<char> *att = doc->allocate_attribute( "fit", (m_fitFor[t] ? "true" : "false") );
      node->append_attribute( att );

      if( is_deprecated )
      {
        att = doc->allocate_attribute( "remark", "Deprecated element" );
        node->append_attribute( att );
      }

      peak_node->append_node( node );
    };//add_node

    if( label )
      add_node( label, false );
    if( extra_label )
      add_node( extra_label, true );
  }//for( loop over coefficients )

  // InterSpec reads "forCalibration"; "useForEnergyCalibration" is written for the future.
  att = doc->allocate_attribute( "forCalibration", (m_useForEnergyCal ? "true" : "false") );
  peak_node->append_attribute( att );

  att = doc->allocate_attribute( "useForEnergyCalibration", (m_useForEnergyCal ? "true" : "false") );
  peak_node->append_attribute( att );

  att = doc->allocate_attribute( "source", (m_useForShieldingSourceFit ? "true" : "false") );
  peak_node->append_attribute( att );

  att = doc->allocate_attribute( "useForManualRelEff", (m_useForManualRelEff ? "true" : "false") );
  peak_node->append_attribute( att );

  // The DRF-fit flags are only written when not their default values
  if( m_useForDrfIntrinsicEffFit != PeakDef::sm_defaultUseForDrfIntrinsicEffFit )
  {
    att = doc->allocate_attribute( "useForDrfIntrinsicEffFit", (m_useForDrfIntrinsicEffFit ? "true" : "false") );
    peak_node->append_attribute( att );
  }

  if( m_useForDrfFwhmFit != PeakDef::sm_defaultUseForDrfFwhmFit )
  {
    att = doc->allocate_attribute( "useForDrfFwhmFit", (m_useForDrfFwhmFit ? "true" : "false") );
    peak_node->append_attribute( att );
  }

  if( m_useForDrfDepthOfInteractionFit != PeakDef::sm_defaultUseForDrfDepthOfInteractionFit )
  {
    att = doc->allocate_attribute( "useForDrfDepthOfInteractionFit", (m_useForDrfDepthOfInteractionFit ? "true" : "false") );
    peak_node->append_attribute( att );
  }

  if( !m_lineColor.isDefault() )
  {
    val = doc->allocate_string( m_lineColor.cssText(false).c_str() );
    node = doc->allocate_node( node_element, "LineColor", val );
    peak_node->append_node( node );
  }//if( !m_lineColor.isDefault() )

  const char *gammaTypeVal = "NormalGamma";
  switch( m_source.gamma_type )
  {
    case PeakDef::NormalGamma:       gammaTypeVal = "NormalGamma";       break;
    case PeakDef::AnnihilationGamma: gammaTypeVal = "AnnihilationGamma"; break;
    case PeakDef::SingleEscapeGamma: gammaTypeVal = "SingleEscapeGamma"; break;
    case PeakDef::DoubleEscapeGamma: gammaTypeVal = "DoubleEscapeGamma"; break;
    case PeakDef::XrayGamma:         gammaTypeVal = "XrayGamma";         break;
  }//switch( m_source.gamma_type )

  switch( m_source.kind )
  {
    case Source::Kind::None:
      break;

    case Source::Kind::Nuclide:
    {
      xml_node<char> *nuc_node = doc->allocate_node( node_element, "Nuclide" );
      peak_node->append_node( nuc_node );

      val = doc->allocate_string( m_source.name.c_str() );
      node = doc->allocate_node( node_element, "Name", val );
      nuc_node->append_node( node );

      // InterSpec finds the transition from its parent, child, and photon energy.  Annihilation
      //  photons are left without a transition (the energy InterSpec would need is the positron's),
      //  which InterSpec accepts.
      const string &trans_parent = m_source.decay_parent.empty() ? m_source.name : m_source.decay_parent;
      if( (m_source.gamma_type != PeakDef::AnnihilationGamma) && !m_source.decay_child.empty() )
      {
        val = doc->allocate_string( trans_parent.c_str() );
        node = doc->allocate_node( node_element, "DecayParent", val );
        nuc_node->append_node( node );

        val = doc->allocate_string( m_source.decay_child.c_str() );
        node = doc->allocate_node( node_element, "DecayChild", val );
        nuc_node->append_node( node );

        snprintf( buffer, sizeof(buffer), "%1.8e", static_cast<double>(m_source.particle_energy) );
        val = doc->allocate_string( buffer );
        node = doc->allocate_node( node_element, "DecayGammaEnergy", val );
        nuc_node->append_node( node );
      }//if( we can specify the transition )

      node = doc->allocate_node( node_element, "DecayGammaType", gammaTypeVal );
      nuc_node->append_node( node );
      break;
    }//case Source::Kind::Nuclide:

    case Source::Kind::Xray:
    {
      xml_node<char> *xray_node = doc->allocate_node( node_element, "XRay" );
      peak_node->append_node( xray_node );

      val = doc->allocate_string( m_source.name.c_str() );
      node = doc->allocate_node( node_element, "Element", val );
      xray_node->append_node( node );

      snprintf( buffer, sizeof(buffer), "%1.8e", static_cast<double>(m_source.particle_energy) );
      val = doc->allocate_string( buffer );
      node = doc->allocate_node( node_element, "Energy", val );
      xray_node->append_node( node );
      break;
    }//case Source::Kind::Xray:

    case Source::Kind::Reaction:
    {
      xml_node<char> *rctn_node = doc->allocate_node( node_element, "Reaction" );
      peak_node->append_node( rctn_node );

      val = doc->allocate_string( m_source.name.c_str() );
      node = doc->allocate_node( node_element, "Name", val );
      rctn_node->append_node( node );

      snprintf( buffer, sizeof(buffer), "%1.8e", static_cast<double>(m_source.particle_energy) );
      val = doc->allocate_string( buffer );
      node = doc->allocate_node( node_element, "Energy", val );
      rctn_node->append_node( node );

      node = doc->allocate_node( node_element, "Type", gammaTypeVal );
      rctn_node->append_node( node );
      break;
    }//case Source::Kind::Reaction:
  }//switch( m_source.kind )

  return peak_node;
}//rapidxml::xml_node<char> *toXml(...)


void PeakDef::fromXml( const rapidxml::xml_node<char> *peak_node,
                       const std::map<int,std::shared_ptr<PeakContinuum> > &continuums,
                       std::vector<std::string> *warnings )
{
  using namespace rapidxml;

  if( !peak_node )
    throw logic_error( "PeakDef::fromXml(...): invalid input node" );

  if( !compare( peak_node->name(), peak_node->name_size(), "Peak", 4, false ) )
    throw std::logic_error( "PeakDef::fromXml(...): invalid input node name" );

  reset();

  int contID;
  xml_attribute<char> *att = peak_node->first_attribute( "continuumID", 11 );
  if( !att )
    throw runtime_error( "No continuum ID" );
  if( sscanf( att->value(), "%i", &contID ) != 1 )
    throw runtime_error( "Non integer continuum ID" );

  const auto contpos = continuums.find( contID );
  if( contpos == continuums.end() )
    throw runtime_error( "Couldnt find valid continuum for peak" );

  m_continuum = contpos->second;

  const auto read_bool_att = [peak_node]( const char *name, bool &value ) -> bool {
    const xml_attribute<char> *att = peak_node->first_attribute( name, strlen(name) );
    if( !att )
      return false;
    value = compare( att->value(), att->value_size(), "true", 4, false );
    if( !value && !compare( att->value(), att->value_size(), "false", 5, false ) )
      throw runtime_error( string("invalid ") + name + " value" );
    return true;
  };//read_bool_att

  if( !read_bool_att( "forCalibration", m_useForEnergyCal ) )
    throw runtime_error( "missing forCalibration attribute" );

  if( !read_bool_att( "source", m_useForShieldingSourceFit ) )
    throw runtime_error( "missing source attribute" );

  // useForManualRelEff was added June 2022, so it is not required
  if( !read_bool_att( "useForManualRelEff", m_useForManualRelEff ) )
    m_useForManualRelEff = true;

  m_useForDrfIntrinsicEffFit = PeakDef::sm_defaultUseForDrfIntrinsicEffFit;
  att = peak_node->first_attribute( "useForDrfIntrinsicEffFit", 24 );
  if( !att )
    att = peak_node->first_attribute( "useForDrfFit", 12 );
  if( att )
    m_useForDrfIntrinsicEffFit = compare( att->value(), att->value_size(), "true", 4, false );

  m_useForDrfFwhmFit = PeakDef::sm_defaultUseForDrfFwhmFit;
  att = peak_node->first_attribute( "useForDrfFwhmFit", 16 );
  if( att )
    m_useForDrfFwhmFit = compare( att->value(), att->value_size(), "true", 4, false );

  m_useForDrfDepthOfInteractionFit = PeakDef::sm_defaultUseForDrfDepthOfInteractionFit;
  att = peak_node->first_attribute( "useForDrfDepthOfInteractionFit", 30 );
  if( att )
    m_useForDrfDepthOfInteractionFit = compare( att->value(), att->value_size(), "true", 4, false );

  att = peak_node->first_attribute( "version", 7 );
  if( !att || !att->value_size() )
    throw runtime_error( "missing version attribute" );

  int majorVersion = 0;
  if( sscanf( att->value(), "%i", &majorVersion ) != 1 )
    throw runtime_error( "Non integer version number" );

  // Accept both version 0.x (pre-Voigt) and version 1.x
  if( (majorVersion != sm_peak_xml_major) && (majorVersion != sm_pre_voigt_peak_xml_major) )
    throw runtime_error( "Invalid peak version" );

  xml_node<char> *node = peak_node->first_node( "UserLabel", 9 );
  if( node && node->value() )
    m_userLabel = node->value();

  node = peak_node->first_node( "Type", 4 );
  if( !node || !node->value() )
    throw runtime_error( "No peak type" );

  if( compare( node->value(), node->value_size(), "GaussianDefined", 15, false ) )
    m_type = GaussianDefined;
  else if( compare( node->value(), node->value_size(), "DataDefined", 11, false ) )
    m_type = DataDefined;
  else
    throw runtime_error( "Invalid peak type" );

  node = peak_node->first_node( "Skew", 4 );
  const auto skew_is = [node]( const char *name ){
    return compare( node->value(), node->value_size(), name, strlen(name), false );
  };

  if( !node || !node->value_size() || skew_is("NoSkew") || skew_is("LandauSkew") )
    m_skewType = NoSkew;  //LandauSkew was retired 20231101; it was probably never used
  else if( skew_is("Bortel") || skew_is("ExGauss") )
    m_skewType = Bortel;
  else if( skew_is("CrystalBall") || skew_is("CB") )
    m_skewType = CrystalBall;
  else if( skew_is("DoubleSidedCrystalBall") || skew_is("DSCB") )
    m_skewType = DoubleSidedCrystalBall;
  else if( skew_is("GaussExp") )
    m_skewType = GaussExp;
  else if( skew_is("ExpGaussExp") )
    m_skewType = ExpGaussExp;
  else if( skew_is("VoigtPlusBortel") || skew_is("VoigtWithExpTail") || skew_is("VoigtExp") || skew_is("VoigtBortel") )
    m_skewType = VoigtPlusBortel;
  else if( skew_is("GaussPlusBortel") || skew_is("GaussBortel") )
    m_skewType = GaussPlusBortel;
  else if( skew_is("DoubleBortel") )
    m_skewType = DoubleBortel;
  else if( skew_is("GadrasGeneric") )
    m_skewType = GadrasGeneric;
  else if( skew_is("GadrasCZT") )
    m_skewType = GadrasCZT;
  else
    throw runtime_error( "Invalid peak skew type" );

  const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skewType );
  for( CoefficientType t = CoefficientType(0); t < NumCoefficientTypes; t = CoefficientType(t+1) )
  {
    const bool is_skew = (t >= SkewPar0) && (t <= SkewPar5);
    if( is_skew && (static_cast<size_t>(t - SkewPar0) >= num_skew_pars) )
    {
      m_coefficients[t] = 0.0;
      m_uncertainties[t] = 0.0;
      m_fitFor[t] = false;
      continue;
    }//if( a skew parameter this skew type doesnt have )

    const char *label = to_string( t );
    node = peak_node->first_node( label );
    if( !node || !node->value() )
      throw runtime_error( "No coefficient " + string(label) );

    float value = 0.0f, uncert = 0.0f;
    const int nread = sscanf( node->value(), "%g %g", &value, &uncert );
    if( (nread != 1) && (nread != 2) )
      throw runtime_error( "unable to read value or uncert for " + string(label) );

    m_coefficients[t] = value;
    m_uncertainties[t] = uncert;

    att = node->first_attribute( "fit", 3 );
    if( !att || !att->value() )
    {
      if( t != Chi2DOF )
        throw runtime_error( "No fit attribute for " + string(label) );
      m_fitFor[t] = false;
    }else
    {
      m_fitFor[t] = compare( att->value(), att->value_size(), "true", 4, false );
      if( !m_fitFor[t] && !compare( att->value(), att->value_size(), "false", 5, false ) )
        throw runtime_error( "invalid fit value" );
    }
  }//for( loop over coefficients )

  const xml_node<char> *line_color_node = peak_node->first_node( "LineColor", 9 );
  if( line_color_node && (line_color_node->value_size() >= 7) )
    m_lineColor = Wt::WColor( string( line_color_node->value(), line_color_node->value_size() ) );
  else
    m_lineColor = Wt::WColor();

  const xml_node<char> *nuc_node = peak_node->first_node( "Nuclide", 7 );
  const xml_node<char> *xray_node = peak_node->first_node( "XRay", 4 );
  const xml_node<char> *rctn_node = peak_node->first_node( "Reaction", 8 );

  try
  {
    Source src;

    if( nuc_node )
    {
      const xml_node<char> *name_node = nuc_node->first_node( "Name", 4 );
      const xml_node<char> *p_node = nuc_node->first_node( "DecayParent", 11 );
      const xml_node<char> *c_node = nuc_node->first_node( "DecayChild", 10 );
      const xml_node<char> *e_node = nuc_node->first_node( "DecayGammaEnergy", 16 );
      const xml_node<char> *type_node = nuc_node->first_node( "DecayGammaType", 14 );

      if( !name_node || !name_node->value_size() )
        throw runtime_error( "Missing nuclide Name node" );

      if( !type_node || !gamma_type_from_str( type_node->value(), type_node->value_size(), src.gamma_type ) )
        throw runtime_error( "Missing or invalid DecayGammaType node" );

      src.kind = Source::Kind::Nuclide;
      src.name = name_node->value();

      const bool isNormalNucTrans = (p_node && c_node && e_node && p_node->value_size() && e_node->value_size());
      if( isNormalNucTrans )
      {
        float energy = 0.0f;
        if( (sscanf( e_node->value(), "%g", &energy ) != 1) || !(energy > 0.0f) )
          throw runtime_error( "Invalid nuclide gamma energy" );

        src.decay_parent = (p_node->value() == src.name) ? string() : string( p_node->value() );
        src.decay_child = c_node->value();
        src.particle_energy = energy;
      }//if( isNormalNucTrans )

      if( src.gamma_type == PeakDef::AnnihilationGamma )
      {
        // The DecayGammaEnergy InterSpec writes for annihilation is the positron's energy
        src.particle_energy = sm_annihilation_energy;
      }else if( !isNormalNucTrans )
      {
        // InterSpec doesnt consider this peak to have a source either (no transition)
        src = Source();
      }
    }else if( xray_node )
    {
      const xml_node<char> *el_node = xray_node->first_node( "Element", 7 );
      const xml_node<char> *energy_node = xray_node->first_node( "Energy", 6 );

      if( !el_node || !el_node->value_size() || !energy_node || !energy_node->value() )
        throw runtime_error( "Ill specified xray" );

      float energy = 0.0f;
      if( sscanf( energy_node->value(), "%g", &energy ) != 1 )
        throw runtime_error( "non numeric xray energy" );

      src.kind = Source::Kind::Xray;
      src.name = el_node->value();
      src.particle_energy = energy;
      src.gamma_type = PeakDef::XrayGamma;
    }else if( rctn_node )
    {
      const xml_node<char> *name_node = rctn_node->first_node( "Name", 4 );
      const xml_node<char> *energy_node = rctn_node->first_node( "Energy", 6 );
      const xml_node<char> *type_node = rctn_node->first_node( "Type", 4 );

      if( !name_node || !name_node->value_size() || !energy_node || !energy_node->value() )
        throw runtime_error( "Ill specified reaction" );

      float energy = 0.0f;
      if( sscanf( energy_node->value(), "%g", &energy ) != 1 )
        throw runtime_error( "non numeric reaction energy" );

      src.kind = Source::Kind::Reaction;
      src.name = name_node->value();
      src.particle_energy = energy;
      src.gamma_type = PeakDef::NormalGamma; //early versions didnt write the type
      if( type_node )
        gamma_type_from_str( type_node->value(), type_node->value_size(), src.gamma_type );
    }//if( nuc_node ) / else if( xray_node ) / else if( rctn_node )

    setSource( src );
  }catch( std::exception &e )
  {
    clearSources();

    char msg[256];
    snprintf( msg, sizeof(msg), "Failed to assign peak at %.2f keV to nuclide/xray/reaction: %s",
              mean(), e.what() );
    if( warnings )
      warnings->push_back( msg );
  }//try / catch
}//void fromXml(...)
