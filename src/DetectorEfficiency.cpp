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

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

const double DetectorEfficiencyUncert::sm_defaultLogEnergyCorrLength = 0.5;
const size_t DetectorEfficiencyUncert::sm_maxCovarianceNodes = 100;

namespace
{
  // Mirrors of file-local helpers in DetectorPeakResponse.cpp, for encoding
  //  float arrays into app-url query values.
  string to_url_flt_array( const std::vector<float> &vals, const size_t num_sig_fig )
  {
    string answer;
    for( size_t i = 0; i < vals.size(); ++i )
      answer += (i ? "*" : "") + SpecUtils::printCompact( vals[i], num_sig_fig );
    return answer;
  }//string to_url_flt_array

  vector<float> from_url_flt_array( const string &val )
  {
    vector<float> answer;
    vector<string> parts;
    SpecUtils::split( parts, val, "*" );
    for( const string &v : parts )
    {
      float vf;
      if( !(stringstream(v) >> vf) )
        throw runtime_error( "Invalid float array entry '" + v + "'" );
      answer.push_back( vf );
    }//for( const string &v : parts )

    return answer;
  }//from_url_flt_array(...)


  /** Appends a node containing the floats, printed with enough digits to
   round-trip exactly, as a space separated list.
   */
  void append_float_list_node( ::rapidxml::xml_node<char> *parent,
                               ::rapidxml::xml_document<char> *doc,
                               const char *name,
                               const vector<float> &values )
  {
    char buffer[64];
    stringstream valstrm;
    for( size_t i = 0; i < values.size(); ++i )
    {
      snprintf( buffer, sizeof(buffer), "%1.8E", values[i] );
      valstrm << (i ? " " : "") << buffer;
    }

    const char *val = doc->allocate_string( valstrm.str().c_str() );
    ::rapidxml::xml_node<char> *node = doc->allocate_node( ::rapidxml::node_element,
                                                           doc->allocate_string(name), val );
    parent->append_node( node );
  }//append_float_list_node(...)


  vector<float> parse_float_list_node( const ::rapidxml::xml_node<char> *parent,
                                       const char *name )
  {
    vector<float> answer;
    const ::rapidxml::xml_node<char> *node = parent->first_node( name );
    if( node && node->value_size() )
    {
      if( !SpecUtils::split_to_floats( node->value(), node->value_size(), answer ) )
        throw runtime_error( string("Invalid float list in <") + name + ">" );
    }
    return answer;
  }//parse_float_list_node(...)


  /** Build the (n_req x n_node) interpolation operator L: linear in
   log-energy between bracketing nodes, constant extrapolation beyond the
   ends.  Rows are convex combinations of node values, so L*C*L^T is PSD
   whenever C is.
   */
  vector<double> interpolation_weights( const vector<double> &energies,
                                        const vector<float> &nodes )
  {
    const size_t nreq = energies.size();
    const size_t nnode = nodes.size();
    vector<double> L( nreq * nnode, 0.0 );

    for( size_t i = 0; i < nreq; ++i )
    {
      const double energy = energies[i];

      if( (nnode == 1) || (energy <= nodes.front()) )
      {
        L[i*nnode + 0] = 1.0;
        continue;
      }

      if( energy >= nodes.back() )
      {
        L[i*nnode + (nnode - 1)] = 1.0;
        continue;
      }

      // Find first node strictly greater than energy
      size_t j = 1;
      while( (j < nnode) && (nodes[j] <= energy) )
        ++j;

      assert( (j >= 1) && (j < nnode) );

      const double lnE = std::log( energy );
      const double lnLow = std::log( static_cast<double>(nodes[j-1]) );
      const double lnUp = std::log( static_cast<double>(nodes[j]) );
      const double w = (lnUp > lnLow) ? ((lnE - lnLow) / (lnUp - lnLow)) : 0.0;

      L[i*nnode + (j-1)] = 1.0 - w;
      L[i*nnode + j] = w;
    }//for( size_t i = 0; i < nreq; ++i )

    return L;
  }//interpolation_weights(...)
}//namespace


bool EffUncertBand::operator==( const EffUncertBand &rhs ) const
{
  return (lowerEnergy == rhs.lowerEnergy)
         && (upperEnergy == rhs.upperEnergy)
         && (fractionalUncert == rhs.fractionalUncert);
}//EffUncertBand::operator==


DetectorEfficiencyUncert::DetectorEfficiencyUncert()
  : m_bands{},
    m_covEnergies{},
    m_covMatrix{},
    m_coefCovMatrix{},
    m_corrLength( -1.0 )
{
}//DetectorEfficiencyUncert()


bool DetectorEfficiencyUncert::hasBands() const
{
  return !m_bands.empty();
}


bool DetectorEfficiencyUncert::hasNodeCovariance() const
{
  return !m_covEnergies.empty();
}


bool DetectorEfficiencyUncert::isEmpty() const
{
  return m_bands.empty() && m_covEnergies.empty() && m_coefCovMatrix.empty();
}


vector<double> DetectorEfficiencyUncert::nodeFracCovariance( const vector<double> &energies ) const
{
  const size_t nreq = energies.size();
  vector<double> answer( nreq * nreq, 0.0 );

  const size_t nnode = m_covEnergies.size();
  if( !nnode )
    return answer;

  assert( m_covMatrix.size() == (nnode * nnode) );

  const vector<double> L = interpolation_weights( energies, m_covEnergies );

  // tmp = L * C, size (nreq x nnode)
  vector<double> tmp( nreq * nnode, 0.0 );
  for( size_t i = 0; i < nreq; ++i )
  {
    for( size_t k = 0; k < nnode; ++k )
    {
      const double l_ik = L[i*nnode + k];
      if( l_ik == 0.0 )
        continue;
      for( size_t j = 0; j < nnode; ++j )
        tmp[i*nnode + j] += l_ik * static_cast<double>( m_covMatrix[k*nnode + j] );
    }
  }//for( size_t i = 0; i < nreq; ++i )

  // answer = tmp * L^T, size (nreq x nreq)
  for( size_t i = 0; i < nreq; ++i )
  {
    for( size_t j = 0; j < nreq; ++j )
    {
      double sum = 0.0;
      for( size_t k = 0; k < nnode; ++k )
        sum += tmp[i*nnode + k] * L[j*nnode + k];
      answer[i*nreq + j] = sum;
    }
  }//for( size_t i = 0; i < nreq; ++i )

  return answer;
}//nodeFracCovariance(...)


vector<double> DetectorEfficiencyUncert::efficiencyFracCovariance( const vector<double> &energies ) const
{
  vector<double> answer = nodeFracCovariance( energies );

  if( m_bands.empty() )
    return answer;

  const size_t nreq = energies.size();

  // Index of the band each energy falls in, or -1 if none.
  vector<int> band_index( nreq, -1 );
  for( size_t i = 0; i < nreq; ++i )
  {
    for( size_t b = 0; b < m_bands.size(); ++b )
    {
      if( (energies[i] >= m_bands[b].lowerEnergy) && (energies[i] < m_bands[b].upperEnergy) )
      {
        band_index[i] = static_cast<int>( b );
        break;
      }
    }
  }//for( size_t i = 0; i < nreq; ++i )

  // 100% correlated within a band, uncorrelated between bands.
  for( size_t i = 0; i < nreq; ++i )
  {
    if( band_index[i] < 0 )
      continue;

    const double u = m_bands[band_index[i]].fractionalUncert;
    for( size_t j = 0; j < nreq; ++j )
    {
      if( band_index[j] == band_index[i] )
        answer[i*nreq + j] += u * u;
    }
  }//for( size_t i = 0; i < nreq; ++i )

  return answer;
}//efficiencyFracCovariance(...)


vector<double> DetectorEfficiencyUncert::fracUncertainties( const vector<double> &energies ) const
{
  const size_t nreq = energies.size();
  const vector<double> cov = efficiencyFracCovariance( energies );

  vector<double> answer( nreq, 0.0 );
  for( size_t i = 0; i < nreq; ++i )
    answer[i] = std::sqrt( std::max( 0.0, cov[i*nreq + i] ) );

  return answer;
}//fracUncertainties(...)


shared_ptr<DetectorEfficiencyUncert> DetectorEfficiencyUncert::fromPointUncerts(
                                            const vector<float> &energies,
                                            const vector<float> &fracUncerts,
                                            const double corrLength )
{
  if( energies.size() != fracUncerts.size() )
    throw runtime_error( "DetectorEfficiencyUncert::fromPointUncerts: number of"
                         " energies and uncertainties must match" );

  if( energies.empty() )
    throw runtime_error( "DetectorEfficiencyUncert::fromPointUncerts: no input points" );

  // Sort by energy, removing exact-duplicate energies (keeping the first).
  vector<pair<float,float>> pts;
  pts.reserve( energies.size() );
  for( size_t i = 0; i < energies.size(); ++i )
  {
    if( (energies[i] <= 0.0f) || IsNan(energies[i]) || IsInf(energies[i]) )
      throw runtime_error( "DetectorEfficiencyUncert::fromPointUncerts: energies must be > 0" );
    if( (fracUncerts[i] < 0.0f) || IsNan(fracUncerts[i]) || IsInf(fracUncerts[i]) )
      throw runtime_error( "DetectorEfficiencyUncert::fromPointUncerts: uncertainties must be >= 0" );
    pts.emplace_back( energies[i], fracUncerts[i] );
  }

  std::stable_sort( begin(pts), end(pts),
    []( const pair<float,float> &lhs, const pair<float,float> &rhs ) -> bool {
      return lhs.first < rhs.first;
  } );

  pts.erase( std::unique( begin(pts), end(pts),
    []( const pair<float,float> &lhs, const pair<float,float> &rhs ) -> bool {
      return lhs.first == rhs.first;
  } ), end(pts) );

  if( pts.size() > sm_maxCovarianceNodes )
    throw runtime_error( "DetectorEfficiencyUncert::fromPointUncerts: too many points" );

  const size_t nnode = pts.size();
  vector<float> node_energies( nnode ), cov( nnode * nnode, 0.0f );
  for( size_t i = 0; i < nnode; ++i )
    node_energies[i] = pts[i].first;

  for( size_t j = 0; j < nnode; ++j )
  {
    for( size_t k = j; k < nnode; ++k )
    {
      double rho = (j == k) ? 1.0 : 0.0;
      if( (j != k) && (corrLength > 0.0) )
      {
        const double dlne = std::log( static_cast<double>(pts[j].first) )
                            - std::log( static_cast<double>(pts[k].first) );
        rho = std::exp( -0.5 * std::pow( dlne / corrLength, 2.0 ) );
      }

      const double c_jk = rho * pts[j].second * pts[k].second;
      cov[j*nnode + k] = static_cast<float>( c_jk );
      cov[k*nnode + j] = static_cast<float>( c_jk );
    }//for( size_t k = j; k < nnode; ++k )
  }//for( size_t j = 0; j < nnode; ++j )

  auto answer = make_shared<DetectorEfficiencyUncert>();
  answer->setNodeCovariance( node_energies, cov );
  answer->m_corrLength = (corrLength > 0.0) ? corrLength : -1.0;

  return answer;
}//fromPointUncerts(...)


void DetectorEfficiencyUncert::setBands( const vector<EffUncertBand> &bands )
{
  for( size_t i = 0; i < bands.size(); ++i )
  {
    const EffUncertBand &band = bands[i];

    if( IsNan(band.lowerEnergy) || IsNan(band.upperEnergy) || IsNan(band.fractionalUncert)
        || IsInf(band.lowerEnergy) || IsInf(band.upperEnergy) || IsInf(band.fractionalUncert) )
      throw runtime_error( "DetectorEfficiencyUncert::setBands: NaN/Inf band value" );

    if( band.upperEnergy <= band.lowerEnergy )
      throw runtime_error( "DetectorEfficiencyUncert::setBands: band upper energy must"
                           " be greater than lower energy" );

    if( band.fractionalUncert < 0.0f )
      throw runtime_error( "DetectorEfficiencyUncert::setBands: band uncertainty must be >= 0" );

    if( i && (band.lowerEnergy < bands[i-1].upperEnergy) )
      throw runtime_error( "DetectorEfficiencyUncert::setBands: bands must be sorted"
                           " and non-overlapping" );
  }//for( size_t i = 0; i < bands.size(); ++i )

  m_bands = bands;
}//setBands(...)


void DetectorEfficiencyUncert::setNodeCovariance( const vector<float> &energies,
                                             const vector<float> &covRowMajor )
{
  const size_t nnode = energies.size();

  if( !nnode )
  {
    m_covEnergies.clear();
    m_covMatrix.clear();
    m_corrLength = -1.0;
    return;
  }

  if( nnode > sm_maxCovarianceNodes )
    throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: too many nodes" );

  if( covRowMajor.size() != (nnode * nnode) )
    throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: covariance"
                         " matrix size must be N*N" );

  for( size_t i = 0; i < nnode; ++i )
  {
    if( (energies[i] <= 0.0f) || IsNan(energies[i]) || IsInf(energies[i]) )
      throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: node energies must be > 0" );

    if( i && (energies[i] <= energies[i-1]) )
      throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: node energies"
                           " must be strictly increasing" );
  }//for( size_t i = 0; i < nnode; ++i )

  // Validate approximate symmetry and non-negative diagonal, then symmetrize.
  vector<float> cov = covRowMajor;
  for( size_t i = 0; i < nnode; ++i )
  {
    const float diag = cov[i*nnode + i];
    if( (diag < 0.0f) || IsNan(diag) || IsInf(diag) )
      throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: diagonal"
                           " elements must be >= 0" );

    for( size_t j = i + 1; j < nnode; ++j )
    {
      const float a = cov[i*nnode + j];
      const float b = cov[j*nnode + i];

      if( IsNan(a) || IsInf(a) || IsNan(b) || IsInf(b) )
        throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: NaN/Inf"
                             " covariance element" );

      const float scale = std::max( fabs(a), fabs(b) );
      if( fabs(a - b) > std::max( 1.0E-4f * scale, 1.0E-12f ) )
        throw runtime_error( "DetectorEfficiencyUncert::setNodeCovariance: covariance"
                             " matrix is not symmetric" );

      const float sym = 0.5f * (a + b);
      cov[i*nnode + j] = sym;
      cov[j*nnode + i] = sym;
    }//for( size_t j = i + 1; j < nnode; ++j )
  }//for( size_t i = 0; i < nnode; ++i )

  m_covEnergies = energies;
  m_covMatrix = cov;
  m_corrLength = -1.0;
}//setNodeCovariance(...)


void DetectorEfficiencyUncert::setCoefficientCovariance( const vector<float> &covRowMajor )
{
  if( covRowMajor.empty() )
  {
    m_coefCovMatrix.clear();
    return;
  }

  const size_t n = static_cast<size_t>( std::lround( std::sqrt( static_cast<double>(covRowMajor.size()) ) ) );
  if( (n * n) != covRowMajor.size() )
    throw runtime_error( "DetectorEfficiencyUncert::setCoefficientCovariance:"
                         " matrix must be square" );

  m_coefCovMatrix = covRowMajor;
}//setCoefficientCovariance(...)


const vector<EffUncertBand> &DetectorEfficiencyUncert::bands() const
{
  return m_bands;
}


const vector<float> &DetectorEfficiencyUncert::covarianceEnergies() const
{
  return m_covEnergies;
}


const vector<float> &DetectorEfficiencyUncert::covarianceMatrix() const
{
  return m_covMatrix;
}


const vector<float> &DetectorEfficiencyUncert::coefficientCovariance() const
{
  return m_coefCovMatrix;
}


double DetectorEfficiencyUncert::correlationLength() const
{
  return m_corrLength;
}


void DetectorEfficiencyUncert::toXml( ::rapidxml::xml_node<char> *parent,
                                 ::rapidxml::xml_document<char> *doc ) const
{
  using namespace rapidxml;

  xml_node<char> *base_node = doc->allocate_node( node_element, "EfficiencyUncert" );
  parent->append_node( base_node );

  if( !m_bands.empty() )
  {
    vector<float> values;
    values.reserve( 3 * m_bands.size() );
    for( const EffUncertBand &band : m_bands )
    {
      values.push_back( band.lowerEnergy );
      values.push_back( band.upperEnergy );
      values.push_back( band.fractionalUncert );
    }
    append_float_list_node( base_node, doc, "Bands", values );
  }//if( !m_bands.empty() )

  if( !m_covEnergies.empty() )
  {
    append_float_list_node( base_node, doc, "CovEnergies", m_covEnergies );
    append_float_list_node( base_node, doc, "CovMatrix", m_covMatrix );
  }

  if( !m_coefCovMatrix.empty() )
    append_float_list_node( base_node, doc, "CoefCovMatrix", m_coefCovMatrix );

  if( m_corrLength > 0.0 )
  {
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%1.8E", m_corrLength );
    const char *val = doc->allocate_string( buffer );
    xml_node<char> *node = doc->allocate_node( node_element, "CorrelationLength", val );
    base_node->append_node( node );
  }
}//DetectorEfficiencyUncert::toXml(...)


void DetectorEfficiencyUncert::fromXml( const ::rapidxml::xml_node<char> *node )
{
  using ::rapidxml::internal::compare;

  if( !node )
    throw runtime_error( "DetectorEfficiencyUncert::fromXml: null node" );

  if( !compare( node->name(), node->name_size(), "EfficiencyUncert", 16, false ) )
    throw runtime_error( "DetectorEfficiencyUncert::fromXml: invalid node name" );

  m_bands.clear();
  m_covEnergies.clear();
  m_covMatrix.clear();
  m_coefCovMatrix.clear();
  m_corrLength = -1.0;

  const vector<float> band_vals = parse_float_list_node( node, "Bands" );
  if( !band_vals.empty() )
  {
    if( (band_vals.size() % 3) != 0 )
      throw runtime_error( "DetectorEfficiencyUncert::fromXml: Bands must hold"
                           " (lower, upper, uncert) triples" );

    vector<EffUncertBand> bands;
    for( size_t i = 0; i < band_vals.size(); i += 3 )
    {
      EffUncertBand band;
      band.lowerEnergy = band_vals[i];
      band.upperEnergy = band_vals[i+1];
      band.fractionalUncert = band_vals[i+2];
      bands.push_back( band );
    }
    setBands( bands );
  }//if( !band_vals.empty() )

  const vector<float> cov_energies = parse_float_list_node( node, "CovEnergies" );
  if( !cov_energies.empty() )
  {
    const vector<float> cov_matrix = parse_float_list_node( node, "CovMatrix" );
    setNodeCovariance( cov_energies, cov_matrix );
  }

  const vector<float> coef_cov = parse_float_list_node( node, "CoefCovMatrix" );
  if( !coef_cov.empty() )
    setCoefficientCovariance( coef_cov );

  const ::rapidxml::xml_node<char> *corr_node = node->first_node( "CorrelationLength" );
  if( corr_node && corr_node->value_size() )
  {
    double val;
    if( !SpecUtils::parse_double( corr_node->value(), corr_node->value_size(), val ) )
      throw runtime_error( "DetectorEfficiencyUncert::fromXml: invalid CorrelationLength" );
    m_corrLength = val;
  }
}//DetectorEfficiencyUncert::fromXml(...)


void DetectorEfficiencyUncert::toUrlParts( map<string,string> &parts, const string &prefix ) const
{
  if( !m_bands.empty() )
  {
    vector<float> values;
    values.reserve( 3 * m_bands.size() );
    for( const EffUncertBand &band : m_bands )
    {
      values.push_back( band.lowerEnergy );
      values.push_back( band.upperEnergy );
      values.push_back( band.fractionalUncert );
    }
    parts[prefix + "EFUB"] = to_url_flt_array( values, 5 );
  }//if( !m_bands.empty() )

  if( !m_covEnergies.empty() )
  {
    const size_t nnode = m_covEnergies.size();

    parts[prefix + "EFUE"] = to_url_flt_array( m_covEnergies, 4 );

    // Upper triangle (including diagonal), N*(N+1)/2 values
    vector<float> upper;
    upper.reserve( (nnode * (nnode + 1)) / 2 );
    for( size_t i = 0; i < nnode; ++i )
      for( size_t j = i; j < nnode; ++j )
        upper.push_back( m_covMatrix[i*nnode + j] );

    parts[prefix + "EFUC"] = to_url_flt_array( upper, 5 );

    if( m_corrLength > 0.0 )
      parts[prefix + "EFUL"] = SpecUtils::printCompact( m_corrLength, 5 );
  }//if( !m_covEnergies.empty() )
}//DetectorEfficiencyUncert::toUrlParts(...)


shared_ptr<DetectorEfficiencyUncert> DetectorEfficiencyUncert::fromUrlParts(
                                            const map<string,string> &parts,
                                            const string &prefix )
{
  const auto band_pos = parts.find( prefix + "EFUB" );
  const auto energies_pos = parts.find( prefix + "EFUE" );
  const auto cov_pos = parts.find( prefix + "EFUC" );
  const auto corr_pos = parts.find( prefix + "EFUL" );

  if( (band_pos == end(parts)) && (energies_pos == end(parts)) )
    return nullptr;

  auto answer = make_shared<DetectorEfficiencyUncert>();

  if( band_pos != end(parts) )
  {
    const vector<float> band_vals = from_url_flt_array( band_pos->second );
    if( band_vals.empty() || ((band_vals.size() % 3) != 0) )
      throw runtime_error( "DetectorEfficiencyUncert::fromUrlParts: invalid EFUB" );

    vector<EffUncertBand> bands;
    for( size_t i = 0; i < band_vals.size(); i += 3 )
    {
      EffUncertBand band;
      band.lowerEnergy = band_vals[i];
      band.upperEnergy = band_vals[i+1];
      band.fractionalUncert = band_vals[i+2];
      bands.push_back( band );
    }
    answer->setBands( bands );
  }//if( band_pos != end(parts) )

  if( energies_pos != end(parts) )
  {
    if( cov_pos == end(parts) )
      throw runtime_error( "DetectorEfficiencyUncert::fromUrlParts: EFUE without EFUC" );

    const vector<float> energies = from_url_flt_array( energies_pos->second );
    const vector<float> upper = from_url_flt_array( cov_pos->second );
    const size_t nnode = energies.size();

    if( upper.size() != ((nnode * (nnode + 1)) / 2) )
      throw runtime_error( "DetectorEfficiencyUncert::fromUrlParts: EFUC has wrong"
                           " number of entries" );

    vector<float> cov( nnode * nnode, 0.0f );
    size_t pos = 0;
    for( size_t i = 0; i < nnode; ++i )
    {
      for( size_t j = i; j < nnode; ++j )
      {
        cov[i*nnode + j] = upper[pos];
        cov[j*nnode + i] = upper[pos];
        ++pos;
      }
    }

    answer->setNodeCovariance( energies, cov );

    if( corr_pos != end(parts) )
    {
      double val;
      if( !(stringstream(corr_pos->second) >> val) || (val <= 0.0) )
        throw runtime_error( "DetectorEfficiencyUncert::fromUrlParts: invalid EFUL" );
      answer->m_corrLength = val;
    }
  }//if( energies_pos != end(parts) )

  return answer;
}//DetectorEfficiencyUncert::fromUrlParts(...)


void DetectorEfficiencyUncert::appendToHash( std::size_t &seed ) const
{
  for( const EffUncertBand &band : m_bands )
  {
    boost::hash_combine( seed, band.lowerEnergy );
    boost::hash_combine( seed, band.upperEnergy );
    boost::hash_combine( seed, band.fractionalUncert );
  }

  for( const float val : m_covEnergies )
    boost::hash_combine( seed, val );

  for( const float val : m_covMatrix )
    boost::hash_combine( seed, val );

  for( const float val : m_coefCovMatrix )
    boost::hash_combine( seed, val );

  if( m_corrLength > 0.0 )
    boost::hash_combine( seed, m_corrLength );
}//DetectorEfficiencyUncert::appendToHash(...)


bool DetectorEfficiencyUncert::operator==( const DetectorEfficiencyUncert &rhs ) const
{
  return (m_bands == rhs.m_bands)
         && (m_covEnergies == rhs.m_covEnergies)
         && (m_covMatrix == rhs.m_covMatrix)
         && (m_coefCovMatrix == rhs.m_coefCovMatrix)
         && (m_corrLength == rhs.m_corrLength);
}//DetectorEfficiencyUncert::operator==


#if( PERFORM_DEVELOPER_CHECKS )
namespace
{
  /** Throws if the two float vectors differ by more than small numerical
   rounding (serialization round-trips may be off in the last ULP).
   */
  void check_float_vectors_close( const vector<float> &lhs, const vector<float> &rhs,
                                  const char *what )
  {
    if( lhs.size() != rhs.size() )
      throw runtime_error( string(what) + ": size of LHS ("
                           + std::to_string(lhs.size()) + ") doesnt match RHS ("
                           + std::to_string(rhs.size()) + ")" );

    for( size_t i = 0; i < lhs.size(); ++i )
    {
      const float a = lhs[i], b = rhs[i];
      const float scale = std::max( fabs(a), fabs(b) );
      if( (fabs(a - b) > (1.0E-5f * scale)) && (fabs(a - b) > 1.0E-9f) )
        throw runtime_error( string(what) + ": element " + std::to_string(i)
                             + " of LHS (" + std::to_string(a) + ") doesnt match RHS ("
                             + std::to_string(b) + ")" );
    }
  }//check_float_vectors_close(...)
}//namespace


void DetectorEfficiencyUncert::equalEnough( const DetectorEfficiencyUncert &lhs,
                                       const DetectorEfficiencyUncert &rhs )
{
  if( lhs.m_bands.size() != rhs.m_bands.size() )
    throw runtime_error( "DetectorEfficiencyUncert: number of bands doesnt match" );

  for( size_t i = 0; i < lhs.m_bands.size(); ++i )
  {
    const EffUncertBand &a = lhs.m_bands[i];
    const EffUncertBand &b = rhs.m_bands[i];
    check_float_vectors_close( { a.lowerEnergy, a.upperEnergy, a.fractionalUncert },
                               { b.lowerEnergy, b.upperEnergy, b.fractionalUncert },
                               "DetectorEfficiencyUncert band" );
  }

  check_float_vectors_close( lhs.m_covEnergies, rhs.m_covEnergies,
                             "DetectorEfficiencyUncert covariance energies" );
  check_float_vectors_close( lhs.m_covMatrix, rhs.m_covMatrix,
                             "DetectorEfficiencyUncert covariance matrix" );
  check_float_vectors_close( lhs.m_coefCovMatrix, rhs.m_coefCovMatrix,
                             "DetectorEfficiencyUncert coefficient covariance" );

  const double corr_diff = fabs( lhs.m_corrLength - rhs.m_corrLength );
  const double corr_scale = std::max( fabs(lhs.m_corrLength), fabs(rhs.m_corrLength) );
  if( (corr_diff > (1.0E-5 * corr_scale)) && (corr_diff > 1.0E-9) )
    throw runtime_error( "DetectorEfficiencyUncert: correlation length of LHS ("
                         + std::to_string(lhs.m_corrLength) + ") doesnt match RHS ("
                         + std::to_string(rhs.m_corrLength) + ")" );
}//DetectorEfficiencyUncert::equalEnough(...)
#endif //PERFORM_DEVELOPER_CHECKS


DetectorEfficiencyCurve::DetectorEfficiencyCurve()
  : m_form( DetectorPeakResponse::kNumEfficiencyFnctForms ),
    m_energyUnits( static_cast<float>(PhysicalUnits::keV) ),
    m_energyEfficiencies{},
    m_formula{},
    m_fcn{},
    m_expOfLogCoeffs{},
    m_expOfLogCoeffUncerts{},
    m_uncert{}
{
}//DetectorEfficiencyCurve()


void DetectorEfficiencyCurve::reset()
{
  m_form = DetectorPeakResponse::kNumEfficiencyFnctForms;
  m_energyUnits = static_cast<float>( PhysicalUnits::keV );
  m_energyEfficiencies.clear();
  m_formula.clear();
  m_fcn = std::function<float(float)>();
  m_expOfLogCoeffs.clear();
  m_expOfLogCoeffUncerts.clear();
  m_uncert.reset();
}//DetectorEfficiencyCurve::reset()


bool DetectorEfficiencyCurve::isValid() const
{
  switch( m_form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
      return (m_energyEfficiencies.size() >= 2);
    case DetectorPeakResponse::kFunctialEfficienyForm:
      return !!m_fcn;
    case DetectorPeakResponse::kExpOfLogPowerSeries:
      return !m_expOfLogCoeffs.empty();
    case DetectorPeakResponse::kNumEfficiencyFnctForms:
      break;
  }//switch( m_form )

  return false;
}//DetectorEfficiencyCurve::isValid()


float DetectorEfficiencyCurve::efficiency( const float energy ) const
{
  const float x = energy / m_energyUnits;

  switch( m_form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
      if( m_energyEfficiencies.size() < 2 )
        break;
      return DetectorPeakResponse::akimaInterpolate( x, m_energyEfficiencies );

    case DetectorPeakResponse::kFunctialEfficienyForm:
      if( !m_fcn )
        break;
      return m_fcn( x );

    case DetectorPeakResponse::kExpOfLogPowerSeries:
      if( m_expOfLogCoeffs.empty() )
        break;
      return DetectorPeakResponse::expOfLogPowerSeriesEfficiency( x, m_expOfLogCoeffs );

    case DetectorPeakResponse::kNumEfficiencyFnctForms:
      break;
  }//switch( m_form )

  throw runtime_error( "DetectorEfficiencyCurve::efficiency: curve has not been initialized" );
}//DetectorEfficiencyCurve::efficiency(...)


void DetectorEfficiencyCurve::setFromPairs(
                  const vector<DetectorPeakResponse::EnergyEfficiencyPair> &pairs,
                  const float energyUnits )
{
  if( pairs.size() < 2 )
    throw runtime_error( "DetectorEfficiencyCurve::setFromPairs: need at least two points" );

  if( (energyUnits <= 0.0f) || IsNan(energyUnits) || IsInf(energyUnits) )
    throw runtime_error( "DetectorEfficiencyCurve::setFromPairs: invalid energy units" );

  vector<DetectorPeakResponse::EnergyEfficiencyPair> sorted_pairs = pairs;
  std::sort( begin(sorted_pairs), end(sorted_pairs) );

  for( const DetectorPeakResponse::EnergyEfficiencyPair &p : sorted_pairs )
  {
    if( (p.energy < 0.0f) || (p.efficiency < 0.0f)
        || IsNan(p.energy) || IsInf(p.energy)
        || IsNan(p.efficiency) || IsInf(p.efficiency) )
      throw runtime_error( "DetectorEfficiencyCurve::setFromPairs: energies and"
                           " efficiencies must be >= 0" );
  }

  m_form = DetectorPeakResponse::kEnergyEfficiencyPairs;
  m_energyUnits = energyUnits;
  m_energyEfficiencies = sorted_pairs;
  m_formula.clear();
  m_fcn = std::function<float(float)>();
  m_expOfLogCoeffs.clear();
  m_expOfLogCoeffUncerts.clear();
}//DetectorEfficiencyCurve::setFromPairs(...)


void DetectorEfficiencyCurve::setFromFormula( const string &formula, const float energyUnits )
{
  if( (energyUnits <= 0.0f) || IsNan(energyUnits) || IsInf(energyUnits) )
    throw runtime_error( "DetectorEfficiencyCurve::setFromFormula: invalid energy units" );

  const bool isMeV = (energyUnits > 10.0f);

  // Throws if the formula is invalid
  std::function<float(float)> fcn
            = DetectorPeakResponse::makeEfficiencyFunctionFromFormula( formula, isMeV );

  m_form = DetectorPeakResponse::kFunctialEfficienyForm;
  m_energyUnits = energyUnits;
  m_energyEfficiencies.clear();
  m_formula = formula;
  m_fcn = fcn;
  m_expOfLogCoeffs.clear();
  m_expOfLogCoeffUncerts.clear();
}//DetectorEfficiencyCurve::setFromFormula(...)


void DetectorEfficiencyCurve::setFromExpOfLogPowerSeries( const vector<float> &coefs,
                                                     const vector<float> &uncerts,
                                                     const float energyUnits )
{
  if( coefs.empty() )
    throw runtime_error( "DetectorEfficiencyCurve::setFromExpOfLogPowerSeries: no coefficients" );

  if( (energyUnits <= 0.0f) || IsNan(energyUnits) || IsInf(energyUnits) )
    throw runtime_error( "DetectorEfficiencyCurve::setFromExpOfLogPowerSeries: invalid energy units" );

  if( !uncerts.empty() && (uncerts.size() != coefs.size()) )
    throw runtime_error( "DetectorEfficiencyCurve::setFromExpOfLogPowerSeries: number of"
                         " uncertainties must be empty or match the coefficients" );

  m_form = DetectorPeakResponse::kExpOfLogPowerSeries;
  m_energyUnits = energyUnits;
  m_energyEfficiencies.clear();
  m_formula.clear();
  m_fcn = std::function<float(float)>();
  m_expOfLogCoeffs = coefs;
  m_expOfLogCoeffUncerts = uncerts;
}//DetectorEfficiencyCurve::setFromExpOfLogPowerSeries(...)


void DetectorEfficiencyCurve::setRawFields(
                  const DetectorPeakResponse::EfficiencyFnctForm form,
                  const float energyUnits,
                  const vector<DetectorPeakResponse::EnergyEfficiencyPair> &pairs,
                  const string &formula,
                  const vector<float> &coefs,
                  const vector<float> &coefUncerts )
{
  // Preserve the rich uncertainty (serialized separately); only the
  //  representation fields are (re)set here.
  switch( form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
      setFromPairs( pairs, energyUnits );
      break;

    case DetectorPeakResponse::kFunctialEfficienyForm:
      setFromFormula( formula, energyUnits );
      break;

    case DetectorPeakResponse::kExpOfLogPowerSeries:
      setFromExpOfLogPowerSeries( coefs, coefUncerts, energyUnits );
      break;

    case DetectorPeakResponse::kNumEfficiencyFnctForms:
    {
      // An invalid/empty efficiency - keep the units, clear the rest.
      std::shared_ptr<const DetectorEfficiencyUncert> uncert = m_uncert;
      reset();
      m_energyUnits = ((energyUnits > 0.0f) && !IsNan(energyUnits) && !IsInf(energyUnits))
                          ? energyUnits : static_cast<float>(PhysicalUnits::keV);
      m_uncert = uncert;
      break;
    }
  }//switch( form )
}//DetectorEfficiencyCurve::setRawFields(...)


DetectorPeakResponse::EfficiencyFnctForm DetectorEfficiencyCurve::form() const
{
  return m_form;
}


float DetectorEfficiencyCurve::energyUnits() const
{
  return m_energyUnits;
}


const vector<DetectorPeakResponse::EnergyEfficiencyPair> &DetectorEfficiencyCurve::energyEfficiencies() const
{
  return m_energyEfficiencies;
}


const string &DetectorEfficiencyCurve::formula() const
{
  return m_formula;
}


const std::function<float(float)> &DetectorEfficiencyCurve::efficiencyFcn() const
{
  return m_fcn;
}


const vector<float> &DetectorEfficiencyCurve::expOfLogPowerSeriesCoeffs() const
{
  return m_expOfLogCoeffs;
}


const vector<float> &DetectorEfficiencyCurve::expOfLogPowerSeriesUncerts() const
{
  return m_expOfLogCoeffUncerts;
}


shared_ptr<const DetectorEfficiencyUncert> DetectorEfficiencyCurve::uncertainty() const
{
  return m_uncert;
}


void DetectorEfficiencyCurve::setUncertainty( shared_ptr<const DetectorEfficiencyUncert> uncert )
{
  m_uncert = uncert;
}


DetectorEfficiencyCurve DetectorEfficiencyCurve::scaledByConstant( const double factor ) const
{
  if( !isValid() )
    throw runtime_error( "DetectorEfficiencyCurve::scaledByConstant: invalid curve" );

  if( (factor <= 0.0) || IsNan(factor) || IsInf(factor) )
    throw runtime_error( "DetectorEfficiencyCurve::scaledByConstant: factor must be > 0" );

  DetectorEfficiencyCurve answer( *this );

  switch( m_form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
      for( DetectorPeakResponse::EnergyEfficiencyPair &p : answer.m_energyEfficiencies )
        p.efficiency = static_cast<float>( factor * p.efficiency );
      break;

    case DetectorPeakResponse::kFunctialEfficienyForm:
    {
      // Wrap the parsed function; leave the stored formula string unchanged
      //  (matches the legacy DetectorPeakResponse::convertFixedGeometryType).
      const std::function<float(float)> old_fcn = m_fcn;
      answer.m_fcn = [old_fcn, factor]( float energy ) -> float {
        return static_cast<float>( factor * old_fcn(energy) );
      };
      break;
    }

    case DetectorPeakResponse::kExpOfLogPowerSeries:
      answer.m_expOfLogCoeffs[0] += static_cast<float>( std::log(factor) );
      break;

    case DetectorPeakResponse::kNumEfficiencyFnctForms:
      break;
  }//switch( m_form )

  return answer;
}//DetectorEfficiencyCurve::scaledByConstant(...)


void DetectorEfficiencyCurve::toXml( ::rapidxml::xml_node<char> *parent,
                                ::rapidxml::xml_document<char> *doc,
                                const char *node_name ) const
{
  using namespace rapidxml;

  xml_node<char> *base_node = doc->allocate_node( node_element,
                                                  doc->allocate_string(node_name) );
  parent->append_node( base_node );

  const char *form_str = "Undefined";
  switch( m_form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs: form_str = "EnergyEfficiencyPairs"; break;
    case DetectorPeakResponse::kFunctialEfficienyForm: form_str = "FunctialEfficienyForm"; break;
    case DetectorPeakResponse::kExpOfLogPowerSeries:   form_str = "ExpOfLogPowerSeries";   break;
    case DetectorPeakResponse::kNumEfficiencyFnctForms:                                    break;
  }//switch( m_form )

  xml_node<char> *node = doc->allocate_node( node_element, "EfficiencyForm", form_str );
  base_node->append_node( node );

  char buffer[64];
  snprintf( buffer, sizeof(buffer), "%1.8E", m_energyUnits );
  const char *val = doc->allocate_string( buffer );
  node = doc->allocate_node( node_element, "EfficiencyEnergyUnits", val );
  base_node->append_node( node );

  if( !m_energyEfficiencies.empty() )
  {
    vector<float> values;
    values.reserve( 2 * m_energyEfficiencies.size() );
    for( const DetectorPeakResponse::EnergyEfficiencyPair &p : m_energyEfficiencies )
    {
      values.push_back( p.energy );
      values.push_back( p.efficiency );
    }
    append_float_list_node( base_node, doc, "EnergyEfficiencies", values );
  }//if( !m_energyEfficiencies.empty() )

  if( !m_formula.empty() )
  {
    val = doc->allocate_string( m_formula.c_str() );
    node = doc->allocate_node( node_element, "EfficiencyFormula", val );
    base_node->append_node( node );
  }

  if( !m_expOfLogCoeffs.empty() )
    append_float_list_node( base_node, doc, "ExpOfLogPowerSeriesCoeffs", m_expOfLogCoeffs );

  if( !m_expOfLogCoeffUncerts.empty() )
    append_float_list_node( base_node, doc, "ExpOfLogPowerSeriesUncerts", m_expOfLogCoeffUncerts );

  if( m_uncert && !m_uncert->isEmpty() )
    m_uncert->toXml( base_node, doc );
}//DetectorEfficiencyCurve::toXml(...)


void DetectorEfficiencyCurve::fromXml( const ::rapidxml::xml_node<char> *node )
{
  using ::rapidxml::internal::compare;

  if( !node )
    throw runtime_error( "DetectorEfficiencyCurve::fromXml: null node" );

  m_form = DetectorPeakResponse::kNumEfficiencyFnctForms;
  m_energyUnits = static_cast<float>( PhysicalUnits::keV );
  m_energyEfficiencies.clear();
  m_formula.clear();
  m_fcn = std::function<float(float)>();
  m_expOfLogCoeffs.clear();
  m_expOfLogCoeffUncerts.clear();
  m_uncert.reset();

  const ::rapidxml::xml_node<char> *form_node = node->first_node( "EfficiencyForm" );
  if( !form_node || !form_node->value_size() )
    throw runtime_error( "DetectorEfficiencyCurve::fromXml: missing EfficiencyForm" );

  DetectorPeakResponse::EfficiencyFnctForm form;
  if( compare(form_node->value(), form_node->value_size(), "EnergyEfficiencyPairs", 21, false) )
    form = DetectorPeakResponse::kEnergyEfficiencyPairs;
  else if( compare(form_node->value(), form_node->value_size(), "FunctialEfficienyForm", 21, false) )
    form = DetectorPeakResponse::kFunctialEfficienyForm;
  else if( compare(form_node->value(), form_node->value_size(), "ExpOfLogPowerSeries", 19, false) )
    form = DetectorPeakResponse::kExpOfLogPowerSeries;
  else
    throw runtime_error( "DetectorEfficiencyCurve::fromXml: invalid EfficiencyForm" );

  float energy_units = static_cast<float>( PhysicalUnits::keV );
  const ::rapidxml::xml_node<char> *units_node = node->first_node( "EfficiencyEnergyUnits" );
  if( units_node && units_node->value_size() )
  {
    if( !SpecUtils::parse_float( units_node->value(), units_node->value_size(), energy_units ) )
      throw runtime_error( "DetectorEfficiencyCurve::fromXml: invalid EfficiencyEnergyUnits" );
  }

  switch( form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
    {
      const vector<float> values = parse_float_list_node( node, "EnergyEfficiencies" );
      if( values.size() < 4 || (values.size() % 2) != 0 )
        throw runtime_error( "DetectorEfficiencyCurve::fromXml: invalid EnergyEfficiencies" );

      vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
      for( size_t i = 0; i < values.size(); i += 2 )
      {
        DetectorPeakResponse::EnergyEfficiencyPair p;
        p.energy = values[i];
        p.efficiency = values[i+1];
        pairs.push_back( p );
      }
      setFromPairs( pairs, energy_units );
      break;
    }//case kEnergyEfficiencyPairs

    case DetectorPeakResponse::kFunctialEfficienyForm:
    {
      const ::rapidxml::xml_node<char> *formula_node = node->first_node( "EfficiencyFormula" );
      if( !formula_node || !formula_node->value_size() )
        throw runtime_error( "DetectorEfficiencyCurve::fromXml: missing EfficiencyFormula" );
      setFromFormula( string(formula_node->value(), formula_node->value() + formula_node->value_size()),
                      energy_units );
      break;
    }//case kFunctialEfficienyForm

    case DetectorPeakResponse::kExpOfLogPowerSeries:
    {
      const vector<float> coefs = parse_float_list_node( node, "ExpOfLogPowerSeriesCoeffs" );
      if( coefs.empty() || (coefs.size() > 25) )
        throw runtime_error( "DetectorEfficiencyCurve::fromXml: invalid ExpOfLogPowerSeriesCoeffs" );

      vector<float> uncerts = parse_float_list_node( node, "ExpOfLogPowerSeriesUncerts" );
      if( !uncerts.empty() && (uncerts.size() != coefs.size()) )
        throw runtime_error( "DetectorEfficiencyCurve::fromXml: ExpOfLogPowerSeriesUncerts"
                             " size doesnt match coefficients" );

      setFromExpOfLogPowerSeries( coefs, uncerts, energy_units );
      break;
    }//case kExpOfLogPowerSeries

    case DetectorPeakResponse::kNumEfficiencyFnctForms:
      break;
  }//switch( form )

  const ::rapidxml::xml_node<char> *uncert_node = node->first_node( "EfficiencyUncert" );
  if( uncert_node )
  {
    auto uncert = make_shared<DetectorEfficiencyUncert>();
    uncert->fromXml( uncert_node );
    m_uncert = uncert;
  }
}//DetectorEfficiencyCurve::fromXml(...)


void DetectorEfficiencyCurve::toUrlParts( map<string,string> &parts, const string &prefix ) const
{
  switch( m_form )
  {
    case DetectorPeakResponse::kEnergyEfficiencyPairs:
    {
      parts[prefix + "EFT"] = "P";

      vector<float> energies, effs;
      for( const DetectorPeakResponse::EnergyEfficiencyPair &p : m_energyEfficiencies )
      {
        energies.push_back( p.energy );
        effs.push_back( p.efficiency );
      }
      parts[prefix + "EFX"] = to_url_flt_array( energies, 4 );
      parts[prefix + "EFY"] = to_url_flt_array( effs, 5 );
      break;
    }//case kEnergyEfficiencyPairs

    case DetectorPeakResponse::kFunctialEfficienyForm:
      parts[prefix + "EFT"] = "F";
      // Note: caller (DetectorPeakResponse::toAppUrl) url-encodes all values
      //  for keys ending in "EFE".
      parts[prefix + "EFE"] = m_formula;
      break;

    case DetectorPeakResponse::kExpOfLogPowerSeries:
      parts[prefix + "EFT"] = "E";
      parts[prefix + "EFC"] = to_url_flt_array( m_expOfLogCoeffs, 7 );
      if( !m_expOfLogCoeffUncerts.empty() )
        parts[prefix + "EFU"] = to_url_flt_array( m_expOfLogCoeffUncerts, 7 );
      break;

    case DetectorPeakResponse::kNumEfficiencyFnctForms:
      return;
  }//switch( m_form )

  const bool notKeV = (fabs(m_energyUnits - PhysicalUnits::keV) > 0.001);
  if( notKeV )
    parts[prefix + "EUNIT"] = SpecUtils::printCompact( m_energyUnits, 4 );

  if( m_uncert && !m_uncert->isEmpty() )
    m_uncert->toUrlParts( parts, prefix );
}//DetectorEfficiencyCurve::toUrlParts(...)


shared_ptr<DetectorEfficiencyCurve> DetectorEfficiencyCurve::fromUrlParts(
                                            const map<string,string> &parts,
                                            const string &prefix )
{
  const auto form_pos = parts.find( prefix + "EFT" );
  if( form_pos == end(parts) )
    return nullptr;

  float energy_units = static_cast<float>( PhysicalUnits::keV );
  const auto units_pos = parts.find( prefix + "EUNIT" );
  if( units_pos != end(parts) )
  {
    if( !(stringstream(units_pos->second) >> energy_units) || (energy_units <= 0.0f) )
      throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: invalid "
                           + prefix + "EUNIT" );
  }

  auto answer = make_shared<DetectorEfficiencyCurve>();

  string form_str = form_pos->second;
  SpecUtils::to_upper_ascii( form_str );

  if( form_str == "P" )
  {
    const auto x_pos = parts.find( prefix + "EFX" );
    const auto y_pos = parts.find( prefix + "EFY" );
    if( (x_pos == end(parts)) || (y_pos == end(parts)) )
      throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: missing "
                           + prefix + "EFX/" + prefix + "EFY" );

    const vector<float> x = from_url_flt_array( x_pos->second );
    const vector<float> y = from_url_flt_array( y_pos->second );

    if( (x.size() < 2) || (x.size() != y.size()) )
      throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: invalid "
                           + prefix + "EFX/" + prefix + "EFY" );

    vector<DetectorPeakResponse::EnergyEfficiencyPair> pairs;
    for( size_t i = 0; i < x.size(); ++i )
    {
      DetectorPeakResponse::EnergyEfficiencyPair p;
      p.energy = x[i];
      p.efficiency = y[i];
      pairs.push_back( p );
    }
    answer->setFromPairs( pairs, energy_units );
  }else if( form_str == "F" )
  {
    const auto eqn_pos = parts.find( prefix + "EFE" );
    if( eqn_pos == end(parts) )
      throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: missing " + prefix + "EFE" );
    answer->setFromFormula( eqn_pos->second, energy_units );
  }else if( form_str == "E" )
  {
    const auto coef_pos = parts.find( prefix + "EFC" );
    if( coef_pos == end(parts) )
      throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: missing " + prefix + "EFC" );

    const vector<float> coefs = from_url_flt_array( coef_pos->second );
    if( coefs.empty() )
      throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: invalid " + prefix + "EFC" );

    vector<float> uncerts;
    const auto uncert_pos = parts.find( prefix + "EFU" );
    if( uncert_pos != end(parts) )
    {
      uncerts = from_url_flt_array( uncert_pos->second );
      if( !uncerts.empty() && (uncerts.size() != coefs.size()) )
        throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: " + prefix
                             + "EFC and " + prefix + "EFU different lengths" );
    }

    answer->setFromExpOfLogPowerSeries( coefs, uncerts, energy_units );
  }else
  {
    throw runtime_error( "DetectorEfficiencyCurve::fromUrlParts: invalid "
                         + prefix + "EFT value '" + form_pos->second + "'" );
  }

  answer->m_uncert = DetectorEfficiencyUncert::fromUrlParts( parts, prefix );

  return answer;
}//DetectorEfficiencyCurve::fromUrlParts(...)


void DetectorEfficiencyCurve::appendToHash( std::size_t &seed ) const
{
  boost::hash_combine( seed, static_cast<int>(m_form) );
  boost::hash_combine( seed, m_energyUnits );

  for( const DetectorPeakResponse::EnergyEfficiencyPair &p : m_energyEfficiencies )
  {
    boost::hash_combine( seed, p.energy );
    boost::hash_combine( seed, p.efficiency );
  }

  boost::hash_combine( seed, m_formula );

  for( const float val : m_expOfLogCoeffs )
    boost::hash_combine( seed, val );

  for( const float val : m_expOfLogCoeffUncerts )
    boost::hash_combine( seed, val );

  if( m_uncert )
    m_uncert->appendToHash( seed );
}//DetectorEfficiencyCurve::appendToHash(...)


bool DetectorEfficiencyCurve::operator==( const DetectorEfficiencyCurve &rhs ) const
{
  return (m_form == rhs.m_form)
         && (m_energyUnits == rhs.m_energyUnits)
         && (m_energyEfficiencies == rhs.m_energyEfficiencies)
         && (m_formula == rhs.m_formula)
         && ((!m_fcn) == (!rhs.m_fcn))
         && (m_expOfLogCoeffs == rhs.m_expOfLogCoeffs)
         && (m_expOfLogCoeffUncerts == rhs.m_expOfLogCoeffUncerts)
         && ((!m_uncert && !rhs.m_uncert)
             || (m_uncert && rhs.m_uncert && (*m_uncert == *rhs.m_uncert)));
}//DetectorEfficiencyCurve::operator==


#if( PERFORM_DEVELOPER_CHECKS )
void DetectorEfficiencyCurve::equalEnough( const DetectorEfficiencyCurve &lhs,
                                      const DetectorEfficiencyCurve &rhs )
{
  if( lhs.m_form != rhs.m_form )
    throw runtime_error( "DetectorEfficiencyCurve: efficiency form of LHS ("
                         + std::to_string(int(lhs.m_form)) + ") doesnt match RHS ("
                         + std::to_string(int(rhs.m_form)) + ")" );

  check_float_vectors_close( { lhs.m_energyUnits }, { rhs.m_energyUnits },
                             "DetectorEfficiencyCurve energy units" );

  if( lhs.m_energyEfficiencies.size() != rhs.m_energyEfficiencies.size() )
    throw runtime_error( "DetectorEfficiencyCurve: number of energy/efficiency pairs doesnt match" );

  for( size_t i = 0; i < lhs.m_energyEfficiencies.size(); ++i )
  {
    check_float_vectors_close(
        { lhs.m_energyEfficiencies[i].energy, lhs.m_energyEfficiencies[i].efficiency },
        { rhs.m_energyEfficiencies[i].energy, rhs.m_energyEfficiencies[i].efficiency },
        "DetectorEfficiencyCurve energy/efficiency pair" );
  }

  if( lhs.m_formula != rhs.m_formula )
    throw runtime_error( "DetectorEfficiencyCurve: formula of LHS ('" + lhs.m_formula
                         + "') doesnt match RHS ('" + rhs.m_formula + "')" );

  if( (!lhs.m_fcn) != (!rhs.m_fcn) )
    throw runtime_error( "DetectorEfficiencyCurve: availability of parsed formula doesnt match" );

  check_float_vectors_close( lhs.m_expOfLogCoeffs, rhs.m_expOfLogCoeffs,
                             "DetectorEfficiencyCurve exp-of-log coefficients" );

  check_float_vectors_close( lhs.m_expOfLogCoeffUncerts, rhs.m_expOfLogCoeffUncerts,
                             "DetectorEfficiencyCurve exp-of-log coefficient uncertainties" );

  if( (!lhs.m_uncert) != (!rhs.m_uncert) )
    throw runtime_error( "DetectorEfficiencyCurve: availability of uncertainty doesnt match" );

  if( lhs.m_uncert && rhs.m_uncert )
    DetectorEfficiencyUncert::equalEnough( *lhs.m_uncert, *rhs.m_uncert );
}//DetectorEfficiencyCurve::equalEnough(...)
#endif //PERFORM_DEVELOPER_CHECKS
