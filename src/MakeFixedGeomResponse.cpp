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
#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include <Eigen/Dense>

#include "SpecUtils/RapidXmlUtils.hpp"

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"
#include "efficiency/EfficiencyCalculator.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MakeFixedGeomResponse.h"

using namespace std;

namespace
{
  /** The source region: the innermost non-generic layer that carries a
   self-attenuating or trace source; -1 for a point source at the center.
   */
  int source_layer_index( const MakeFixedGeomResponse::Setup &setup )
  {
    for( size_t i = 0; i < setup.shieldings.size(); ++i )
    {
      const ShieldingSourceFitCalc::ShieldingInfo &info = setup.shieldings[i];
      if( !info.m_traceSources.empty() || !info.m_nuclideFractions_.empty() )
        return static_cast<int>( i );
    }
    return -1;
  }//source_layer_index(...)
}//namespace


std::string MakeFixedGeomResponse::Setup::toXmlString() const
{
  rapidxml::xml_document<char> doc;
  rapidxml::xml_node<char> *base_node
        = doc.allocate_node( rapidxml::node_element, "ActShieldSetup" );
  doc.append_node( base_node );
  base_node->append_attribute( doc.allocate_attribute( "version", "0" ) );

  const char *geom_val = doc.allocate_string( GammaInteractionCalc::to_str(geometry) );
  base_node->append_node( doc.allocate_node( rapidxml::node_element, "Geometry", geom_val ) );

  char buffer[64];
  snprintf( buffer, sizeof(buffer), "%.9g", (distance / PhysicalUnits::cm) );
  const char *dist_val = doc.allocate_string( buffer );
  base_node->append_node( doc.allocate_node( rapidxml::node_element, "DistanceCm", dist_val ) );

  rapidxml::xml_node<char> *shields_node
        = doc.allocate_node( rapidxml::node_element, "Shieldings" );
  base_node->append_node( shields_node );
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : shieldings )
    info.serialize( shields_node );

  string answer;
  rapidxml::print( std::back_inserter(answer), doc, rapidxml::print_no_indenting );
  return answer;
}//Setup::toXmlString()


void MakeFixedGeomResponse::Setup::fromXmlString( const std::string &xml )
{
  geometry = GammaInteractionCalc::GeometryType::Spherical;
  distance = 0.0;
  shieldings.clear();

  vector<char> xml_buf( xml.begin(), xml.end() );
  xml_buf.push_back( '\0' );

  rapidxml::xml_document<char> doc;
  doc.parse<rapidxml::parse_trim_whitespace>( xml_buf.data() );

  const rapidxml::xml_node<char> *base_node = doc.first_node( "ActShieldSetup" );
  if( !base_node )
    throw runtime_error( "Setup::fromXmlString: no ActShieldSetup node" );

  const rapidxml::xml_node<char> *geom_node = base_node->first_node( "Geometry" );
  const string geom_str = SpecUtils::xml_value_str( geom_node );
  bool found_geom = false;
  for( GammaInteractionCalc::GeometryType type :
          { GammaInteractionCalc::GeometryType::Spherical,
            GammaInteractionCalc::GeometryType::CylinderEndOn,
            GammaInteractionCalc::GeometryType::CylinderSideOn,
            GammaInteractionCalc::GeometryType::Rectangular } )
  {
    if( geom_str == GammaInteractionCalc::to_str(type) )
    {
      geometry = type;
      found_geom = true;
    }
  }
  if( !found_geom )
    throw runtime_error( "Setup::fromXmlString: invalid Geometry '" + geom_str + "'" );

  const rapidxml::xml_node<char> *dist_node = base_node->first_node( "DistanceCm" );
  const string dist_str = SpecUtils::xml_value_str( dist_node );
  if( !(stringstream(dist_str) >> distance) )
    throw runtime_error( "Setup::fromXmlString: invalid DistanceCm" );
  distance *= PhysicalUnits::cm;

  const rapidxml::xml_node<char> *shields_node = base_node->first_node( "Shieldings" );
  if( shields_node )
  {
    for( const rapidxml::xml_node<char> *node = shields_node->first_node();
         node; node = node->next_sibling() )
    {
      ShieldingSourceFitCalc::ShieldingInfo info;
      info.deSerialize( node );
      shieldings.push_back( std::move(info) );
    }
  }//if( shields_node )
}//Setup::fromXmlString(...)


bool MakeFixedGeomResponse::sceneRepresentable( const Setup &setup, std::string *reason )
{
  const auto fail = [reason]( const string &why ) -> bool {
    if( reason )
      *reason = why;
    return false;
  };

  if( setup.distance <= 0.0 )
    return fail( "Distance must be positive." );

  const int src_layer = source_layer_index( setup );

  for( size_t i = 0; i < setup.shieldings.size(); ++i )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &info = setup.shieldings[i];

    if( info.m_isGenericMaterial )
      return fail( "Generic (atomic number / areal density) shieldings have no"
                   " physical extent the Monte Carlo can transport through." );

    if( !info.m_material )
      return fail( "A shielding layer has no material defined." );

    const bool is_src = (!info.m_traceSources.empty() || !info.m_nuclideFractions_.empty());
    if( is_src && (static_cast<int>(i) != src_layer) )
      return fail( "Only a single source layer is supported." );

    if( is_src && (i != 0) )
      return fail( "The source must be the innermost layer for the Monte-Carlo"
                   " scene (inner non-source layers are not supported yet)." );

    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      if( (trace.m_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution)
          && (setup.geometry == GammaInteractionCalc::GeometryType::Spherical) )
        return fail( "Exponentially-distributed trace sources are not supported"
                     " for spherical geometry." );
    }
  }//for( shieldings )

  return true;
}//sceneRepresentable(...)


std::shared_ptr<DetectorPeakResponse> MakeFixedGeomResponse::computeFixedGeomDrf(
                  const std::shared_ptr<const DetectorPeakResponse> &base_drf,
                  const Setup &setup,
                  const std::vector<double> &extra_energies_keV,
                  const double fep_precision,
                  const std::function<void(double)> &progress,
                  const std::shared_ptr<std::atomic<bool>> &cancel )
{
  using GammaInteractionCalc::GeometryType;
  using GammaInteractionCalc::TraceActivityType;

  if( !base_drf || !base_drf->isValid() )
    throw runtime_error( "No valid detector response to base the computation on." );

  const shared_ptr<const ceelo::DetectorResponse> mc_resp = base_drf->ceeloResponse();
  if( !mc_resp )
    throw runtime_error( "The detector response has no Monte-Carlo model of the"
        " detector attached; add one via the detector editor first." );

  string why;
  if( !sceneRepresentable( setup, &why ) )
    throw runtime_error( "Scene not representable: " + why );

  // --- Assemble the CeeLo scene: detector from the stored descriptor -------
  ceelo::EfficiencyCalculator calc;
  vector<unique_ptr<ceelo::Material>> owned_mats;
  ceelo::ResponseGenerator::configure_calculator( calc, mc_resp->descriptor, owned_mats );
  calc.set_air_attenuation( ceelo::AirAttenuation::AnalyticNoScatter );

  const auto add_material = [&owned_mats]( const Material &mat ) -> const ceelo::Material * {
    const ceelo::MaterialSpec spec = CeeLoUtils::to_ceelo_material( mat );
    owned_mats.push_back( make_unique<ceelo::Material>( spec.to_material() ) );
    return owned_mats.back().get();
  };

  // z = 0 is the crystal face; DRF distances are from the detector face
  //  (front of the outermost attenuator).
  double front_off_cm = 0.0;
  for( const ceelo::LayerSpec &layer : mc_resp->descriptor.layers )
    front_off_cm += layer.front_thickness_cm;

  const double center_z_cm = -( front_off_cm + (setup.distance / PhysicalUnits::cm) );
  const Eigen::Vector3d center( 0.0, 0.0, center_z_cm );

  const int src_layer = source_layer_index( setup );
  DetectorPeakResponse::EffGeometryType eff_geom_type
                = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;

  if( src_layer < 0 )
  {
    calc.set_point_source( center );
  }else
  {
    const ShieldingSourceFitCalc::ShieldingInfo &src_info = setup.shieldings[src_layer];
    const double dim0_cm = src_info.m_dimensions[0] / PhysicalUnits::cm;
    const double dim1_cm = src_info.m_dimensions[1] / PhysicalUnits::cm;
    const double dim2_cm = src_info.m_dimensions[2] / PhysicalUnits::cm;

    switch( setup.geometry )
    {
      case GeometryType::Spherical:
        calc.set_spherical_source( center, dim0_cm );
        break;

      case GeometryType::CylinderEndOn:
        calc.set_cylindrical_source( center, dim0_cm, dim1_cm );
        break;

      case GeometryType::CylinderSideOn:
      {
        Eigen::Matrix3d rot;
        rot << 0, 0, 1,   0, 1, 0,   -1, 0, 0;  //local z <- detector x
        calc.set_cylindrical_source( center, dim0_cm, dim1_cm, rot );
        break;
      }

      case GeometryType::Rectangular:
        calc.set_rectangular_source( center,
                            Eigen::Vector3d( dim0_cm, dim1_cm, dim2_cm ) );
        break;

      case GeometryType::NumGeometryType:
        throw runtime_error( "Invalid geometry type" );
    }//switch( setup.geometry )

    calc.set_source_material( add_material( *src_info.m_material ) );

    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : src_info.m_traceSources )
    {
      switch( trace.m_type )
      {
        case TraceActivityType::TotalActivity:
        case TraceActivityType::ActivityPerCm3:
          break;  //uniform emission; per-decay efficiency
        case TraceActivityType::ActivityPerGram:
          eff_geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram;
          break;
        case TraceActivityType::ExponentialDistribution:
          eff_geom_type = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;
          calc.set_exponential_depth_distribution( trace.m_relaxationDistance / PhysicalUnits::cm );
          break;
        case TraceActivityType::NumTraceActivityType:
          throw runtime_error( "Invalid trace activity type" );
      }//switch( trace.m_type )
    }//for( trace sources )
  }//if( point source ) / else

  // Shield layers, innermost first, skipping the source layer itself.
  for( size_t i = std::max(0, src_layer + (src_layer >= 0 ? 1 : 0));
       i < setup.shieldings.size(); ++i )
  {
    // For a point source, every layer is a shield; for a volumetric source,
    //  layers after the source.
    if( (src_layer >= 0) && (static_cast<int>(i) <= src_layer) )
      continue;

    const ShieldingSourceFitCalc::ShieldingInfo &info = setup.shieldings[i];
    const ceelo::Material * const mat = add_material( *info.m_material );

    switch( setup.geometry )
    {
      case GeometryType::Spherical:
        calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm );
        break;

      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        if( src_layer >= 0 )
          calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm,
                                  info.m_dimensions[1] / PhysicalUnits::cm );
        else
          calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm );
        break;

      case GeometryType::Rectangular:
        if( src_layer >= 0 )
          calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm,
                                  info.m_dimensions[1] / PhysicalUnits::cm,
                                  info.m_dimensions[2] / PhysicalUnits::cm );
        else
          calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm );
        break;

      case GeometryType::NumGeometryType:
        break;
    }//switch( setup.geometry )
  }//for( shield layers )

  // --- Energy grid: log grid over the DRF's valid range + the fit lines ----
  double e_lo = base_drf->lowerEnergy(), e_hi = base_drf->upperEnergy();
  if( (e_lo <= 10.0) || (e_hi <= e_lo) )
  {
    e_lo = 45.0;
    e_hi = 3000.0;
  }

  set<double> energy_set;
  const int n_grid = 36;
  for( int i = 0; i < n_grid; ++i )
    energy_set.insert( e_lo * std::pow( e_hi/e_lo, static_cast<double>(i)/(n_grid-1) ) );
  for( const double energy : extra_energies_keV )
  {
    if( (energy >= e_lo) && (energy <= e_hi) )
      energy_set.insert( energy );
  }

  const vector<double> energies( begin(energy_set), end(energy_set) );

  // --- Per-energy precision-targeted MC -------------------------------------
  vector<DetectorPeakResponse::EnergyEffPoint> fep_points;
  vector<DetectorPeakResponse::EnergyEfficiencyPair> tot_pairs;

  for( size_t i = 0; i < energies.size(); ++i )
  {
    if( cancel && cancel->load() )
      throw runtime_error( "cancelled" );

    ceelo::SimulationConfig cfg;
    cfg.energy_keV = energies[i];
    cfg.termination.target_fep_rel_precision = fep_precision;
    cfg.termination.max_events = 40000000;
    cfg.termination.max_wall_seconds = 20.0;
    cfg.termination.min_events = 20000;
    cfg.seed = 7000 + i;  //deterministic
    const ceelo::EfficiencyResult res = calc.compute( cfg );

    DetectorPeakResponse::EnergyEffPoint fep;
    fep.energy = static_cast<float>( energies[i] );
    fep.efficiency = static_cast<float>( std::max( 0.0, res.full_energy_peak_efficiency ) );
    if( res.full_energy_peak_efficiency > 0.0 )
      fep.efficiencyUncert = static_cast<float>( res.fep_uncertainty
                                                 / res.full_energy_peak_efficiency );
    fep_points.push_back( fep );

    DetectorPeakResponse::EnergyEfficiencyPair tot;
    tot.energy = static_cast<float>( energies[i] );
    tot.efficiency = static_cast<float>( std::max( 0.0, res.total_efficiency ) );
    tot_pairs.push_back( tot );

    if( progress )
      progress( static_cast<double>(i + 1) / energies.size() );
  }//for( energies )

  // --- Package the DRF -------------------------------------------------------
  auto drf = make_shared<DetectorPeakResponse>( *base_drf );
  drf->setParentHashValue( base_drf->hashValue() );

  {// FEP curve (+ per-point MC uncertainties -> DRF efficiency uncertainty)
    stringstream csv;
    csv << "energy,efficiency\n";
    for( const DetectorPeakResponse::EnergyEffPoint &p : fep_points )
      csv << p.energy << "," << p.efficiency << "\n";
    drf->fromEnergyEfficiencyCsv( csv, base_drf->detectorDiameter(), 0.0,
                                  static_cast<float>(PhysicalUnits::keV), eff_geom_type );
  }

  {
    vector<float> uncert_energies, uncerts;
    for( const DetectorPeakResponse::EnergyEffPoint &p : fep_points )
    {
      if( p.efficiencyUncert.has_value() && (*p.efficiencyUncert > 0.0f) )
      {
        uncert_energies.push_back( p.energy );
        uncerts.push_back( *p.efficiencyUncert );
      }
    }
    if( !uncert_energies.empty() )
      drf->setEfficiencyUncert( DetectorEfficiencyUncert::fromPointUncerts(
                                                uncert_energies, uncerts ) );
  }

  {// Total-efficiency curve (for cascade-summing corrections)
    auto tot_curve = make_shared<DetectorEfficiencyCurve>();
    tot_curve->setFromPairs( tot_pairs, static_cast<float>(PhysicalUnits::keV) );
    drf->setTotalEfficiencyCurve( tot_curve );
  }

  // The far-field MC characterization does not describe this scene - the new
  //  curves do; ditto the raw calibration points.
  drf->setCeeloResponse( nullptr );
  drf->setFixedGeometrySetupXml( setup.toXmlString() );
  drf->setName( base_drf->name() + " (fixed-geom MC)" );

  return drf;
}//computeFixedGeomDrf(...)
