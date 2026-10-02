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
#include <algorithm>

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

  // Always written, so a per-mass/area DRF whose blob lacks it (older, or unknown) is detectable.
  snprintf( buffer, sizeof(buffer), "%.9g", fep_scale );
  const char *scale_val = doc.allocate_string( buffer );
  base_node->append_node( doc.allocate_node( rapidxml::node_element, "FepScale", scale_val ) );

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
  fep_scale = 1.0;
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

  const rapidxml::xml_node<char> *scale_node = base_node->first_node( "FepScale" );
  if( scale_node )
  {
    const string scale_str = SpecUtils::xml_value_str( scale_node );
    if( !(stringstream(scale_str) >> fep_scale) || !(fep_scale > 0.0) || std::isinf(fep_scale) )
      throw runtime_error( "Setup::fromXmlString: invalid FepScale" );
  }//if( scale_node )

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

    // One DRF has one activity convention (total / per gram / per m^2) and CeeLo one emission
    //  distribution, so every source in the layer must share them.
    const auto convention = []( const GammaInteractionCalc::TraceActivityType type ) -> int {
      switch( type )
      {
        case GammaInteractionCalc::TraceActivityType::ActivityPerGram:         return 1;
        case GammaInteractionCalc::TraceActivityType::ExponentialDistribution: return 2;
        default:                                                               return 0;  //total, per cm3
      }
    };

    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      const ShieldingSourceFitCalc::TraceSourceInfo &first = info.m_traceSources.front();  //loop implies non-empty
      if( convention(trace.m_type) != convention(first.m_type) )
        return fail( "All trace sources in the source layer must use the same activity"
                     " convention (total, per gram, or per m^2) for a fixed-geometry response." );
      if( (convention(trace.m_type) == 2) && (trace.m_relaxationDistance != first.m_relaxationDistance) )
        return fail( "All exponentially-distributed trace sources must share one relaxation length"
                     " for a fixed-geometry response." );
      if( !info.m_nuclideFractions_.empty() && (convention(trace.m_type) != 0) )
        return fail( "A self-attenuating source can't share its layer with a per-gram or per-m^2"
                     " trace source for a fixed-geometry response." );
    }//for( trace sources )

    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      // CeeLo's exponential depth profile runs along the source's own axis (local z), which is the
      //  depth InterSpec integrates for an end-on cylinder or a box.  A sphere's and a side-on
      //  cylinder's in-situ depth is RADIAL (eval_spherical / eval_cylinder), which the scene
      //  cannot represent.
      if( trace.m_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
      {
        if( setup.geometry == GammaInteractionCalc::GeometryType::Spherical )
          return fail( "Exponentially-distributed trace sources are not supported"
                       " for spherical geometry." );
        if( setup.geometry == GammaInteractionCalc::GeometryType::CylinderSideOn )
          return fail( "Exponentially-distributed trace sources are not supported"
                       " for side-on cylindrical geometry (the depth profile is radial)." );
      }
    }
  }//for( shieldings )

  return true;
}//sceneRepresentable(...)


std::shared_ptr<DetectorPeakResponse> MakeFixedGeomResponse::computeFixedGeomDrf(
                  const std::shared_ptr<const DetectorPeakResponse> &base_drf,
                  const Setup &setup,
                  const std::vector<double> &extra_energies_keV,
                  const std::vector<double> &partner_energies_keV,
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
  //  For an imported efficiency grid, the user's (possibly edited) description of the detector rather
  //  than the imported geometry the grid is tied to - see DetectorPeakResponse::monteCarloGeometry.
  const shared_ptr<const ceelo::GeometryDescriptor> det_geom = base_drf->monteCarloGeometry();
  const ceelo::GeometryDescriptor &det_gd = det_geom ? *det_geom : mc_resp->descriptor;

  ceelo::EfficiencyCalculator calc;
  vector<unique_ptr<ceelo::Material>> owned_mats;
  ceelo::ResponseGenerator::configure_calculator( calc, det_gd, owned_mats );
  calc.set_air_attenuation( ceelo::AirAttenuation::AnalyticNoScatter );

  const auto add_material = [&owned_mats]( const Material &mat ) -> const ceelo::Material * {
    const ceelo::MaterialSpec spec = CeeLoUtils::to_ceelo_material( mat );
    owned_mats.push_back( make_unique<ceelo::Material>( spec.to_material() ) );
    return owned_mats.back().get();
  };

  // z = 0 is the crystal face; DRF distances are from the detector face
  //  (front of the outermost attenuator).
  const double front_off_cm = det_gd.endcap_front_offset_cm();

  const double center_z_cm = -( front_off_cm + (setup.distance / PhysicalUnits::cm) );
  const Eigen::Vector3d center( 0.0, 0.0, center_z_cm );

  const int src_layer = source_layer_index( setup );
  DetectorPeakResponse::EffGeometryType eff_geom_type
                = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;
  double fep_scale = 1.0;  //see Setup::fep_scale

  // Side-on cylinder: its axis is the detector x axis.
  Eigen::Matrix3d side_on_rot;
  side_on_rot << 0, 0, 1,   0, 1, 0,   -1, 0, 0;  //local z <- detector x

  if( src_layer < 0 )
  {
    // A point source's shields take the geometry's shape.  CeeLo only builds per-axis shells
    //  around an extended source, so use a vanishingly small, non-attenuating one of that shape;
    //  a bare CeeLo point would get spherical shells of dims[0] instead.
    const double tiny_cm = 1.0E-4;
    switch( setup.geometry )
    {
      case GeometryType::Spherical:
        calc.set_point_source( center );
        break;
      case GeometryType::CylinderEndOn:
        calc.set_cylindrical_source( center, tiny_cm, tiny_cm );
        break;
      case GeometryType::CylinderSideOn:
        calc.set_cylindrical_source( center, tiny_cm, tiny_cm, side_on_rot );
        break;
      case GeometryType::Rectangular:
        calc.set_rectangular_source( center, Eigen::Vector3d( tiny_cm, tiny_cm, tiny_cm ) );
        break;
      case GeometryType::NumGeometryType:
        throw runtime_error( "Invalid geometry type" );
    }//switch( setup.geometry )
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
        calc.set_cylindrical_source( center, dim0_cm, dim1_cm, side_on_rot );
        break;

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

    // Per-gram / per-m^2 activities: the FEP curve carries the source mass / emitting area, as
    //  DetectorPeakResponse::convertFixedGeometryType does (the MC itself is per decay).  The
    //  emitting area matches ShieldingSourceChi2Fcn::inSituEmittingArea for the shapes
    //  sceneRepresentable allows an exponential distribution in.
    const double pi = PhysicalUnits::pi;
    const double r = src_info.m_dimensions[0], d1 = src_info.m_dimensions[1], d2 = src_info.m_dimensions[2];
    if( eff_geom_type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram )
    {
      double volume = 0.0;
      switch( setup.geometry )
      {
        case GeometryType::Spherical:      volume = (4.0/3.0)*pi*r*r*r;   break;
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn: volume = pi*r*r*(2.0*d1);      break;  //d1 = half-length
        case GeometryType::Rectangular:    volume = 8.0*r*d1*d2;          break;  //half-dims
        case GeometryType::NumGeometryType: break;
      }
      fep_scale = volume * static_cast<double>(src_info.m_material->density) / PhysicalUnits::gram;
    }else if( eff_geom_type == DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2 )
    {
      const double area = (setup.geometry == GeometryType::Rectangular) ? (2.0*r)*(2.0*d1) : (pi*r*r);
      fep_scale = area / PhysicalUnits::m2;
    }

    if( !(fep_scale > 0.0) || std::isinf(fep_scale) )
      throw runtime_error( "The source has no mass or emitting area." );
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
        calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm,
                                info.m_dimensions[1] / PhysicalUnits::cm );
        break;

      case GeometryType::Rectangular:
        calc.add_source_shield( mat, info.m_dimensions[0] / PhysicalUnits::cm,
                                info.m_dimensions[1] / PhysicalUnits::cm,
                                info.m_dimensions[2] / PhysicalUnits::cm );
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

  // Cascade-summing partners (x-rays especially) can sit outside the DRF's range; the total curve
  //  clamps outside its nodes, so extend the grid to cover them rather than add a node each
  //  (every node is up to ~20 s of MC).
  for( const double energy : partner_energies_keV )
  {
    if( energy >= 15.0 )
      e_lo = std::min( e_lo, energy );
    if( energy <= 3500.0 )
      e_hi = std::max( e_hi, energy );
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

  // --- Correction of the detector model to the detector, k(E) -----------------
  //  A curve-transfer response - an EFFTRAN of a measured curve, or an imported efficiency grid - IS
  //  the measured efficiency; there is no Monte-Carlo model behind it whose grounding says how far
  //  the model is off.  So measure that here: a bare point source on axis at the scene's distance,
  //  in vacuum like the response, and k(E) = response / MC there - on ~10 energies, interpolated
  //  linearly in (ln E, ln k).  Any other response carries its own k(E) in its grounding.
  const bool ground_to_response = (mc_resp->provenance.method == ceelo::ProductionMethod::CurveTransfer);
  const size_t n_k = ground_to_response ? 10 : 0;
  const size_t n_steps = n_k + energies.size();
  vector<double> k_ln_energies, ln_k;

  if( ground_to_response )
  {
    ceelo::EfficiencyCalculator point_calc;
    vector<unique_ptr<ceelo::Material>> point_mats;
    ceelo::ResponseGenerator::configure_calculator( point_calc, det_gd, point_mats );

    // Not closer than the response is characterized (its floor is from the crystal-face origin).
    const double min_face_cm = CeeLoUtils::faceDistanceFromCrystalOrigin( mc_resp->descriptor,
                                                              mc_resp->provenance.min_distance_cm );
    const double d_face_cm = std::max( setup.distance / PhysicalUnits::cm, min_face_cm );
    point_calc.set_point_source( Eigen::Vector3d( 0.0, 0.0, -(front_off_cm + d_face_cm) ) );
    const Eigen::Vector3d resp_pos = CeeLoUtils::sourcePositionFromFace( mc_resp->descriptor,
                                                                         0.0, 0.0, d_face_cm );

    for( size_t i = 0; i < n_k; ++i )
    {
      if( cancel && cancel->load() )
        throw runtime_error( "cancelled" );

      const double energy = e_lo * std::pow( e_hi/e_lo, static_cast<double>(i)/(n_k - 1) );

      ceelo::SimulationConfig cfg;
      cfg.energy_keV = energy;
      cfg.termination.target_fep_rel_precision = fep_precision;
      cfg.termination.max_events = 40000000;
      cfg.termination.max_wall_seconds = 20.0;
      cfg.termination.min_events = 20000;
      cfg.seed = 9000 + i;  //deterministic
      const ceelo::EfficiencyResult res = point_calc.compute( cfg );
      const double resp_eff = mc_resp->eps_fep_at( energy, resp_pos ).value;

      if( (res.full_energy_peak_efficiency > 0.0) && (resp_eff > 0.0) )
      {
        k_ln_energies.push_back( std::log( energy ) );
        ln_k.push_back( std::log( resp_eff / res.full_energy_peak_efficiency ) );
      }

      if( progress )
        progress( static_cast<double>(i + 1) / n_steps );
    }//for( size_t i = 0; i < n_k; ++i )

    if( ln_k.empty() )
      throw runtime_error( "Could not relate the detector model to the detector's efficiency." );
  }//if( ground_to_response )

  const auto k_at = [&k_ln_energies,&ln_k]( const double energy ) -> double {
    const double x = std::log( energy );
    if( x <= k_ln_energies.front() )
      return std::exp( ln_k.front() );
    if( x >= k_ln_energies.back() )
      return std::exp( ln_k.back() );
    const size_t hi = std::upper_bound( begin(k_ln_energies), end(k_ln_energies), x ) - begin(k_ln_energies);
    const double frac = (x - k_ln_energies[hi-1]) / (k_ln_energies[hi] - k_ln_energies[hi-1]);
    return std::exp( ln_k[hi-1] + frac*(ln_k[hi] - ln_k[hi-1]) );
  };//k_at

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

    // The detector model alone is not the detector: the response the DRF carries is grounded to
    //  the measured efficiency by k(E), so this scene's MC must be too, or the fixed-geometry DRF
    //  disagrees with the DRF it was made from by k.  The total is scaled by the same k - a
    //  detector that is more (or less) efficient than its model is so for any deposit too - which
    //  matches the response's eta-table total tier; its other tiers model the total differently.
    double k_ground = 1.0;
    if( ground_to_response )
    {
      k_ground = k_at( energies[i] );
    }else if( !mc_resp->grounding.empty() )
    {
      bool clamped = false;
      k_ground = std::exp( mc_resp->grounding.eval_ln_k( energies[i], clamped ) );
    }

    DetectorPeakResponse::EnergyEffPoint fep;
    fep.energy = static_cast<float>( energies[i] );
    fep.efficiency = static_cast<float>( fep_scale * k_ground * std::max( 0.0, res.full_energy_peak_efficiency ) );
    if( res.full_energy_peak_efficiency > 0.0 )
      fep.efficiencyUncert = static_cast<float>( res.fep_uncertainty
                                                 / res.full_energy_peak_efficiency );
    fep_points.push_back( fep );

    DetectorPeakResponse::EnergyEfficiencyPair tot;
    tot.energy = static_cast<float>( energies[i] );
    tot.efficiency = static_cast<float>( k_ground * std::max( 0.0, res.total_efficiency ) );  //per decay, never fep_scale'd
    tot_pairs.push_back( tot );

    if( progress )
      progress( static_cast<double>(n_k + i + 1) / n_steps );
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
  Setup embedded = setup;
  embedded.fep_scale = fep_scale;
  drf->setFixedGeometrySetupXml( embedded.toXmlString() );
  drf->setName( base_drf->name() + " (fixed-geom MC)" );
  if( base_drf->hasImportedGrid() )
    drf->setDescription( base_drf->description()
                         + "  Fixed-geometry Monte Carlo, corrected to the imported efficiency grid." );

  return drf;
}//computeFixedGeomDrf(...)


double MakeFixedGeomResponse::perDecayFepScale( const DetectorPeakResponse &drf )
{
  if( !drf.isFixedGeometry()
      || (drf.geometryType() == DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct) )
    return 1.0;

  // Only the scale is needed, so read just that element (this runs per energy in some callers, and
  //  a full Setup parse would also look up every material).
  const string &xml = drf.fixedGeometrySetupXml();
  if( xml.empty() )
    return 0.0;

  try
  {
    vector<char> xml_buf( xml.begin(), xml.end() );
    xml_buf.push_back( '\0' );
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( xml_buf.data() );

    const rapidxml::xml_node<char> *base_node = doc.first_node( "ActShieldSetup" );
    const rapidxml::xml_node<char> *scale_node = base_node ? base_node->first_node( "FepScale" ) : nullptr;
    double scale = 0.0;
    if( scale_node && (stringstream( SpecUtils::xml_value_str(scale_node) ) >> scale)
        && (scale > 0.0) && !std::isinf(scale) )
      return scale;
  }catch( std::exception & )
  {
  }

  return 0.0;
}//perDecayFepScale(...)
