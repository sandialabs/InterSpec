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
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <Eigen/Core>

#include "SandiaDecay.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/AngleOutxImport.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;


namespace CeeLoUtils
{

ceelo::MaterialSpec to_ceelo_material( const Material &mat )
{
  ceelo::MaterialSpec spec;
  spec.name = mat.name;
  spec.density_g_per_cm3 = mat.density * PhysicalUnits::cm3 / PhysicalUnits::g;

  map<int,double> frac_by_z;
  for( const Material::ElementFractionPair &efp : mat.elements )
    frac_by_z[efp.first->atomicNumber] += efp.second;
  for( const Material::NuclideFractionPair &nfp : mat.nuclides )
    frac_by_z[nfp.first->atomicNumber] += nfp.second;

  double sum = 0.0;
  for( const auto &zf : frac_by_z )
    sum += zf.second;
  if( sum <= 0.0 )
    throw runtime_error( "Material '" + mat.name + "' has no composition." );

  for( const auto &zf : frac_by_z )
  {
    if( (zf.first < 1) || (zf.first > 92) )
      throw runtime_error( "Material '" + mat.name + "' contains an element"
                           " (Z=" + std::to_string(zf.first) + ") the Monte"
                           " Carlo has no cross-section data for." );
    ceelo::MaterialComponent c;
    c.Z = static_cast<uint8_t>( zf.first );
    c.mass_fraction = zf.second / sum;
    spec.composition.push_back( c );
  }

  return spec;
}//to_ceelo_material(...)


TransferAnchor transferAnchorForDrf(
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const ceelo::GeometryDescriptor &geom,
                      const double override_ref_distance_cm )
{
  if( !drf || !drf->isValid() )
    throw runtime_error( "transferAnchorForDrf: invalid detector response." );

  if( drf->isFixedGeometry() )
    throw runtime_error( "transferAnchorForDrf: a fixed-geometry efficiency"
                         " has no source-detector geometry to transfer." );

  TransferAnchor answer;

  // Preferred: the raw measured points, when they were all taken at one
  //  distance - they are the data the curve was fit to (design doc: anchor on
  //  raw points, not a pre-fit curve, when possible).
  const shared_ptr<const MeasuredDrfPoints> raw = drf->measuredPoints();
  if( raw && !raw->empty() )
  {
    vector<MeasuredEffPoint> pts;
    for( const MeasuredEffPoint &p : raw->points() )
    {
      if( (p.distance > 0.0f) && (p.efficiency > 0.0f) && (p.energy > 0.0f) )
        pts.push_back( p );
    }

    double min_d = std::numeric_limits<double>::max(), max_d = 0.0, sum_d = 0.0;
    for( const MeasuredEffPoint &p : pts )
    {
      min_d = std::min( min_d, double(p.distance) );
      max_d = std::max( max_d, double(p.distance) );
      sum_d += p.distance;
    }

    if( (pts.size() >= 2) && (max_d < 1.01*min_d) )
    {
      sort( begin(pts), end(pts),
            []( const MeasuredEffPoint &a, const MeasuredEffPoint &b ){
              return a.energy < b.energy;
            } );

      // Merge points at (nearly) the same energy - the anchor needs strictly
      //  ascending energies.  Inverse-variance weight when uncertainties are
      //  available, plain average otherwise.
      for( size_t i = 0; i < pts.size(); /* advanced inside */ )
      {
        size_t j = i + 1;
        while( (j < pts.size()) && (pts[j].energy < pts[i].energy + 0.01f) )
          ++j;

        double sum_w = 0.0, sum_weff = 0.0, num_zero_sig = 0.0, sum_eff = 0.0;
        for( size_t k = i; k < j; ++k )
        {
          const double fs = std::sqrt( pts[k].fracStatUncert*pts[k].fracStatUncert
                                       + pts[k].fracCertUncert*pts[k].fracCertUncert );
          const double sigma = fs * pts[k].efficiency;
          sum_eff += pts[k].efficiency;
          if( sigma > 0.0 )
          {
            sum_w += 1.0 / (sigma*sigma);
            sum_weff += pts[k].efficiency / (sigma*sigma);
          }else
          {
            num_zero_sig += 1.0;
          }
        }//for( size_t k = i; k < j; ++k )

        double eff, frac_sigma;
        if( num_zero_sig > 0.0 )
        {
          // Any zero-uncertainty point dominates an inverse-variance average.
          eff = sum_eff / double(j - i);
          frac_sigma = 0.0;
        }else
        {
          eff = sum_weff / sum_w;
          frac_sigma = 1.0 / ( std::sqrt(sum_w) * eff );
        }

        answer.curve.energies_keV.push_back( pts[i].energy );
        answer.curve.eff.push_back( eff );
        answer.curve.frac_sigma.push_back( frac_sigma );

        i = j;
      }//for( size_t i = 0; i < pts.size(); )

      if( answer.curve.energies_keV.size() >= 2 )
      {
        answer.ref_distance_cm = (sum_d / pts.size()) / PhysicalUnits::cm;
        answer.curve_derived = false;
        return answer;
      }

      answer.curve = ceelo::AnchorCurve{};
    }//if( >= 2 points, all at one distance )
  }//if( raw && !raw->empty() )

  // Fallback: sample the fitted intrinsic curve and convert it to absolute
  //  efficiency at a single reference distance (mirrors
  //  MakeMcResponseForDrf::groundingPointsForDrf, the vetted recipe -
  //  intrinsicEfficiency() already backs any air attenuation out, so the
  //  reconstructed absolute efficiencies are in-vacuum, consistent with the
  //  transfer kernel).
  answer.curve_derived = true;

  const double a_cm = geom.transverse_half_extent();
  double d_ref_cm = std::max( 50.0, 10.0 * a_cm );
  if( override_ref_distance_cm > 0.0 )
    d_ref_cm = override_ref_distance_cm;
  else if( (drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute)
           && (drf->absoluteEfficiencyDistance() > 0.0) )
    d_ref_cm = drf->absoluteEfficiencyDistance() / PhysicalUnits::cm;

  double e_lo = drf->lowerEnergy(), e_hi = drf->upperEnergy();
  if( (e_lo <= 0.0) || (e_hi <= e_lo) )
  {
    e_lo = 59.0;
    e_hi = 2614.0;
  }
  e_lo = std::max( e_lo, 20.0 );

  // 16 log-spaced samples, plus flanking samples just either side of each
  //  crystal K-edge: eta(E) = eff/K jumps at an edge even for a smooth
  //  curve (K does), and the transfer segments - but does not add nodes at -
  //  the edges, so un-flanked edges would leave 0/1-node segments.
  const int n_samples = 16;
  vector<double> energies;
  for( int i = 0; i < n_samples; ++i )
    energies.push_back( e_lo * std::pow( e_hi/e_lo, double(i)/(n_samples-1) ) );

  for( const double edge : geom.crystal_k_edges( e_lo, e_hi ) )
  {
    const double below = edge * (1.0 - 1.0e-3), above = edge * (1.0 + 1.0e-3);
    if( (below > e_lo) && (above < e_hi) )
    {
      energies.push_back( below );
      energies.push_back( above );
    }
  }

  sort( begin(energies), end(energies) );
  energies.erase( unique( begin(energies), end(energies) ), end(energies) );

  const vector<double> cov = drf->efficiencyFracCovariance( energies );
  const size_t ne = energies.size();

  const double diam = drf->detectorDiameter();
  const double dist = d_ref_cm * PhysicalUnits::cm;
  const double frac_solid_angle = DetectorPeakResponse::fractionalSolidAngle(
                                      diam, dist + drf->detectorSetback() );

  for( size_t i = 0; i < ne; ++i )
  {
    const double intrinsic = drf->intrinsicEfficiency( static_cast<float>(energies[i]) );
    if( intrinsic <= 0.0 )
      continue;

    double frac_sigma = 0.05;
    if( cov.size() == ne*ne )
      frac_sigma = std::max( 0.01, std::sqrt( std::max( 0.0, cov[i*ne + i] ) ) );

    answer.curve.energies_keV.push_back( energies[i] );
    answer.curve.eff.push_back( intrinsic * frac_solid_angle );
    answer.curve.frac_sigma.push_back( frac_sigma );
  }//for( size_t i = 0; i < ne; ++i )

  if( answer.curve.energies_keV.size() < 2 )
    throw runtime_error( "transferAnchorForDrf: fewer than two usable"
                         " efficiency points could be constructed." );

  answer.ref_distance_cm = d_ref_cm;

  return answer;
}//transferAnchorForDrf(...)


ceelo::AnchorCurve totalTransferAnchorForDrf(
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const TransferAnchor &anchor )
{
  ceelo::AnchorCurve tot_curve;

  if( !drf || !drf->isValid() || !drf->totalEfficiencyCurve() )
    return tot_curve;

  const double diam = drf->detectorDiameter();
  const double dist = anchor.ref_distance_cm * PhysicalUnits::cm;
  const double frac_solid_angle = DetectorPeakResponse::fractionalSolidAngle(
                                      diam, dist + drf->detectorSetback() );
  for( const double energy : anchor.curve.energies_keV )
  {
    const float tot = drf->totalIntrinsicEfficiency( static_cast<float>(energy) );
    if( tot > 0.0f )
    {
      tot_curve.energies_keV.push_back( energy );
      tot_curve.eff.push_back( tot * frac_solid_angle );
    }
  }

  return tot_curve;
}//totalTransferAnchorForDrf(...)


std::shared_ptr<ceelo::DetectorResponse> makeTransferResponse(
                      const ceelo::GeometryDescriptor &geom,
                      const TransferAnchor &anchor,
                      const ceelo::AnchorCurve &tot_curve,
                      const std::string &detector_name )
{
  if( anchor.curve.energies_keV.size() < 2 )
    throw runtime_error( "makeTransferResponse: need at least two anchor points." );

  // The kernel positions are in the crystal-face frame (z = 0 at the crystal
  //  face, source in front at negative z); the anchor distance is in the
  //  descriptor's reference-point convention - convert exactly the way
  //  DetectorResponse::query_position() does for evaluation-time queries.
  double z_cm = anchor.ref_distance_cm;
  if( geom.reference_point == ceelo::ReferencePoint::EndcapFront )
    z_cm += geom.endcap_front_offset_cm();
  const Eigen::Vector3d ref_pos( 0.0, 0.0, -z_cm );

  ceelo::TransferResponseOptions opts;
  opts.detector_name = detector_name;

  const ceelo::AnchorCurve * const tot_anchor
                       = (tot_curve.energies_keV.size() >= 2) ? &tot_curve : nullptr;

  return ceelo::make_transfer_response( geom, anchor.curve, ref_pos, tot_anchor, opts );
}//makeTransferResponse(...)


ceelo::GeometryDescriptor buildAngleGeometry( const AngleOutxContents &contents,
                                              std::vector<std::string> &warnings )
{
  if( !contents.hasGeometry )
    throw runtime_error( "buildAngleGeometry: ANGLE file has no detector geometry." );

  const double radius_cm = contents.crystalRadius / PhysicalUnits::cm;
  const double length_cm = contents.crystalLength / PhysicalUnits::cm;
  if( (radius_cm <= 0.0) || (length_cm <= 0.0) )
    throw runtime_error( "buildAngleGeometry: non-positive crystal dimensions." );

  const shared_ptr<const MaterialDB> matDb = MaterialDB::initialized()
                                              ? MaterialDB::instance() : nullptr;
  const SandiaDecay::SandiaDecayDataBase * const nucDb = DecayDataBaseServer::database();

  // Resolve an ANGLE material name to a CeeLo material spec: try MaterialDB by
  //  name, then a chemical-formula fallback (covers ANGLE trade names, density
  //  inline).  Notes a warning and returns false when unresolved.
  auto resolve = [&]( const string &name, ceelo::MaterialSpec &out ) -> bool {
    if( name.empty() || !matDb )
      return false;

    shared_ptr<const Material> mat;
    try { mat = matDb->material( name ); }catch( std::exception & ){}

    if( !mat )
    {
      // Chemical-formula fallbacks (each element needs an explicit count for
      //  Material::parseChemicalFormula; density inline as "d=").
      static const std::map<string,string> fallbacks = {
        { "Aluminium",    "Al1 d=2.699" },   // British spelling ANGLE uses
        { "Mylar/PET",    "C10H8O4 d=1.38" },
        { "Mylar",        "C10H8O4 d=1.38" },
        { "PET",          "C10H8O4 d=1.38" },
        { "Brass",        "Cu0.63Zn0.37 d=8.5" },
        { "Carbon fiber", "C1 d=1.8" },
        { "Carbon Fiber", "C1 d=1.8" },
      };
      const auto it = fallbacks.find( name );
      const string formula = (it != end(fallbacks)) ? it->second : name;
      try { mat = MaterialDB::materialFromChemicalFormula( formula, nucDb ); }
      catch( std::exception & ){}
    }//if( !mat )

    if( !mat )
    {
      warnings.push_back( "Could not resolve material '" + name + "'; skipped." );
      return false;
    }

    try { out = to_ceelo_material( *mat ); }
    catch( std::exception &e )
    {
      warnings.push_back( "Material '" + name + "': " + string(e.what()) );
      return false;
    }
    return true;
  };//resolve lambda

  ceelo::GeometryDescriptor gd;
  gd.shape = ceelo::DetectorShape::Cylinder;
  gd.dimensions_cm = { radius_cm, length_cm };
  gd.symmetry = ceelo::ResponseSymmetry::Axial;
  gd.reference_point = ceelo::ReferencePoint::EndcapFront;

  // Crystal: ANGLE HPGe detectors are germanium.
  ceelo::MaterialSpec ge;
  if( resolve( "Ge", ge ) )
  {
    gd.crystal_material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( ge );
  }else
  {
    gd.crystal_material_index = -1;
    warnings.push_back( "Could not resolve crystal germanium material." );
  }

  // Coaxial bore from the core (bulletizing / rounding not represented by CeeLo).
  if( contents.hasCore && (contents.coreRadius > 0.0) && (contents.coreDepth > 0.0) )
    gd.bore = ceelo::BoreHoleConfig{ contents.coreRadius / PhysicalUnits::cm,
                                     contents.coreDepth / PhysicalUnits::cm };

  // Dead layer (inactive Ge); the thin germanium contact reads as dead Ge, so
  //  fold it into the side dead layer.
  {
    const double front = contents.deadLayerFront / PhysicalUnits::cm;
    double side = contents.deadLayerSide / PhysicalUnits::cm;
    if( contents.contactSide > 0.0 )
      side += contents.contactSide / PhysicalUnits::cm;
    if( (front > 0.0) || (side > 0.0) )
      gd.dead_layer = ceelo::DeadLayerConfig{ front, side, 0.0 };
  }

  // Endcap / window / housing layers, innermost -> outermost.  The vacuum gap
  //  is a spacing (attenuation ~ 0), not a layer, so it is dropped; the user
  //  reviews / corrects the seeded geometry before generating a response.
  for( const AngleOutxContents::Layer &layer : contents.layers )
  {
    const double front = layer.front / PhysicalUnits::cm;
    const double side = layer.side / PhysicalUnits::cm;
    if( (front <= 0.0) && (side <= 0.0) )
      continue;

    ceelo::MaterialSpec mspec;
    if( !resolve( layer.material, mspec ) )
      continue;

    ceelo::LayerSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( mspec );
    spec.front_thickness_cm = front;
    spec.side_thickness_cm = side;
    spec.z_start_cm = 0.0;
    spec.z_end_cm = length_cm;
    gd.layers.push_back( spec );
  }//for( const AngleOutxContents::Layer &layer : contents.layers )

  return gd;
}//buildAngleGeometry(...)


std::shared_ptr<DetectorPeakResponse> buildAngleSeedDrf(
                      const AngleOutxContents &contents )
{
  if( contents.referencePoints.size() < 2 )
    throw runtime_error( "buildAngleSeedDrf: need at least two reference points." );

  const double diameter = 2.0 * contents.crystalRadius;      // PhysicalUnits
  const double dist = contents.referenceDistanceCm * PhysicalUnits::cm;
  if( (diameter <= 0.0) || (dist <= 0.0) )
    throw runtime_error( "buildAngleSeedDrf: non-positive crystal diameter or"
                         " reference distance." );

  vector<DetectorPeakResponse::EnergyEffPoint> effpts;
  vector<MeasuredEffPoint> measpts;
  for( const pair<float,float> &ep : contents.referencePoints )
  {
    if( (ep.first <= 0.0f) || (ep.second <= 0.0f) )
      continue;

    DetectorPeakResponse::EnergyEffPoint e;
    e.energy = ep.first;
    e.efficiency = ep.second;
    effpts.push_back( e );

    MeasuredEffPoint m;
    m.energy = ep.first;
    m.efficiency = ep.second;
    m.distance = static_cast<float>( dist );
    measpts.push_back( m );
  }//for( const pair<float,float> &ep : contents.referencePoints )

  if( effpts.size() < 2 )
    throw runtime_error( "buildAngleSeedDrf: fewer than two positive reference"
                         " points." );

  auto drf = make_shared<DetectorPeakResponse>();
  drf->setEfficiencyPoints( effpts, static_cast<float>(diameter), dist,
                        DetectorPeakResponse::EffGeometryType::FarFieldAbsolute );
  drf->setName( contents.detName.empty() ? string("ANGLE detector") : contents.detName );
  if( !contents.detDescription.empty() )
    drf->setDescription( contents.detDescription );

  auto meas = make_shared<MeasuredDrfPoints>();
  meas->setPoints( measpts );
  drf->setMeasuredPoints( meas );

  return drf;
}//buildAngleSeedDrf(...)

}//namespace CeeLoUtils
