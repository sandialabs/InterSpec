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
#include <array>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>
#include <cassert>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <Eigen/Core>

#include "io/SolidAngle.h"

#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/AngleOutxImport.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/GadrasDetectorDat.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/MassAttenuationTool.h"

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


ceelo::MaterialSpec toCeeloMaterial( const AngleMaterial &mat,
                                     const SandiaDecay::SandiaDecayDataBase *db )
{
  if( !db )
    throw runtime_error( "toCeeloMaterial: no nuclide database." );

  if( mat.density_g_cm3 <= 0.0 )
    throw runtime_error( "Material '" + mat.name + "' has no density." );

  // Accumulate mass fractions by Z.  <elements> gives them directly; a
  //  <compound> gives atom counts, which become masses through the natural
  //  atomic mass; <compounds> are those, merged by their own mass fractions.
  map<int,double> frac_by_z;

  auto add = [&]( const string &symbol, const double mass ) {
    const SandiaDecay::Element * const el = db->element( symbol );
    if( !el )
      throw runtime_error( "Material '" + mat.name + "': '" + symbol
                           + "' is not an element." );
    frac_by_z[el->atomicNumber] += mass;
  };//add lambda

  for( const AngleElementFraction &ef : mat.elements )
    add( ef.symbol, ef.massFraction );

  for( const AngleCompound &cmp : mat.compounds )
  {
    map<int,double> cmp_mass;
    double cmp_total = 0.0;
    for( const AngleCompoundElement &ce : cmp.elements )
    {
      const SandiaDecay::Element * const el = db->element( ce.symbol );
      if( !el )
        throw runtime_error( "Material '" + mat.name + "': '" + ce.symbol
                             + "' is not an element." );
      const double mass = ce.atoms * el->atomicMass();
      cmp_mass[el->atomicNumber] += mass;
      cmp_total += mass;
    }//for( const AngleCompoundElement &ce : cmp.elements )

    if( cmp_total <= 0.0 )
      continue;

    for( const auto &zm : cmp_mass )
      frac_by_z[zm.first] += cmp.massFraction * (zm.second / cmp_total);
  }//for( const AngleCompound &cmp : mat.compounds )

  double sum = 0.0;
  for( const auto &zf : frac_by_z )
    sum += zf.second;
  if( sum <= 0.0 )
    throw runtime_error( "Material '" + mat.name + "' has no composition." );

  ceelo::MaterialSpec spec;
  spec.name = mat.name;
  spec.density_g_per_cm3 = mat.density_g_cm3;
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
}//toCeeloMaterial(...)


double pinnedEfficiencyDistanceCm( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                   const ceelo::GeometryDescriptor &geom )
{
  if( !drf || !drf->isValid() )
    return -1.0;

  // Raw measured points all taken at one distance: that IS the characterization geometry.
  const std::shared_ptr<const MeasuredDrfPoints> raw = drf->measuredPoints();
  if( raw )
  {
    double min_d = std::numeric_limits<double>::max(), max_d = 0.0, sum_d = 0.0;
    size_t n = 0;
    for( const MeasuredEffPoint &p : raw->points() )
    {
      if( (p.distance > 0.0f) && (p.efficiency > 0.0f) && (p.energy > 0.0f) )
      {
        min_d = std::min( min_d, double(p.distance) );
        max_d = std::max( max_d, double(p.distance) );
        sum_d += p.distance;
        n += 1;
      }
    }

    if( (n >= 2) && (max_d < 1.01*min_d) )
      return (sum_d / n) / PhysicalUnits::cm;
  }//if( raw )

  if( (drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute)
      && (drf->absoluteEfficiencyDistance() > 0.0) )
    return drf->absoluteEfficiencyDistance() / PhysicalUnits::cm;

  const double a_cm = geom.transverse_half_extent();
  return std::max( 50.0, 10.0 * a_cm );
}//pinnedEfficiencyDistanceCm(...)


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


bool attachCurveTransferResponse( DetectorPeakResponse &drf )
{
  if( drf.ceeloResponse() )
    return false;   //already has support; a generated response is better than this one

  const shared_ptr<const ceelo::GeometryDescriptor> geom = drf.geometry();
  if( !geom || !drf.isValid() || drf.isFixedGeometry() )
    return false;

  try
  {
    // `transferAnchorForDrf` takes a shared_ptr; alias one that does not own `drf` - it only reads
    //  it, and does not outlive this call.
    const shared_ptr<const DetectorPeakResponse> drf_ptr( shared_ptr<const void>(), &drf );

    const TransferAnchor anchor = transferAnchorForDrf( drf_ptr, *geom, -1.0 );
    const ceelo::AnchorCurve tot_curve = totalTransferAnchorForDrf( drf_ptr, anchor );
    const shared_ptr<ceelo::DetectorResponse> resp
                    = makeTransferResponse( *geom, anchor, tot_curve, drf.name() );
    if( !resp )
      return false;

    drf.setCeeloResponse( resp );
    return true;
  }catch( std::exception & )
  {
    // No usable anchor (too few points, a geometry the ray-tracer rejects, ...).  The detector is
    //  still perfectly good; it just keeps the flat-disk treatment it would have had anyway.
    return false;
  }
}//attachCurveTransferResponse(...)


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


const char *geometryProblemMsgId( const ceelo::GeometryProblem problem )
{
  switch( problem )
  {
    case ceelo::GeometryProblem::DimensionsMissing:     return "dgi-err-dims";
    case ceelo::GeometryProblem::DeadLayerTooThick:    return "dgi-err-dead-thick";
    case ceelo::GeometryProblem::BulletOnNonCylinder:  return "dgi-err-bullet-shape";
    case ceelo::GeometryProblem::BulletNotFinite:       return "dgi-err-bullet-invalid";
    case ceelo::GeometryProblem::BulletTooWide:         return "dgi-err-bullet-wide";
    case ceelo::GeometryProblem::BulletTooLong:         return "dgi-err-bullet-long";
    case ceelo::GeometryProblem::BulletNoDeadLayerRoom: return "dgi-err-bullet-dead";
    case ceelo::GeometryProblem::BoreOnNonCylinder:     return "dgi-err-bore-shape";
    case ceelo::GeometryProblem::BoreNotFinite:         return "dgi-err-bore-invalid";
    case ceelo::GeometryProblem::BoreTooWide:           return "dgi-err-bore-size";
    case ceelo::GeometryProblem::BoreTooDeep:           return "dgi-err-bore-size";
    case ceelo::GeometryProblem::BoreTipTooBlunt:       return "dgi-err-bore-round";
    case ceelo::GeometryProblem::BoreOutsideFillet:     return "dgi-err-bore-fillet";
    case ceelo::GeometryProblem::BoreInsideDeadLayer:   return "dgi-err-bore-dead";
  }//switch( problem )

  assert( 0 );  //a GeometryProblem was added without a message for it
  return "dgi-err-bore-size";
}//geometryProblemMsgId(...)


void relaxGeometryFeatures( ceelo::GeometryDescriptor &gd,
                            std::vector<std::string> &warnings )
{
  const std::vector<ceelo::GeometryProblem> problems = gd.problems();

  bool drop_fillet = false, drop_round_tip = false;
  for( const ceelo::GeometryProblem problem : problems )
  {
    switch( problem )
    {
      case ceelo::GeometryProblem::BulletOnNonCylinder:
      case ceelo::GeometryProblem::BulletNotFinite:
      case ceelo::GeometryProblem::BulletTooWide:
      case ceelo::GeometryProblem::BulletTooLong:
      case ceelo::GeometryProblem::BulletNoDeadLayerRoom:
      case ceelo::GeometryProblem::BoreOutsideFillet:
        drop_fillet = true;
        break;

      case ceelo::GeometryProblem::BoreTipTooBlunt:
        drop_round_tip = true;
        break;

      // Everything else is about the crystal, bore or dead layer themselves,
      //  which we won't silently rewrite - left for the caller to report.
      case ceelo::GeometryProblem::DimensionsMissing:
      case ceelo::GeometryProblem::DeadLayerTooThick:
      case ceelo::GeometryProblem::BoreOnNonCylinder:
      case ceelo::GeometryProblem::BoreNotFinite:
      case ceelo::GeometryProblem::BoreTooWide:
      case ceelo::GeometryProblem::BoreTooDeep:
      case ceelo::GeometryProblem::BoreInsideDeadLayer:
        break;
    }//switch( problem )
  }//for( const ceelo::GeometryProblem problem : problems )

  if( drop_fillet && (gd.bullet_radius_cm != 0.0) )
  {
    char msg[256];
    snprintf( msg, sizeof(msg), "The crystal's %.3g cm front-edge bulletizing radius"
              " doesn't fit this crystal, so a sharp front edge was used instead.",
              gd.bullet_radius_cm );
    warnings.push_back( msg );
    gd.bullet_radius_cm = 0.0;
  }//if( drop_fillet && ... )

  if( drop_round_tip && gd.bore && gd.bore->rounded_tip )
  {
    warnings.push_back( "The bore hole is shallower than it is wide, so its rounded"
                        " tip couldn't be modeled; a flat-bottomed bore was used." );
    gd.bore->rounded_tip = false;
  }//if( drop_round_tip && ... )

#if( PERFORM_DEVELOPER_CHECKS )
  // Dropping the fillet can only remove problems, never add one, so one
  //  re-check confirms we're done; anything left is about the crystal, bore or
  //  dead layer, which we deliberately don't rewrite.
  {
    const std::vector<ceelo::GeometryProblem> remaining = gd.problems();
    for( const ceelo::GeometryProblem problem : remaining )
    {
      const bool is_feature = (problem == ceelo::GeometryProblem::BoreTipTooBlunt)
                              || (problem == ceelo::GeometryProblem::BoreOutsideFillet)
                              || (problem == ceelo::GeometryProblem::BulletTooWide)
                              || (problem == ceelo::GeometryProblem::BulletTooLong)
                              || (problem == ceelo::GeometryProblem::BulletNotFinite)
                              || (problem == ceelo::GeometryProblem::BulletOnNonCylinder)
                              || (problem == ceelo::GeometryProblem::BulletNoDeadLayerRoom);
      if( is_feature )
        log_developer_error( __func__, ("relaxGeometryFeatures left a droppable problem: "
                                        + std::string(ceelo::to_string(problem))).c_str() );
    }//for( const ceelo::GeometryProblem problem : remaining )
  }
#endif //PERFORM_DEVELOPER_CHECKS
}//relaxGeometryFeatures(...)


ceelo::GeometryDescriptor buildAngleGeometry( const AngleOutxContents &contents,
                                              std::vector<std::string> &warnings )
{
  if( !contents.hasGeometry )
    throw runtime_error( "buildAngleGeometry: ANGLE file has no detector geometry." );

  // Some ANGLE detector types have no CeeLo shape at all (a well cavity, or a
  //  true-coaxial through-bore).  Approximating them as a solid or blind-bored
  //  cylinder is wrong exactly on-axis, where a transfer is anchored, so refuse
  //  and carry the reason - the fixed-geometry curve import still works.
  if( !contents.modeASupported )
    throw runtime_error( contents.modeAObstruction.empty()
                         ? string("buildAngleGeometry: this ANGLE detector cannot be"
                                  " represented as a Monte-Carlo geometry.")
                         : contents.modeAObstruction );

  const double radius_cm = contents.crystalRadius / PhysicalUnits::cm;
  const double length_cm = contents.crystalLength / PhysicalUnits::cm;
  if( (radius_cm <= 0.0) || (length_cm <= 0.0) )
    throw runtime_error( "buildAngleGeometry: non-positive crystal dimensions." );

  const shared_ptr<const MaterialDB> matDb = MaterialDB::initialized()
                                              ? MaterialDB::instance() : nullptr;
  const SandiaDecay::SandiaDecayDataBase * const nucDb = DecayDataBaseServer::database();

  // Resolve an ANGLE material to a CeeLo material spec, in order of authority:
  //
  //   1. an inline definition in the file (density + composition) is used as
  //      given - it needs no MaterialDB at all, and is what a user-defined
  //      ANGLE material always carries;
  //   2. otherwise the name is one of ANGLE's 19 predefined materials, or a
  //      trade name, and is mapped onto a MaterialDB entry;
  //   3. failing that, the name is tried as a chemical formula.
  //
  // The resulting MaterialSpec::name must be something the geometry form can
  //  re-resolve, so step 2/3 keep the resolved MaterialDB name.
  auto resolve_named = [&]( const string &angle_name, ceelo::MaterialSpec &out ) -> bool {
    if( angle_name.empty() || !matDb )
      return false;

    // ANGLE spelling -> MaterialDB name (see data/MaterialDataBase.txt).  ANGLE
    //  ships 19 predefined materials; the ones not listed here already resolve
    //  under their own name (Air, Beryllium, Brass, Copper, Magnesium, Nickel,
    //  Water, and the entries added for this import).
    // `MaterialDB::material()` already matches case-insensitively on the name, on
    //  a "Sym (name)" entry's symbol or parenthesised name, and on the
    //  description - so ANGLE's "Germanium", "Gold", "Air", "Beryllium",
    //  "Copper", "Magnesium", "Nickel" and "Water" all resolve unaided.  Only
    //  the spellings that differ by more than case need mapping.
    static const std::map<string,string> name_map = {
      { "aluminium",         "Aluminum" },      // the British spelling ANGLE uses
      { "aluminium oxide",   "Aluminum oxide" },
      { "carbon fiber",      "CFRP" },
      { "mylar/pet",         "Mylar" },
      { "pet",               "Mylar" },
      { "nai",               "sodium_iodide" },
      { "calcium carbonate", "calcium_carbonate" },
      { "magnesium oxide",   "magnesium_oxide" },
      // ANGLE's two generic trade names.  It does not publish their
      //  compositions, so stand in a common representative and say so; the
      //  user can correct the material in the geometry form.
      { "glass",             "pyrex_glass" },
      { "plastic",           "Polyethylene" },
    };

    const string key = SpecUtils::to_lower_ascii_copy( SpecUtils::trim_copy(angle_name) );
    const auto nit = name_map.find( key );
    const string name = (nit != end(name_map)) ? nit->second : angle_name;

    shared_ptr<const Material> mat;
    try { mat = matDb->material( name ); }catch( std::exception & ){}

    if( !mat )
    {
      // Last-resort chemical-formula fallback (each element an explicit count;
      //  density inline as "d=").
      try { mat = MaterialDB::materialFromChemicalFormula( name, nucDb ); }
      catch( std::exception & ){}
    }//if( !mat )

    if( !mat )
      return false;

    if( (key == "glass") || (key == "plastic") )
      warnings.push_back( "ANGLE does not publish a composition for its '" + angle_name
                          + "' material; '" + name + "' was used instead.  Change the"
                          " material if that is not what the sample was made of." );

    try { out = to_ceelo_material( *mat ); }
    catch( std::exception &e )
    {
      warnings.push_back( "Material '" + angle_name + "': " + string(e.what()) );
      return false;
    }
    return true;
  };//resolve_named lambda

  auto resolve = [&]( const AngleMaterial &material, ceelo::MaterialSpec &out ) -> bool {
    if( !material.isBareReference() )
    {
      try
      {
        out = toCeeloMaterial( material, nucDb );
        return true;
      }catch( std::exception &e )
      {
        warnings.push_back( "Material '" + material.name + "' as defined in the file could"
                            " not be used (" + string(e.what()) + "); falling back to its"
                            " name." );
      }
    }//if( !material.isBareReference() )

    return resolve_named( material.name, out );
  };//resolve lambda

  // A near-transparent spacer (rho ~ 1e-25, single H): preserves a layer's
  //  physical extent (and thus the crystal recess) while contributing ~zero
  //  attenuation.  MaterialDB-independent, so it works even with no DB.
  auto make_transparent = []() -> ceelo::MaterialSpec {
    ceelo::MaterialSpec s;
    s.name = "unresolved-spacer";
    s.density_g_per_cm3 = 1.0e-25;
    s.composition = { ceelo::MaterialComponent{ 1, 1.0 } };
    return s;
  };//make_transparent lambda

  ceelo::GeometryDescriptor gd;
  gd.shape = ceelo::DetectorShape::Cylinder;
  gd.dimensions_cm = { radius_cm, length_cm };
  gd.symmetry = ceelo::ResponseSymmetry::Axial;
  gd.reference_point = ceelo::ReferencePoint::EndcapFront;

  // Crystal material from the ANGLE detector type, using the same CeeLo built-ins
  //  the geometry form offers (so its Crystal combo matches the seeded name by
  //  `.name`).  ANGLE's eight types are all germanium or NaI; the substring
  //  match is kept only for a `type` string this version does not know, so a
  //  hand-edited CZT / LaBr3 file still lands on the right crystal.
  ceelo::Material crystalMat = ceelo::make_HPGe();
  switch( contents.detTypeEnum )
  {
    case AngleDetectorType::NaI:
    case AngleDetectorType::NaIWell:
      crystalMat = ceelo::make_NaI();
      break;

    case AngleDetectorType::ClosedEndCoaxHPGe:
    case AngleDetectorType::TrueCoaxHPGe:
    case AngleDetectorType::ClosedEndCoaxGeLi:
    case AngleDetectorType::OpenEndCoaxGeLi:
    case AngleDetectorType::PlanarLEPD:
    case AngleDetectorType::Well:
      break;  //HPGe

    case AngleDetectorType::Unknown:
      if( SpecUtils::icontains( contents.detType, "NaI" ) )
        crystalMat = ceelo::make_NaI();
      else if( SpecUtils::icontains( contents.detType, "CZT" )
               || SpecUtils::icontains( contents.detType, "CdZnTe" ) )
        crystalMat = ceelo::make_CZT();
      else if( SpecUtils::icontains( contents.detType, "LaBr" ) )
        crystalMat = ceelo::make_LaBr3();
      break;
  }//switch( contents.detTypeEnum )

  gd.crystal_material_index = static_cast<int>( gd.materials.size() );
  gd.materials.push_back( ceelo::MaterialSpec::from( crystalMat ) );

  // Front-edge bulletization: a quarter-torus fillet rounding off the outer
  //  front edge of the crystal.  Zero is a sharp 90-degree edge.
  gd.bullet_radius_cm = contents.bulletizingRadius / PhysicalUnits::cm;

  // Coaxial bore from the core; `rounded` means a round-tipped drill, so the
  //  closed end is a hemisphere of the bore radius (depth is to the apex).
  if( contents.hasCore && (contents.coreRadius > 0.0) && (contents.coreDepth > 0.0) )
  {
    gd.bore = ceelo::BoreHoleConfig{ contents.coreRadius / PhysicalUnits::cm,
                                     contents.coreDepth / PhysicalUnits::cm,
                                     contents.coreRounded };

    // The contact pin runs down the bore; CeeLo's bore is empty.  Usually a thin
    //  rod whose effect is negligible, but say so when it fills much of the bore.
    if( contents.contactPin.radius > (0.5 * contents.coreRadius) )
      warnings.push_back( "The detector's contact pin fills much of the crystal's bore"
                          " (radius " + SpecUtils::printCompact( contents.contactPin.radius
                                                        / PhysicalUnits::mm, 3 )
                          + " mm of " + SpecUtils::printCompact( contents.coreRadius
                                                        / PhysicalUnits::mm, 3 )
                          + " mm); the bore is modeled as empty, so its attenuation is"
                          " not accounted for." );
  }//if( a bore )

  // Dead layer (inactive Ge).  A thin implanted / Li-diffused contact reads as
  //  dead crystal, so fold it in; a thick contact of some other material is one
  //  of `layers` instead.  The back dead layer is a planar-detector thing.
  {
    double front = contents.deadLayerFront / PhysicalUnits::cm;
    double side = contents.deadLayerSide / PhysicalUnits::cm;
    const double back = contents.deadLayerBack / PhysicalUnits::cm;

    if( contents.contact.readsAsDeadCrystal )
    {
      front += contents.contact.front / PhysicalUnits::cm;
      side += contents.contact.side / PhysicalUnits::cm;
    }

    if( (front > 0.0) || (side > 0.0) || (back > 0.0) )
      gd.dead_layer = ceelo::DeadLayerConfig{ front, side, back };
  }

  // The concentric stack, innermost -> outermost, exactly as the file gives it:
  //  reflecting layer, housing, cryostat gap, antimicrophonic shield, endcap,
  //  window and coatings.  The user reviews / corrects the seeded geometry
  //  before generating a response.
  //
  //  A vacuum gap carries no material: CeeLo has no dedicated gap field, so it
  //  is represented by a near-transparent spacer ("galactic vacuum",
  //  rho ~ 1e-25).  That adds to endcap_front_offset_cm() - correctly recessing
  //  the crystal - while contributing ~zero attenuation.
  //
  //  Pieces behind the crystal (a NaI detector's optical coupling and
  //  photomultiplier tube) have nowhere to go: CeeLo attenuators cover the front
  //  and sides only.  They are dropped with a warning; they affect backscatter,
  //  not the full-energy peak.
  for( const AngleOutxContents::Layer &layer : contents.layers )
  {
    const double front = layer.front / PhysicalUnits::cm;
    const double side = layer.side / PhysicalUnits::cm;
    if( (front <= 0.0) && (side <= 0.0) )
      continue;

    if( layer.behindCrystal )
    {
      warnings.push_back( "The detector's " + string(to_string(layer.role)) + " sits behind the"
                          " crystal, which is not modeled; it affects backscatter, not the"
                          " full-energy peak." );
      continue;
    }//if( layer.behindCrystal )

    const bool is_vacuum = (layer.role == AngleLayerRole::VacuumInner)
                           || (layer.role == AngleLayerRole::VacuumOuter);

    ceelo::MaterialSpec mspec;
    if( is_vacuum )
    {
      if( !resolve_named( "galactic vacuum", mspec ) )
        mspec = make_transparent();
    }else if( !resolve( layer.material, mspec ) )
    {
      // For an unresolvable material still insert a layer so the physical extent
      //  (and hence the crystal recess / endcap-front offset) is preserved; use a
      //  near-transparent spacer (rho ~ 1e-25) so it contributes ~zero
      //  attenuation.  This is MaterialDB-independent, so it works even when the
      //  named material cannot be resolved.
      warnings.push_back( "Material '" + layer.material.name + "' (the "
                          + to_string(layer.role) + ") unresolved; inserting a transparent"
                          " spacer to preserve its physical extent (its attenuation is NOT"
                          " modeled - please set the material)." );
      mspec = make_transparent();
    }//if( is_vacuum ) / else if( !resolve(...) )

    ceelo::LayerSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( mspec );
    spec.front_thickness_cm = front;
    spec.side_thickness_cm = side;
    spec.z_start_cm = 0.0;
    spec.z_end_cm = length_cm;
    gd.layers.push_back( spec );
  }//for( const AngleOutxContents::Layer &layer : contents.layers )

  // The fillet and the rounded bore tip are refinements on top of a plain
  //  cylinder; if this file's values can't be represented, drop them and keep
  //  the import rather than failing it.  Done last, since the fillet check
  //  needs the dead layer.
  relaxGeometryFeatures( gd, warnings );

#if( PERFORM_DEVELOPER_CHECKS )
  // Self-consistency check on the assembled endcap-front offset: since every
  //  front-stack element (vacuum gap, each endcap/window layer) is now always
  //  represented (unresolved materials become transparent spacers rather than
  //  being dropped), the descriptor's endcap_front_offset_cm() must equal the
  //  sum of front thicknesses parsed from the file.  A mismatch means a layer
  //  was silently lost - the bug class that collapsed the crystal recess and
  //  produced a flat ~20-30% efficiency error.
  //
  //  The inactive-Ge thickness is NOT part of this sum: it is carved out of the
  //  inside of the crystal, so it moves the active face, not the crystal face.
  //  The cryostat vacuum gap IS part of it, and arrives as one of `layers`.
  {
    double expect_front = 0.0;
    for( const AngleOutxContents::Layer &layer : contents.layers )
    {
      const double f = layer.front / PhysicalUnits::cm;
      if( (f > 0.0) && !layer.behindCrystal )
        expect_front += f;
    }
    const double got = gd.endcap_front_offset_cm();
    if( std::fabs( got - expect_front ) > (1.0e-4 + 1.0e-3*expect_front) )
    {
      char msg[256];
      snprintf( msg, sizeof(msg), "endcap_front_offset (%.4f cm) != sum of parsed"
                " front thicknesses (%.4f cm) - a front-stack layer was dropped.",
                got, expect_front );
      log_developer_error( __func__, msg );
    }
  }
#endif //PERFORM_DEVELOPER_CHECKS

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


//=========================== GADRAS Detector.dat =============================

namespace
{
  /** Minimum geometric thickness for a synthesized generic attenuator, cm.

   Only reached when GADRAS gives no setback to divide up (4 of the 13 shipped
   generic detectors).  Kept small so it perturbs the solid angle negligibly;
   the cost is that the resulting shell is very dense, and therefore too opaque
   for rays crossing it at a grazing angle.
   */
  const double sm_min_generic_layer_cm = 1.0e-3;

  /** Parses one element symbol + optional atom count out of `text` at `pos`.
   Handles "Cs", "I2", "Cd0.96" and "[H2]" (an isotope, folded onto H).
   Returns false at the end of the string.
   */
  bool next_gadras_term( const string &text, size_t &pos,
                         string &symbol, double &atoms )
  {
    while( (pos < text.size()) && std::isspace( static_cast<unsigned char>(text[pos]) ) )
      pos += 1;

    if( pos >= text.size() )
      return false;

    symbol.clear();
    atoms = 1.0;

    if( text[pos] == '[' )
    {
      // An isotope, e.g. "[H2]" - the element symbol, then a mass number we
      //  drop (the cross sections are per-element).
      const size_t close = text.find( ']', pos );
      if( close == string::npos )
        throw runtime_error( "unmatched '[' in '" + text + "'" );

      const string inside = text.substr( pos + 1, close - pos - 1 );
      size_t i = 0;
      while( (i < inside.size()) && std::isalpha( static_cast<unsigned char>(inside[i]) ) )
        i += 1;
      symbol = inside.substr( 0, i );
      pos = close + 1;
    }else
    {
      if( !std::isupper( static_cast<unsigned char>(text[pos]) ) )
        throw runtime_error( "expected an element symbol at offset "
                             + std::to_string(pos) + " of '" + text + "'" );
      symbol += text[pos];
      pos += 1;
      while( (pos < text.size()) && std::islower( static_cast<unsigned char>(text[pos]) ) )
      {
        symbol += text[pos];
        pos += 1;
      }
    }//if( isotope bracket ) / else

    if( symbol.empty() )
      throw runtime_error( "empty element symbol in '" + text + "'" );

    // Optional atom count; may be fractional ("Cd0.96").  Absent means 1.
    const size_t num_start = pos;
    while( (pos < text.size())
           && (std::isdigit( static_cast<unsigned char>(text[pos]) ) || (text[pos] == '.')) )
      pos += 1;

    if( pos > num_start )
    {
      const string numstr = text.substr( num_start, pos - num_start );
      if( !(stringstream(numstr) >> atoms) || (atoms <= 0.0) )
        throw runtime_error( "'" + numstr + "' is not a valid atom count in '" + text + "'" );
    }

    return true;
  }//next_gadras_term(...)


  /** Packs a Z -> mass fraction map into a MaterialSpec, normalizing and
   rejecting elements the Monte Carlo has no cross sections for. */
  ceelo::MaterialSpec spec_from_mass_by_z( const map<int,double> &mass_by_z,
                                           const double density_g_per_cm3,
                                           const string &name )
  {
    double sum = 0.0;
    for( const auto &zm : mass_by_z )
      sum += zm.second;

    if( sum <= 0.0 )
      throw runtime_error( "Material '" + name + "' has no composition." );

    ceelo::MaterialSpec spec;
    spec.name = name;
    spec.density_g_per_cm3 = density_g_per_cm3;
    for( const auto &zm : mass_by_z )
    {
      if( (zm.first < 1) || (zm.first > 92) )
        throw runtime_error( "Material '" + name + "' contains an element (Z="
                             + std::to_string(zm.first) + ") the Monte Carlo has"
                             " no cross-section data for." );
      ceelo::MaterialComponent c;
      c.Z = static_cast<uint8_t>( zm.first );
      c.mass_fraction = zm.second / sum;
      spec.composition.push_back( c );
    }

    return spec;
  }//spec_from_mass_by_z(...)
}//namespace


ceelo::MaterialSpec materialFromGadrasFormula( const std::string &formula,
                                               const double density_g_per_cm3,
                                               const std::string &name )
{
  if( density_g_per_cm3 <= 0.0 )
    throw runtime_error( "materialFromGadrasFormula: '" + name
                         + "' has a non-positive density." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "materialFromGadrasFormula: no nuclide database." );

  map<int,double> mass_by_z;
  size_t pos = 0;
  string symbol;
  double atoms = 0.0;

  try
  {
    while( next_gadras_term( formula, pos, symbol, atoms ) )
    {
      const SandiaDecay::Element * const el = db->element( symbol );
      if( !el )
        throw runtime_error( "'" + symbol + "' is not an element" );
      mass_by_z[el->atomicNumber] += atoms * el->atomicMass();
    }
  }catch( std::exception &e )
  {
    throw runtime_error( "Could not parse GADRAS formula '" + formula + "' for '"
                         + name + "': " + string(e.what()) );
  }

  return spec_from_mass_by_z( mass_by_z, density_g_per_cm3, name );
}//materialFromGadrasFormula(...)


ceelo::MaterialSpec gadrasCrystalMaterial( const std::string &gadras_material_name )
{
  if( gadras_material_name.empty() )
    throw runtime_error( "gadrasCrystalMaterial: no crystal material given." );

  // The CeeLo built-ins spell their names exactly like the GADRAS table, which
  //  is what lets the geometry form's crystal combo match by `.name`.
  string name = gadras_material_name;
  SpecUtils::trim( name );

  if( SpecUtils::iequals_ascii( name, "NaI" ) )
    return ceelo::MaterialSpec::from( ceelo::make_NaI() );
  if( SpecUtils::iequals_ascii( name, "LaBr3" ) )
    return ceelo::MaterialSpec::from( ceelo::make_LaBr3() );
  if( SpecUtils::iequals_ascii( name, "CeBr3" ) )
    return ceelo::MaterialSpec::from( ceelo::make_CeBr3() );
  if( SpecUtils::iequals_ascii( name, "HPGe" ) )
    return ceelo::MaterialSpec::from( ceelo::make_HPGe() );
  if( SpecUtils::iequals_ascii( name, "CZT" ) )
    return ceelo::MaterialSpec::from( ceelo::make_CZT() );
  if( SpecUtils::istarts_with( name, "CLYC" ) )
    return ceelo::MaterialSpec::from( ceelo::make_CLYC() );

  const GadrasDetectorDat::MaterialInfo * const info
                                    = GadrasDetectorDat::materialByName( name );
  if( !info )
    throw runtime_error( "gadrasCrystalMaterial: '" + gadras_material_name
                         + "' is not a GADRAS detector material." );

  return materialFromGadrasFormula( info->formula, info->density, info->name );
}//gadrasCrystalMaterial(...)


ceelo::MaterialSpec genericAttenuatorMaterial( const double atomic_number,
                                               const double areal_density_g_cm2,
                                               const double thickness_cm )
{
  if( (areal_density_g_cm2 <= 0.0) || (thickness_cm <= 0.0) )
    throw runtime_error( "genericAttenuatorMaterial: areal density and thickness"
                         " must both be positive." );

  if( (atomic_number < 1.0) || (atomic_number > 92.0) )
    throw runtime_error( "genericAttenuatorMaterial: atomic number "
                         + std::to_string(atomic_number) + " is outside the"
                         " range the Monte Carlo has cross sections for." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "genericAttenuatorMaterial: no nuclide database." );

  const int z_lo = static_cast<int>( std::floor( atomic_number ) );
  const int z_hi = std::min( 92, z_lo + 1 );
  const double f_hi = atomic_number - z_lo;   //0 when the AN is an integer

  map<int,double> mass_by_z;
  mass_by_z[z_lo] += (1.0 - f_hi);
  if( f_hi > 0.0 )
    mass_by_z[z_hi] += f_hi;

  const string namebuf = genericAttenuatorName( atomic_number, areal_density_g_cm2 );

  // The density is whatever makes this thickness carry the file's areal
  //  density; only their product enters the attenuation.
  const double density = areal_density_g_cm2 / thickness_cm;

  const ceelo::MaterialSpec spec = spec_from_mass_by_z( mass_by_z, density, namebuf );

#if( PERFORM_DEVELOPER_CHECKS )
  // The claim this rests on: a two-element mass mix attenuates like InterSpec's
  //  own fractional-atomic-number model, because a mixture's mu/rho is the
  //  MASS-weighted sum of its components'.  Checked directly against
  //  MassAttenuation - never by building a ceelo::Material, whose mu cache is
  //  keyed on (address, energy, density) and would happily return a previous
  //  material's numbers for a stack temporary at a reused address.
  if( f_hi > 0.0 )
  {
    for( const double energy : { 30.0, 60.0, 122.0, 662.0, 1332.0, 3000.0 } )
    {
      try
      {
        const double lo = MassAttenuation::massAttenuationCoefficientElement(
                                      z_lo, static_cast<float>(energy) );
        const double hi = MassAttenuation::massAttenuationCoefficientElement(
                                      z_hi, static_cast<float>(energy) );
        const double mix = (1.0 - f_hi)*lo + f_hi*hi;
        const double ref = MassAttenuation::massAttenuationCoefficientFracAN(
                                      static_cast<float>(atomic_number),
                                      static_cast<float>(energy) );

        // Skip energies where an absorption edge falls BETWEEN the two
        //  bracketing elements: their mu/rho then differ by a large factor and
        //  no interpolation across atomic number is meaningful there (e.g. 30
        //  keV for AN 50.5, between the Sn and Sb K-edges).  A real mixture
        //  shows both edges, which is what the transport should see; the two
        //  models are simply answering different questions.
        const double ratio = (lo > 0.0 && hi > 0.0) ? std::max(lo/hi, hi/lo) : 1.0;
        if( ratio > 1.5 )
          continue;

        if( (ref > 0.0) && (std::fabs(mix - ref) > 0.05*ref) )
        {
          char msg[256];
          snprintf( msg, sizeof(msg), "Generic attenuator AN=%.3f at %.0f keV:"
                    " two-element mass mix gives mu/rho=%.5g, but"
                    " massAttenuationCoefficientFracAN gives %.5g (%.1f%% off).",
                    atomic_number, energy, mix, ref, 100.0*(mix-ref)/ref );
          log_developer_error( __func__, msg );
        }
      }catch( std::exception & )
      {
        //no cross-section data at this energy/Z - nothing to check against
      }
    }//for( const double energy : ... )
  }//if( f_hi > 0.0 )
#endif //PERFORM_DEVELOPER_CHECKS

  return spec;
}//genericAttenuatorMaterial(...)


ceelo::GeometryDescriptor buildGadrasGeometry( const GadrasDetectorDat &dat,
                                               std::vector<std::string> &warnings )
{
  string crystal_name = dat.materialName();
  if( crystal_name.empty() )
  {
    // Not fatal.  A dropped file is usually just named "Detector.dat", so there
    //  is no directory name to recover the material from, and refusing would
    //  leave the user with nothing - when the geometry form has a crystal
    //  dropdown they can simply set.  So pick a stand-in and say plainly that it
    //  is one.
    //
    //  The resolution the file DOES state at least separates germanium from the
    //  scintillators, which is the difference that matters most here.
    const float fwhm661 = dat.resFWHM661();
    crystal_name = ((fwhm661 > 0.0f) && (fwhm661 < 1.5f)) ? "HPGe" : "NaI";

    warnings.push_back( "Detector.dat does not say what the crystal is made of,"
                        " and there is no detector name to infer it from, so "
                        + crystal_name + " was assumed."
                        + ((fwhm661 > 0.0f)
                            ? (" (Its " + SpecUtils::printCompact(fwhm661,3)
                               + "% resolution at 661 keV is the only clue.)")
                            : string())
                        + "  Set the crystal material in the geometry form - it"
                        " drives the whole efficiency, and the crystal shape with it." );
  }//if( crystal_name.empty() )

  // GADRAS stores an *effective* box (a fitted shape its analytic chord-length
  //  model reproduces the measured efficiency with), not the true crystal; the
  //  material name is what tells inferShape() which solid that box stands for.
  const GadrasDetectorDat::InferredShape shape = dat.inferShape( crystal_name );

  ceelo::GeometryDescriptor gd;
  gd.reference_point = ceelo::ReferencePoint::EndcapFront;

  switch( shape.shape )
  {
    case GadrasDetectorDat::Shape::Cylinder:
    case GadrasDetectorDat::Shape::CoaxialCylinder:
      gd.shape = ceelo::DetectorShape::Cylinder;
      gd.symmetry = ceelo::ResponseSymmetry::Axial;
      // The radius is the EQUAL-AREA one, not inferShape()'s nominal width/2.
      //  GADRAS stores an effective box and reports efficiency per photon
      //  striking a face of area width*height; the circle that intercepts the
      //  same photons is the equal-area one (which is also the diameter
      //  InterSpec's DRF carries, and what fractionalSolidAngle is evaluated
      //  over).  Using width/2 instead shrinks the face by a factor of pi/4 and
      //  the computed intrinsic efficiency comes out ~21% low across the board -
      //  measured, before this was fixed, as MC/GADRAS ratios clustered at 0.78.
      gd.dimensions_cm = { 0.5*dat.equivalentCircularDiameterCm(), shape.dimB };
      break;

    case GadrasDetectorDat::Shape::Box:
    case GadrasDetectorDat::Shape::Rectangular:
      gd.shape = ceelo::DetectorShape::Box;
      gd.symmetry = ceelo::ResponseSymmetry::Quadrant;
      // A box needs no equal-area correction: its face already IS width x
      //  height.  InferredShape gives depth x width x height; CeeLo wants the
      //  two transverse half-widths, then the axial length.
      gd.dimensions_cm = { 0.5*shape.dimB, 0.5*shape.dimC, shape.dimA };
      break;

    case GadrasDetectorDat::Shape::Unknown:
      throw runtime_error( "Could not work out the crystal shape from Detector.dat." );
  }//switch( shape.shape )

  for( const double d : gd.dimensions_cm )
  {
    if( (d <= 0.0) || std::isinf(d) || std::isnan(d) )
      throw runtime_error( "Detector.dat gives a non-positive crystal dimension"
                           " (crystal length " + std::to_string(dat.length())
                           + " cm, width " + std::to_string(dat.width())
                           + " cm), so its geometry cannot be modeled." );
  }

  gd.crystal_material_index = static_cast<int>( gd.materials.size() );
  gd.materials.push_back( gadrasCrystalMaterial( crystal_name ) );

  // GADRAS gives ONE dead-layer thickness with no direction; apply it to the
  //  front and the sides.  It is carved out of the crystal, so it moves where
  //  the active volume starts - never where the crystal face is.
  const double dead_cm = 0.1 * dat.deadLayerMm();
  if( dead_cm > 0.0 )
  {
    gd.dead_layer = ceelo::DeadLayerConfig{ dead_cm, dead_cm, 0.0 };
    warnings.push_back( "Detector.dat gives a single "
                        + SpecUtils::printCompact(dat.deadLayerMm(), 3)
                        + " mm dead-layer thickness with no direction; it was"
                        " applied to both the front face and the sides of the"
                        " crystal." );
  }//if( dead_cm > 0.0 )

  // GADRAS has no bore geometry at all, so a coaxial HPGe comes through as a
  //  solid cylinder.  This is the single largest accuracy caveat of a
  //  Detector.dat-only import, so say so rather than quietly modeling the
  //  wrong solid.
  if( shape.shape == GadrasDetectorDat::Shape::CoaxialCylinder )
    warnings.push_back( "Detector.dat carries no bore-hole geometry, so this"
                        " coaxial HPGe was modeled as a solid cylinder.  If you"
                        " know the bore diameter and depth, enter them in the"
                        " Geometry form - they matter most above ~200 keV." );

  // --- the attenuator stack ------------------------------------------------
  // GADRAS specifies each attenuator only as an effective atomic number plus an
  //  areal density; there is no thickness in the file and CeeLo has no generic
  //  areal-density layer.  The one length the file DOES give for the front
  //  stack is the setback (the crystal's recess behind the reference point), so
  //  divide that between the layers in proportion to their areal densities: the
  //  areal densities stay exact, endcap_front_offset_cm() reproduces the
  //  setback, and nothing is invented beyond how a known total splits.
  struct GenericLayer { double an; double ad; };
  vector<GenericLayer> generic;

  const GadrasDetectorDat::Attenuator inner = dat.innerAttenuator();
  const GadrasDetectorDat::Attenuator outer = dat.outerAttenuator();
  for( const GadrasDetectorDat::Attenuator *att : { &inner, &outer } )
  {
    if( att->arealDensity <= 0.0f )       //every shipped generic detector has outer AD == 0
      continue;

    if( (att->atomicNumber < 1.0f) || (att->atomicNumber > 92.0f) )
    {
      warnings.push_back( "An attenuator with an effective atomic number of "
                          + SpecUtils::printCompact(att->atomicNumber, 3)
                          + " was dropped; only 1 through 92 can be modeled." );
      continue;
    }

    if( att->porosityPercent > 0.0f )
      warnings.push_back( "The attenuator porosity ("
                          + SpecUtils::printCompact(att->porosityPercent, 3)
                          + "%) is not modeled; only its areal density is." );

    generic.push_back( GenericLayer{ att->atomicNumber, att->arealDensity } );
  }//for( const GadrasDetectorDat::Attenuator *att : { &inner, &outer } )

  const double setback_cm = dat.setbackCm();

  double total_ad = 0.0;
  for( const GenericLayer &g : generic )
    total_ad += g.ad;

  double used_front = 0.0;
  for( const GenericLayer &g : generic )
  {
    const double share = (total_ad > 0.0) ? (g.ad / total_ad) : 0.0;
    const double t = std::max( sm_min_generic_layer_cm, setback_cm * share );

    ceelo::LayerSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( genericAttenuatorMaterial( g.an, g.ad, t ) );
    spec.front_thickness_cm = t;
    spec.side_thickness_cm = t;
    spec.z_start_cm = 0.0;
    spec.z_end_cm = gd.dimensions_cm.back();
    gd.layers.push_back( spec );

    used_front += t;
  }//for( const GenericLayer &g : generic )

  // --- the side shield -----------------------------------------------------
  // CeeLo attenuators and collimators are concentric, segmented only along the
  //  axis; there is no azimuthal coverage anywhere in Geometry.  So a shield
  //  that wraps the crystal completely can be modeled, and a partial one cannot.
  //  Never average the areal density by the coverage: transmission is
  //  exp(-mu*AD), so a coverage-weighted AD understates the shielded fraction
  //  and overstates the bare one.
  const GadrasDetectorDat::SideShield side = dat.sideShield();
  if( side.arealDensity > 0.0f )
  {
    const bool full = (side.covPlusX >= 99.0f) && (side.covMinusX >= 99.0f)
                      && (side.covPlusY >= 99.0f) && (side.covMinusY >= 99.0f);
    if( full && (side.atomicNumber >= 1.0f) && (side.atomicNumber <= 92.0f) )
    {
      // Unlike the front stack there is no length in the file to divide up - the
      //  setback says nothing about a side wall - so this is a nominal wall
      //  thickness carrying the file's areal density.  It changes only the outer
      //  radius, never the attenuation.
      const double t = 0.1;
      ceelo::LayerSpec spec;
      spec.material_index = static_cast<int>( gd.materials.size() );
      gd.materials.push_back( genericAttenuatorMaterial( side.atomicNumber,
                                                         side.arealDensity, t ) );
      spec.front_thickness_cm = 0.0;   //a side shield, not part of the front stack
      spec.side_thickness_cm = t;
      spec.z_start_cm = 0.0;
      spec.z_end_cm = gd.dimensions_cm.back();
      gd.layers.push_back( spec );
    }else
    {
      warnings.push_back( "The side shield covers only part of the detector"
                          " (+X " + SpecUtils::printCompact(side.covPlusX,3)
                          + "%, -X " + SpecUtils::printCompact(side.covMinusX,3)
                          + "%, +Y " + SpecUtils::printCompact(side.covPlusY,3)
                          + "%, -Y " + SpecUtils::printCompact(side.covMinusY,3)
                          + "%), which cannot be modeled - shielding layers wrap"
                          " the crystal completely - so it was left out." );
    }
  }//if( side.arealDensity > 0.0f )

  // --- the back shield -----------------------------------------------------
  const GadrasDetectorDat::BackShield back = dat.backShield();
  if( (back.arealDensity > 0.0f) && (back.coverage > 0.0f) )
    warnings.push_back( "The back shield sits behind the crystal, which is not"
                        " modeled; it affects backscatter, not the full-energy"
                        " peak." );

  // --- pin the crystal recess to the GADRAS setback ------------------------
  // Whatever the attenuators did not account for becomes a near-transparent
  //  spacer, so endcap_front_offset_cm() == setbackCm() and a distance measured
  //  from the endcap face means the same thing to us as it does to GADRAS.
  if( setback_cm > (used_front + 1.0e-6) )
  {
    ceelo::MaterialSpec spacer;
    spacer.name = "endcap-standoff";
    spacer.density_g_per_cm3 = 1.0e-25;
    spacer.composition = { ceelo::MaterialComponent{ 1, 1.0 } };

    ceelo::LayerSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( spacer );
    spec.front_thickness_cm = setback_cm - used_front;
    spec.side_thickness_cm = 0.0;
    spec.z_start_cm = 0.0;
    spec.z_end_cm = gd.dimensions_cm.back();
    gd.layers.push_back( spec );
  }else if( used_front > (setback_cm + 1.0e-6) )
  {
    // Only when the minimum-thickness floor overran a tiny/absent setback.
    warnings.push_back( "The detector's attenuators had to be given a small"
                        " thickness, which puts the crystal "
                        + SpecUtils::printCompact(used_front - setback_cm, 3)
                        + " cm further back than Detector.dat's setback." );
  }

  // GADRAS scales its own efficiency by this; a Monte Carlo of the geometry
  //  will not reproduce that scaling.
  if( (dat.effScalar() > 0.0f) && (std::fabs(dat.effScalar() - 1.0f) > 1.0e-3f) )
    warnings.push_back( "Detector.dat applies an efficiency scalar of "
                        + SpecUtils::printCompact(dat.effScalar(), 4)
                        + ", which a Monte Carlo of the geometry does not"
                        " reproduce; expect the computed efficiency to differ by"
                        " about that factor." );

  relaxGeometryFeatures( gd, warnings );

#if( PERFORM_DEVELOPER_CHECKS )
  // The whole distance convention rests on this: GADRAS measures its source
  //  distance from the reference point, and the setback is how far the crystal
  //  sits behind it.  If the front stack does not sum to the setback, every
  //  endcap-referenced distance is silently wrong.
  {
    const double got = gd.endcap_front_offset_cm();
    const double expect = std::max( setback_cm, used_front );
    if( std::fabs( got - expect ) > (1.0e-4 + 1.0e-3*expect) )
    {
      char msg[256];
      snprintf( msg, sizeof(msg), "endcap_front_offset (%.5f cm) != GADRAS setback"
                " (%.5f cm) - the crystal recess is wrong.", got, expect );
      log_developer_error( __func__, msg );
    }

    const vector<ceelo::GeometryProblem> problems = gd.problems();
    for( const ceelo::GeometryProblem problem : problems )
      log_developer_error( __func__, ("buildGadrasGeometry produced a descriptor"
                    " with problem: " + std::string(ceelo::to_string(problem))).c_str() );
  }
#endif //PERFORM_DEVELOPER_CHECKS

  return gd;
}//buildGadrasGeometry(...)


void setLegacyEfficiencyFromResponse( DetectorPeakResponse &drf,
                    const std::shared_ptr<const ceelo::DetectorResponse> &response,
                    const size_t num_points )
{
  if( !response )
    throw runtime_error( "setLegacyEfficiencyFromResponse: null response." );

  const double a_cm = response->transverse_half_extent();
  if( a_cm <= 0.0 )
    throw runtime_error( "setLegacyEfficiencyFromResponse: response has no extent." );

  double e_lo = response->provenance.valid_e_min_keV;
  double e_hi = response->provenance.valid_e_max_keV;
  if( !(e_hi > e_lo) || (e_lo <= 0.0) )
  {
    e_lo = 35.0;      //the generator's own defaults, when provenance is unset
    e_hi = 3000.0;
  }

  // Log-spaced, plus a pair either side of each crystal K-edge: the stored
  //  response is segmented at the edges, and a curve fit through points that
  //  straddle one would smooth away a real discontinuity.
  vector<double> energies;
  energies.reserve( num_points + 8 );
  const size_t n = std::max<size_t>( 8, num_points );
  for( size_t i = 0; i < n; ++i )
  {
    const double f = static_cast<double>(i) / (n - 1);
    energies.push_back( std::exp( std::log(e_lo) + f*(std::log(e_hi) - std::log(e_lo)) ) );
  }

  for( const double edge : response->descriptor.crystal_k_edges( e_lo, e_hi ) )
  {
    energies.push_back( 0.995*edge );
    energies.push_back( 1.005*edge );
  }
  std::sort( begin(energies), end(energies) );

  // Exactly what DetectorPeakResponse::intrinsicEfficiencyEval does for a
  //  CeeLo-backed DRF: on axis, in the response's own far field, over the same
  //  disk solid angle.  Any divergence here would put the stored curve and the
  //  live query at odds.
  const double d_cm = std::max( 1000.0 * a_cm, 100.0 );
  const double omega = ceelo::disk_solid_angle_fraction( d_cm, a_cm );
  if( omega <= 0.0 )
    throw runtime_error( "setLegacyEfficiencyFromResponse: degenerate solid angle." );

  vector<DetectorPeakResponse::EnergyEffPoint> points;
  points.reserve( energies.size() );
  for( const double energy : energies )
  {
    const ceelo::EffResult res = response->eps_fep( energy, 0.0, 0.0, d_cm );
    const double eff = res.value / omega;
    if( (eff <= 0.0) || IsInf(eff) || IsNan(eff) )
      continue;

    DetectorPeakResponse::EnergyEffPoint p;
    p.energy = static_cast<float>( energy );
    p.efficiency = static_cast<float>( eff );
    if( res.sigma > 0.0 )
      p.efficiencyUncert = static_cast<float>( res.sigma / omega );
    points.push_back( p );
  }//for( const double energy : energies )

  if( points.size() < 2 )
    throw runtime_error( "setLegacyEfficiencyFromResponse: the response gave fewer"
                         " than two positive efficiencies." );

  // setEfficiencyPoints resets the setback, the flags and the DRF source along
  //  with the curve, so carry them across ourselves.
  const double setback = drf.detectorSetback();
  const DetectorPeakResponse::DrfSource source = drf.drfSource();

  drf.setEfficiencyPoints( points, static_cast<float>( 2.0 * a_cm * PhysicalUnits::cm ),
                          -1.0, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic );

  if( setback > 0.0 )
    drf.setDetectorSetback( setback );
  drf.setDrfSource( source );
}//setLegacyEfficiencyFromResponse(...)


std::string genericAttenuatorName( const double atomic_number,
                                   const double areal_density_g_cm2 )
{
  char buf[64] = { '\0' };
  snprintf( buf, sizeof(buf), "AN=%.4g, AD=%.4g g/cm2",
            atomic_number, areal_density_g_cm2 );
  return buf;
}//genericAttenuatorName(...)


bool parseGenericAttenuatorName( const std::string &text, double &atomic_number,
                                 double &areal_density_g_cm2 )
{
  // "AN=13.2, AD=1.35 g/cm2" - tolerant of spacing, case, and a missing unit.
  string t = text;
  SpecUtils::trim( t );
  SpecUtils::to_lower_ascii( t );

  const size_t an_pos = t.find( "an" );
  const size_t ad_pos = t.find( "ad" );
  if( (an_pos == string::npos) || (ad_pos == string::npos) || (ad_pos < an_pos) )
    return false;

  auto value_after = []( const string &s, size_t from, double &out ) -> bool {
    const size_t eq = s.find( '=', from );
    if( eq == string::npos )
      return false;
    size_t i = eq + 1;
    while( (i < s.size()) && std::isspace( static_cast<unsigned char>(s[i]) ) )
      i += 1;
    const size_t start = i;
    while( (i < s.size())
           && (std::isdigit( static_cast<unsigned char>(s[i]) ) || (s[i] == '.')
               || (s[i] == 'e') || (s[i] == '+') || (s[i] == '-')) )
      i += 1;
    if( i == start )
      return false;
    return !!(stringstream( s.substr(start, i - start) ) >> out);
  };//value_after lambda

  double an = 0.0, ad = 0.0;
  if( !value_after( t, an_pos, an ) || !value_after( t, ad_pos, ad ) )
    return false;

  if( (an < 1.0) || (an > 92.0) || (ad <= 0.0) )
    return false;

  atomic_number = an;
  areal_density_g_cm2 = ad;
  return true;
}//parseGenericAttenuatorName(...)

}//namespace CeeLoUtils
