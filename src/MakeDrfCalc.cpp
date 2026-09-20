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
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace MakeDrfCalc
{

std::vector<MakeDrfFit::EffFitPoint> intrinsicFitPoints( const MeasuredDrfPoints &points,
                                                         const GeometryChoice &geom,
                                                         const bool inMeV )
{
  const bool use_kernel = (!geom.fixedGeometry && geom.geometry);
  if( !geom.fixedGeometry && !use_kernel && (geom.diameter <= 0.0) )
    throw runtime_error( "A detector diameter is needed to interpret the measured points." );

  std::unique_ptr<CeeLoUtils::GeometryKernel> kernel;
  if( use_kernel )
    kernel.reset( new CeeLoUtils::GeometryKernel( *geom.geometry ) );

  // Flat disk: geometric factor and its distance slope (per cm), central difference
  auto flat_factor = [&geom]( const double dist ) -> double {
    return DetectorPeakResponse::fractionalSolidAngle( geom.diameter, dist + geom.setback );
  };

  vector<MakeDrfFit::EffFitPoint> answer;
  for( const MeasuredEffPoint &p : points.points() )
  {
    if( (p.energy <= 0.0f) || (p.efficiency <= 0.0f) )
      continue;

    double g = 1.0, slope_per_cm = 0.0;
    if( !geom.fixedGeometry )
    {
      if( p.distance < 0.0f )
        continue;  //no geometry to divide out

      if( use_kernel )
      {
        const double d_cm = p.distance / PhysicalUnits::cm;
        g = kernel->intrinsicFactor( p.energy, d_cm );
        if( p.distanceUncert > 0.0f )
          slope_per_cm = kernel->intrinsicFactorDistanceSlope( p.energy, d_cm );
      }else
      {
        g = flat_factor( p.distance );
        if( p.distanceUncert > 0.0f )
        {
          const double delta = std::max( 0.05*PhysicalUnits::cm, 1.0E-3*p.distance );
          const double lo = flat_factor( std::max( 0.0, p.distance - delta ) );
          const double hi = flat_factor( p.distance + delta );
          if( (lo > 0.0) && (hi > 0.0) )
            slope_per_cm = log( hi / lo ) / ((2.0*delta) / PhysicalUnits::cm);
        }
      }//if( use_kernel ) / else
    }//if( !geom.fixedGeometry )

    if( (g <= 0.0) || std::isnan(g) || std::isinf(g) )
      continue;

    MakeDrfFit::EffFitPoint fp;
    fp.energy = inMeV ? (p.energy / 1000.0f) : p.energy;
    fp.efficiency = static_cast<float>( p.efficiency / g );
    fp.fracStatUncert = p.fracStatUncert;
    fp.fracCertUncert = p.fracCertUncert;
    fp.fracDistUncert = static_cast<float>( fabs(slope_per_cm) * (p.distanceUncert / PhysicalUnits::cm) );
    fp.sourceKey = p.sourceKey;
    answer.push_back( std::move(fp) );
  }//for( const MeasuredEffPoint &p : points.points() )

  return answer;
}//intrinsicFitPoints(...)


std::shared_ptr<DetectorPeakResponse> assembleDrf( const std::string &name,
                                                   const std::string &description,
                                                   const MeasuredDrfPoints &points,
                                                   const GeometryChoice &geom,
                                                   const FitResults &fit,
                                                   std::shared_ptr<ceelo::DetectorResponse> mcResponse )
{
  if( fit.eff.coefs.empty() )
    throw runtime_error( "Equation coefficients are empty." );

  for( const float val : fit.eff.coefs )
  {
    if( std::isnan(val) || std::isinf(val) )
      throw runtime_error( "An equation coefficient is invalid." );
  }

  const bool far_field = ((fit.geometryType == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic)
                          || (fit.geometryType == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute));
  if( far_field && ((geom.diameter <= 0.0) || std::isnan(geom.diameter) || std::isinf(geom.diameter)) )
    throw runtime_error( "Detector diameter entered is not a valid distance." );

  auto drf = make_shared<DetectorPeakResponse>( name, description );

  const float eqnEnergyUnits = fit.effInMeV ? 1000.0f : 1.0f;
  drf->fromExpOfLogPowerSeries( fit.eff.coefs, fit.eff.uncerts, 0.0,
                                far_field ? static_cast<float>(geom.diameter) : 0.0f,
                                eqnEnergyUnits, fit.lowerEnergy, fit.upperEnergy, fit.geometryType );
  drf->setDrfSource( DetectorPeakResponse::DrfSource::UserCreatedDrf );

  // The flat-disk setback: entered directly, or implied by the geometry's front layers
  double setback = geom.setback;
  if( far_field && geom.geometry )
    setback = geom.geometry->endcap_front_offset_cm() * PhysicalUnits::cm;
  if( far_field && (setback > 0.0) )
    drf->setDetectorSetback( setback );

  if( !fit.fwhm.coefs.empty() && (fit.fwhmForm != DetectorPeakResponse::kNumResolutionFnctForm) )
    drf->setFwhmCoefficients( fit.fwhm.coefs, fit.fwhmForm, fit.fwhm.uncerts );

  if( fit.peakFitPrefs )
    drf->setPeakFitDetPrefs( fit.peakFitPrefs );

  // Raw points (+ source table) for provenance and re-grounding, and the efficiency uncertainty:
  //  the fit's coefficient covariance (what the curve's uncertainty actually is), with the points'
  //  own block covariance kept alongside as the fallback for other representations.
  shared_ptr<DetectorEfficiencyUncert> uncert;
  if( !points.empty() )
  {
    auto pts = make_shared<MeasuredDrfPoints>( points );
    drf->setMeasuredPoints( pts );
    uncert = pts->toEfficiencyUncert();
  }
  if( !uncert )
    uncert = make_shared<DetectorEfficiencyUncert>();

  const size_t ncoef = fit.eff.coefs.size();
  if( fit.eff.covRowMajor.size() == ncoef*ncoef )
  {
    // A covariance the setter refuses (not a possible set of errors) must not take the whole DRF
    //  down with it - the points' own block covariance is still a usable statement.
    try
    {
      uncert->setCoefficientCovariance( fit.eff.covRowMajor );
    }catch( std::exception &e )
    {
      cerr << "MakeDrfCalc::assembleDrf: not storing the fit covariance: " << e.what() << endl;
    }
  }//if( the fit produced a coefficient covariance )

  if( !uncert->isEmpty() )
    drf->setEfficiencyUncert( uncert );

  if( geom.geometry && far_field )
    drf->setGeometry( geom.geometry );

  // The CeeLo response: a generated (MC) one grounded to the raw points, else the instant
  //  curve-transfer anchored on the fitted curve with its covariance.  Never both - the MC
  //  grounding carries the raw-point covariance, the transfer carries the fit covariance.
  if( mcResponse && far_field )
  {
    try
    {
      // Re-ground only a response that ARRIVED grounded: the MC tool grounds (or deliberately does
      //  not, per its "correct to this detector's measured efficiency" checkbox) before handing the
      //  response over, and re-grounding an ungrounded one would override that choice.  We still
      //  redo it, because the measured points may have changed since it was generated.
      const bool was_grounded = !mcResponse->grounding.empty();

      // k(E) is the ratio of measurement to the UNGROUNDED model, so the old grounding has to go
      //  before any model efficiency is evaluated - `ground_to_points` clears it too, but only
      //  after honoring a pre-filled `model_eff`, which would otherwise already carry exp(ln_k)
      //  and drive the new k(E) to 1 (silently cancelling the grounding).  Clearing here also
      //  makes this whole block idempotent, which matters because one save calls assembleDrf
      //  several times (store, N42 export, reference sheet) on this same shared response.
      mcResponse->grounding = ceelo::GroundingBlock();

      const ceelo::GeometryDescriptor &gd = mcResponse->descriptor;
      vector<ceelo::GroundingPoint> ground_pts;
      for( const MeasuredEffPoint &p : points.points() )
      {
        if( (p.distance < 0.0f) || (p.efficiency <= 0.0f) )
          continue;
        ceelo::GroundingPoint gp;
        gp.energy_keV = p.energy;
        gp.measured_eff = p.efficiency;
        gp.frac_stat_sigma = p.fracStatUncert;
        gp.frac_cert_sigma = p.fracCertUncert;
        gp.source_key = p.sourceKey;
        gp.distance_cm = p.distance / PhysicalUnits::cm;
        gp.cos_theta = 1.0;  //Create DRF sources are on-axis
        // Model efficiency at the point's own (face-referenced) distance
        gp.model_eff = mcResponse->eps_fep_at( p.energy,
                          CeeLoUtils::sourcePositionFromFace( gd, 0.0, 0.0, gp.distance_cm ) ).value;
        ground_pts.push_back( std::move(gp) );
      }//for( const MeasuredEffPoint &p : points.points() )

      if( was_grounded && !ground_pts.empty() )
        ceelo::ResponseGenerator::ground_to_points( *mcResponse, ground_pts, /*curve_derived=*/false );

      drf->setCeeloResponse( mcResponse );
    }catch( std::exception &e )
    {
      cerr << "MakeDrfCalc::assembleDrf: failed to ground MC response: " << e.what() << endl;
      drf->setCeeloResponse( mcResponse );  //attach ungrounded
    }
  }else if( geom.geometry && far_field )
  {
    try
    {
      const CeeLoUtils::TransferAnchor anchor
                  = CeeLoUtils::curveAnchorWithCovarianceForDrf( drf, *geom.geometry, -1.0 );
      const ceelo::AnchorCurve tot_curve = CeeLoUtils::totalTransferAnchorForDrf( drf, anchor );
      drf->setCeeloResponse( CeeLoUtils::makeTransferResponse( *geom.geometry, anchor, tot_curve, name ) );
    }catch( std::exception &e )
    {
      cerr << "MakeDrfCalc::assembleDrf: could not build the curve-transfer response: " << e.what() << endl;
    }
  }//if( mcResponse ) / else if( geometry )

  if( !drf->isValid() )
    throw runtime_error( "DRF wasnt valid after creation" );

  return drf;
}//assembleDrf(...)


GeometryChoice geometryChoiceForDrf( const DetectorPeakResponse &drf )
{
  GeometryChoice answer;
  answer.fixedGeometry = drf.isFixedGeometry();
  if( answer.fixedGeometry )
    return answer;

  answer.geometry = drf.geometry();
  answer.diameter = drf.detectorDiameter();
  answer.setback = drf.detectorSetback();

  return answer;
}//geometryChoiceForDrf(...)


void refitEfficiencyFromPoints( DetectorPeakResponse &drf,
                                const MeasuredDrfPoints &points,
                                const int nterms,
                                std::string &warnings )
{
  warnings.clear();

  if( points.empty() )
    throw runtime_error( "There are no measured points to fit." );

  // The fit is of the INTRINSIC efficiency (intrinsicFitPoints divides each point's source geometry
  //  out), so installing it on a DRF whose curve is read as absolute efficiency would be wrong by a
  //  solid angle.  Such a detector has to be re-characterized, not re-fit in place.
  if( drf.geometryType() == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute )
    throw runtime_error( "This detector states ABSOLUTE efficiency at a fixed distance, which a"
                         " re-fit of intrinsic efficiency cannot express." );

  const shared_ptr<const DetectorEfficiencyCurve> curve = drf.efficiencyCurve();

  // The equation keeps the units it was written in, so its coefficients stay comparable with the
  //  ones the user (or a previous fit) can see.
  const float units = (curve && curve->isValid() && (curve->energyUnits() > 0.0f))
                      ? curve->energyUnits() : static_cast<float>(PhysicalUnits::keV);
  const bool inMeV = (units > 10.0f);

  int order = nterms;
  if( order <= 0 )
    order = curve ? static_cast<int>( curve->expOfLogPowerSeriesCoeffs().size() ) : 0;
  if( order <= 0 )
    order = 4;  //a DRF that was not an equation before (a pairs curve being re-fit) gets the default

  const GeometryChoice geom = geometryChoiceForDrf( drf );
  const vector<MakeDrfFit::EffFitPoint> effpts = intrinsicFitPoints( points, geom, inMeV );

  if( static_cast<int>(effpts.size()) < 2 )
    throw runtime_error( "At least two usable measured points are needed to fit an efficiency"
                         " equation (points need a positive energy, efficiency, and - for a"
                         " far-field detector - a source distance)." );

  // More terms than points would be an under-determined fit; fall back rather than refuse, since the
  //  term count is not something the user chose here.
  if( static_cast<int>(effpts.size()) < order )
    order = static_cast<int>( effpts.size() );

  const MakeDrfFit::EffFitResult fit = MakeDrfFit::performEfficiencyFit( effpts, order );
  warnings = fit.warnings;

  if( fit.coefs.empty() )
    throw runtime_error( "The efficiency fit returned no coefficients." );

  for( const float val : fit.coefs )
  {
    if( std::isnan(val) || std::isinf(val) )
      throw runtime_error( "The efficiency fit returned an invalid coefficient." );
  }

  // The energy range the points actually cover (they are sorted by energy).
  const float lower = points.points().front().energy;
  const float upper = points.points().back().energy;

  // Build the whole new state before touching `drf`, so a throw leaves it as it was.
  auto new_curve = make_shared<DetectorEfficiencyCurve>();
  new_curve->setFromExpOfLogPowerSeries( fit.coefs, fit.uncerts, units );

  auto pts = make_shared<MeasuredDrfPoints>( points );

  // The two uncertainty stores this DRF ends up with, from this one fit: the coefficient covariance
  //  (what an equation curve's uncertainty IS - see DetectorEfficiencyCurve::fracCovariance), and
  //  the points' own block covariance beside it as the fallback for other representations.
  shared_ptr<DetectorEfficiencyUncert> uncert = pts->toEfficiencyUncert();
  if( !uncert )
    uncert = make_shared<DetectorEfficiencyUncert>();

  const size_t ncoef = fit.coefs.size();
  bool have_coef_cov = (fit.covRowMajor.size() == (ncoef*ncoef));
  if( have_coef_cov )
  {
    try
    {
      uncert->setCoefficientCovariance( fit.covRowMajor );
    }catch( std::exception & )
    {
      have_coef_cov = false;  //not a possible set of errors; the points' covariance still stands
    }
  }//if( the fit produced a coefficient covariance )

  if( !have_coef_cov )
  {
    warnings += (warnings.empty() ? "" : "  ");
    warnings += "The fit did not produce a usable coefficient covariance, so the equation carries"
                " only the uncertainty the measured points imply.";
  }

  new_curve->setUncertainty( uncert->isEmpty() ? nullptr : uncert );

  drf.replaceEfficiencyCurve( new_curve );
  drf.setMeasuredPoints( pts );
  if( upper > lower )
    drf.setEnergyRange( lower, upper );

  // A point the fitted curve cannot get near is almost always a mis-typed cell rather than a bad
  //  fit, and it drags the whole curve with it - so name it.  A fit does not pass through its
  //  points, hence the generous factor before anything is said.
  {
    double worst_ratio = 1.0;
    float worst_energy = 0.0f;
    for( const MakeDrfFit::EffFitPoint &fitpt : effpts )
    {
      if( (fitpt.efficiency <= 0.0f) || (fitpt.energy <= 0.0f) )
        continue;

      const float energy_kev = inMeV ? (1000.0f*fitpt.energy) : fitpt.energy;
      const double curve_eff = new_curve->efficiency( energy_kev );
      if( curve_eff <= 0.0 )
        continue;

      const double ratio = curve_eff / static_cast<double>( fitpt.efficiency );
      const double off_by = (ratio > 1.0) ? ratio : (1.0/ratio);
      if( off_by > worst_ratio )
      {
        worst_ratio = off_by;
        worst_energy = energy_kev;
      }
    }//for( const MakeDrfFit::EffFitPoint &fitpt : effpts )

    if( worst_ratio > 2.0 )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "the fitted curve is a factor of %.3g from the point at"
                " %.1f keV - check that value", worst_ratio, worst_energy );
      warnings += (warnings.empty() ? "" : "  ");
      warnings += buffer;
    }
  }
}//refitEfficiencyFromPoints(...)

}//namespace MakeDrfCalc
