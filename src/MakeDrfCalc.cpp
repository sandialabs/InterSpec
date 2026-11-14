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
#include <set>
#include <cmath>
#include <cfloat>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "SandiaDecay/SandiaDecay.h"

#include "io/DetectorResponse.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/MakeDrfCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

using namespace std;

namespace
{
  /** The energy a peak's gammas are attributed at: its assigned gamma's, or with none (e.g., a peak
   of a source the user gave a branching ratio for) its mean. */
  double peak_gamma_energy( const PeakDef &peak )
  {
    if( peak.parentNuclide() || peak.xrayElement() || peak.reaction() )
    {
      try
      {
        return peak.gammaParticleEnergy();
      }catch( std::exception & )
      {
        //A source without the gamma resolved; use the mean
      }
    }//if( the peak has a source )

    return peak.mean();
  }//peak_gamma_energy(...)
}//namespace


namespace MakeDrfCalc
{

std::vector<MakeDrfFit::EffFitPoint> intrinsicFitPoints( const MeasuredDrfPoints &points,
                                                         const GeometryChoice &geom,
                                                         const bool inMeV )
{
  const bool use_kernel = (!geom.fixedGeometry && geom.geometry);
  if( !geom.fixedGeometry && !use_kernel && (geom.diameter <= 0.0) )
    throw runtime_error( "A detector diameter is needed to interpret the measured points." );

  // The kernel quotes the far field against the diameter the fitted curve will be stored with, so
  //  the DRF's flat-disk evaluation reads the points back (for Make Detector Response's geometry
  //  mode that is the object disk; for a re-fit of an existing DRF, whatever diameter it states).
  std::unique_ptr<CeeLoUtils::GeometryKernel> kernel;
  if( use_kernel )
    kernel.reset( new CeeLoUtils::GeometryKernel( *geom.geometry, 2048, geom.diameter / PhysicalUnits::cm ) );

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
    fp.efficiency = static_cast<float>( p.efficiency / (g * geom.efficiencyMultiple) );
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

  if( !(geom.efficiencyMultiple > 0.0) || !std::isfinite(geom.efficiencyMultiple) )
    throw runtime_error( "The efficiency multiple must be greater than zero." );

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

  // Before any response is made: they are given it (the fitted equation is one element's).
  drf->setEfficiencyMultiple( geom.efficiencyMultiple );

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
  // Only a response for the geometry being recorded: one generated for an earlier geometry, or
  //  before switching to the diameter-only flat disk, describes some other detector.
  if( mcResponse && geom.geometry && far_field )
  {
    // Re-ground only a response that ARRIVED grounded: the MC tool grounds (or deliberately does
    //  not, per its "correct to this detector's measured efficiency" checkbox) before handing the
    //  response over, and re-grounding an ungrounded one would override that choice.  We still
    //  redo it, because the measured points may have changed since.
    const bool was_grounded = !mcResponse->grounding.empty();

    // The multiple may have been changed since the Monte Carlo ran; k(E) is fit with the scale in
    //  the model, so the response has to carry the DRF's before grounding (and before attaching,
    //  either way - the DRF already holds the multiple, and attaching a response adopts its scale).
    const bool rescale = (mcResponse->efficiency_scale != geom.efficiencyMultiple);

    // Ground/rescale a COPY: the caller's response is attached to every DRF stored or exported from
    //  it, so changing it in place would change DRFs already in use (and race threads evaluating
    //  them).  DetectorResponse is not copyable, so copy through its codec (as
    //  CeeLoUtils::withoutGrounding does) - only when there is something to change, as that is a
    //  full serialize and parse.
    shared_ptr<ceelo::DetectorResponse> response;
    if( was_grounded || rescale )
    {
      try
      {
        response = ceelo::DetectorResponse::from_xml_string( mcResponse->to_xml_string() );
      }catch( std::exception &e )
      {
        cerr << "MakeDrfCalc::assembleDrf: failed to copy MC response: " << e.what() << endl;
      }

      if( response )
        response->efficiency_scale = geom.efficiencyMultiple;
    }//if( was_grounded || rescale )

    if( !response )
    {
      // As handed over, unmodified - unless it carries another multiple, which would make every
      //  efficiency off by their ratio; then the DRF goes without a response (the curve stands).
      if( !rescale )
        drf->setCeeloResponse( mcResponse );
      else
        cerr << "MakeDrfCalc::assembleDrf: not attaching the MC response, as it could not be given"
                " the efficiency multiple." << endl;
    }else if( was_grounded )
    {
      try
      {
        // k(E) is the ratio of measurement to the UNGROUNDED model, so the old grounding has to go
        //  before any model efficiency is evaluated - `ground_to_points` clears it too, but only
        //  after honoring a pre-filled `model_eff`, which would otherwise already carry exp(ln_k)
        //  and drive the new k(E) to 1 (silently cancelling the grounding).
        response->grounding = ceelo::GroundingBlock();

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
          gp.distance_cm = p.distance / PhysicalUnits::cm;  //from the face; see groundResponseToFacePoints
          gp.cos_theta = 1.0;  //Create DRF sources are on-axis
          ground_pts.push_back( std::move(gp) );
        }//for( const MeasuredEffPoint &p : points.points() )

        if( !ground_pts.empty() )
          CeeLoUtils::groundResponseToFacePoints( *response, ground_pts, /*curve_derived=*/false );
      }catch( std::exception &e )
      {
        cerr << "MakeDrfCalc::assembleDrf: failed to ground MC response: " << e.what() << endl;
        response->grounding = ceelo::GroundingBlock();  //attach ungrounded
      }

      drf->setCeeloResponse( response );
    }else
    {
      drf->setCeeloResponse( response );  //only rescaled
    }//if( could not copy the response ) / else if( re-ground ) / else
  }else if( geom.geometry && far_field )
  {
    try
    {
      const CeeLoUtils::TransferAnchor anchor
                  = CeeLoUtils::curveAnchorWithCovarianceForDrf( drf, *geom.geometry, -1.0 );
      const ceelo::AnchorCurve tot_curve = CeeLoUtils::totalTransferAnchorForDrf( drf, anchor );
      drf->setCeeloResponse( CeeLoUtils::makeTransferResponse( *geom.geometry, anchor, tot_curve, name,
                                                               drf->efficiencyMultiple() ) );
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
  answer.efficiencyMultiple = drf.efficiencyMultiple();

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


SampleResult samplePoints( const SampleInput &sample,
                           const std::vector<BackgroundPeak> &background,
                           const bool airAttenuation )
{
  typedef GammaInteractionCalc::ShieldingSourceChi2Fcn Chi2Fcn;

  const vector<CalSource> &sources = sample.sources;
  const size_t nsources = sources.size();

  SampleResult answer;
  vector<bool> usable( nsources, false );
  for( size_t i = 0; i < nsources; ++i )
    usable[i] = sources[i].usable;

  // Each gamma line counts toward the nearest of all the spectrum's peaks, used or not.  As
  //  Chi2Fcn::observedPeakEnergyWidths, but a peak without a source assigned is kept, at its mean.
  vector<pair<double,double>> energy_widths;
  for( const SamplePeak &p : sample.peaks )
  {
    if( p.peak )
      energy_widths.emplace_back( peak_gamma_energy( *p.peak ),
                                  p.peak->gausPeak() ? p.peak->sigma() : 0.25*p.peak->roiWidth() );
  }
  std::sort( begin(energy_widths), end(energy_widths) );

  // Per source with a nuclide: the gammas per second (before any attenuation) it puts in each
  //  peak, keyed by the peak's gamma energy (as `energy_widths` is)
  vector<map<double,double>> rates( nsources );
  for( size_t i = 0; i < nsources; ++i )
  {
    const CalSource &src = sources[i];
    if( !usable[i] || !src.nuclide )
      continue;

    try
    {
      const bool decay_during_meas = (sample.realTime > sm_decayDuringMeasFraction*src.nuclide->halfLife);
      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity( src.nuclide, Chi2Fcn::sm_activityUnits );
      Chi2Fcn::cluster_peak_activities( rates[i], energy_widths, mix, src.activity, src.age,
                                        sm_photopeakClusterSigma, -1.0, decay_during_meas,
                                        sample.realTime, nullptr );
    }catch( std::exception &e )
    {
      usable[i] = false;
      answer.sourceErrors.emplace_back( i, e.what() );
    }//try / catch
  }//for( size_t i = 0; i < nsources; ++i )

  // The lines of the unusable sources (at unit activity, and the default age, as theirs may be
  //  what is wrong), and the strongest of each - a peak they add more than a trace to is left out
  vector<vector<SandiaDecay::EnergyRatePair>> unusable_lines( nsources );
  vector<double> unusable_max_rate( nsources, 0.0 );
  for( size_t i = 0; i < nsources; ++i )
  {
    const SandiaDecay::Nuclide * const nuc = sources[i].nuclide;
    if( usable[i] || !nuc )
      continue;

    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, 1.0, PeakDef::defaultDecayTime( nuc, nullptr ) );
    unusable_lines[i] = mix.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
    for( const SandiaDecay::EnergyRatePair &line : unusable_lines[i] )
      unusable_max_rate[i] = std::max( unusable_max_rate[i], line.numPerSecond );
  }//for( size_t i = 0; i < nsources; ++i )

  // The fraction of a source's gammas at `energy` that get through its shielding and the air
  auto transmission = [airAttenuation]( const CalSource &src, const double energy ) -> double {
    double trans = 1.0;
    if( airAttenuation && (src.distance > 0.0) )
      trans *= exp( -GammaInteractionCalc::transmission_length_coefficient_air( energy ) * src.distance );

    if( src.shieldMaterial )
    {
      // A named material attenuates by its actual composition (an effective atomic number misses,
      //  e.g., hydrogen's extra electrons per gram - ~10% of mu for water)
      trans *= exp( -GammaInteractionCalc::transmition_coefficient_material( src.shieldMaterial.get(),
                                                                             energy, src.shieldThickness ) );
    }else if( (src.shieldAtomicNumber >= 1.0) && (src.shieldArealDensity > 0.0) )
    {
      const double mu = MassAttenuation::massAttenuationCoefficientFracAN( src.shieldAtomicNumber, energy );
      trans *= exp( -mu * src.shieldArealDensity );
    }
    return trans;
  };//transmission lambda

  for( size_t peak_index = 0; peak_index < sample.peaks.size(); ++peak_index )
  {
    const SamplePeak &sp = sample.peaks[peak_index];
    if( !sp.peak || (sp.source < 0) )
      continue;

    assert( static_cast<size_t>(sp.source) < nsources );
    if( (static_cast<size_t>(sp.source) >= nsources) || !usable[sp.source] )
      continue;  //the source's own problems say why

    const PeakDef &peak = *sp.peak;
    const CalSource &src = sources[sp.source];

    PeakPoint point;
    point.peak = peak_index;
    point.energy = peak_gamma_energy( peak );  //also the key `energy_widths` (so `rates`) use

    // A peak with no nuclide uses the branching ratio the user entered for it
    if( !src.nuclide && (sp.userBr < 0.0) )
    {
      point.status = PeakPoint::Status::NoUserBranchRatio;
      answer.points.push_back( point );
      continue;
    }

    // A data-defined peak has no sigma; a quarter of its ROI stands in, as in `energy_widths`
    const double peak_sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();

    // Leave out peaks an unusable source would add more than a trace to
    const double match_width = sm_photopeakClusterSigma * peak_sigma;
    for( size_t i = 0; (i < nsources) && (point.interferer < 0); ++i )
    {
      for( const SandiaDecay::EnergyRatePair &line : unusable_lines[i] )
      {
        if( (fabs(line.energy - point.energy) < match_width)
            && (line.numPerSecond > 1.0E-3*unusable_max_rate[i]) )
        {
          point.interferer = static_cast<int>( i );
          break;
        }
      }
    }//for( loop over unusable sources )

    if( point.interferer >= 0 )
    {
      point.status = PeakPoint::Status::NeedsUnusableSource;
      answer.points.push_back( point );
      continue;
    }

    // The background rate is the matching peaks' summed area over the summed live time of the
    //  background spectra they are in - each spectrum counted once, even when two of its peaks
    //  match (e.g., a doublet).
    double bkg_area = 0.0, bkg_area_uncert_sq = 0.0, bkg_lt = 0.0;
    set<size_t> bkg_spectra_used;
    for( const BackgroundPeak &b : background )
    {
      if( !b.peak )
        continue;

      const PeakDef &bpeak = *b.peak;
      if( (fabs(bpeak.mean() - peak.mean()) < 1.5*peak_sigma)
         || (bpeak.parentNuclide() && peak.parentNuclide()
             && (fabs(peak_gamma_energy(bpeak) - point.energy) < 1.0)) )
      {
        bkg_area += bpeak.peakArea();
        bkg_area_uncert_sq += bpeak.peakAreaUncert() * bpeak.peakAreaUncert();
        if( bkg_spectra_used.insert( b.spectrum ).second )
          bkg_lt += b.liveTime;
      }
    }//for( const BackgroundPeak &b : background )

    point.peakArea = peak.peakArea();
    point.peakAreaUncert = peak.peakAreaUncert();

    const bool sub_back = ((bkg_area > DBL_EPSILON) && (bkg_lt > DBL_EPSILON));
    double scaled_back = -1.0, scaled_back_uncert = 0.0;
    if( sub_back )
    {
      point.bkgArea = bkg_area;
      point.bkgLiveTime = bkg_lt;

      scaled_back = sample.liveTime * bkg_area / bkg_lt;
      scaled_back_uncert = scaled_back * sqrt( bkg_area_uncert_sq ) / bkg_area;
      point.peakArea -= scaled_back;
      point.peakAreaUncert = sqrt( point.peakAreaUncert*point.peakAreaUncert
                                   + scaled_back_uncert*scaled_back_uncert );
    }//if( sub_back )

    // Annihilation photons are Doppler broadened, whatever made the positron, so a 511 keV peak is
    //  wider than the detector resolution
    const bool is_annihilation = (peak.sourceGammaType() == PeakDef::SourceGammaType::AnnihilationGamma)
                                 || (peak.gausPeak() && (fabs(peak.mean() - 510.9989) < 0.5*peak.fwhm()));
    point.useForFwhm = (peak.gausPeak() && !is_annihilation);

    // The gammas per second that land in this peak, from each of the spectrum's sources, after
    //  their shielding and the air; their activity uncertainties add in quadrature.
    double rate_sq_uncert = 0.0;
    for( size_t i = 0; i < nsources; ++i )
    {
      if( !usable[i] )
        continue;

      const CalSource &s = sources[i];
      double rate = 0.0;
      if( s.nuclide )
      {
        const auto pos = rates[i].find( point.energy );
        if( pos != end(rates[i]) )
          rate = pos->second;
      }else if( i == static_cast<size_t>(sp.source) )
      {
        rate = sp.userBr * s.activity;
      }

      if( rate <= 0.0 )
        continue;

      rate *= transmission( s, point.energy );
      point.sourceRate += rate;
      rate_sq_uncert += pow( rate * s.fracActivityUncert, 2.0 );
    }//for( size_t i = 0; i < nsources; ++i )

    point.fracSourceUncert = (point.sourceRate > 0.0) ? (sqrt(rate_sq_uncert) / point.sourceRate) : 0.0;

    // Absolute efficiency at this source's own distance
    const double expected = point.sourceRate * sample.liveTime;
    const double eff = point.peakArea / expected;

    if( !(expected > 0.0) || std::isnan(eff) || std::isinf(eff) )
    {
      point.status = PeakPoint::Status::NoSourceGammas;
    }else if( !(point.peakArea > 0.0) )
    {
      point.status = PeakPoint::Status::NoNetArea;
    }else
    {
      // The statistical (peak area) and certificate (source activity) parts are kept separate, plus
      //  the provenance - the fit, a Monte-Carlo grounding, and a later re-fit from the DRF alone all
      //  start from exactly this.
      MeasuredEffPoint &raw = point.measured;
      raw.energy = static_cast<float>( point.energy );
      raw.efficiency = static_cast<float>( eff );
      raw.fracStatUncert = static_cast<float>( point.peakAreaUncert / point.peakArea );
      raw.fracCertUncert = static_cast<float>( point.fracSourceUncert );
      raw.distance = static_cast<float>( src.distance );
      raw.distanceUncert = (src.distance > 0.0) ? static_cast<float>( src.distanceUncert ) : 0.0f;
      raw.peakArea = static_cast<float>( point.peakArea );
      raw.peakAreaUncert = static_cast<float>( point.peakAreaUncert );
      raw.liveTime = static_cast<float>( sample.liveTime );
      if( sub_back )
      {
        raw.bkgPeakArea = static_cast<float>( scaled_back );
        raw.bkgPeakAreaUncert = static_cast<float>( scaled_back_uncert );
      }
    }//if( no gammas ) / else if( no net area ) / else

    answer.points.push_back( point );
  }//for( loop over peaks )

  return answer;
}//samplePoints(...)


int defaultEquationTerms( const int npoints, const int maxTerms )
{
  if( (npoints < 1) || (maxTerms < 1) )
    return 0;
  return std::min( maxTerms, std::max( std::min(npoints, 2), npoints - 2 ) );
}//defaultEquationTerms(...)

}//namespace MakeDrfCalc
