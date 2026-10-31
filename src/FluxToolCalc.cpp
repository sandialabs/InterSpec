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

#define _USE_MATH_DEFINES

#include "InterSpec_config.h"

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

#include "io/DetectorResponse.h"

#include "SpecUtils/SpecFile.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/FluxToolCalc.h"
#include "InterSpec/DrfModifyCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace FluxToolCalc
{

double crystal_setback( const DetectorPeakResponse &det )
{
  if( det.isFixedGeometry() )
    return 0.0;

  const shared_ptr<const ceelo::DetectorResponse> response = det.ceeloResponse();
  if( response )
    return response->descriptor.endcap_front_offset_cm() * PhysicalUnits::cm;

  return det.detectorSetback();
}//crystal_setback(...)


double crystal_diameter( const DetectorPeakResponse &det )
{
  if( det.isFixedGeometry() )
    return 0.0;

  // The crystal disk is what DetectorPeakResponse::intrinsicEfficiencyEval quotes a responses
  //  intrinsic efficiency against, so the two decompositions agree.
  const shared_ptr<const ceelo::DetectorResponse> response = det.ceeloResponse();
  if( response )
    return 2.0 * CeeLoUtils::crystalHalfExtent( response->descriptor ) * PhysicalUnits::cm;

  return det.detectorDiameter();
}//crystal_diameter(...)


double geometric_efficiency( const DetectorPeakResponse &det, const double distance )
{
  if( det.isFixedGeometry() )
    return 1.0;

  return DetectorPeakResponse::fractionalSolidAngle( crystal_diameter(det), distance + crystal_setback(det) );
}//geometric_efficiency(...)


bool outside_energy_range( const DetectorPeakResponse &det, const double energy )
{
  const double lower = det.lowerEnergy(), upper = det.upperEnergy();
  if( !(upper > lower) )
    return false;

  return (energy < lower) || (energy > upper);
}//outside_energy_range(...)


PeakFlux peak_flux( const PeakDef &peak, const DetectorPeakResponse &det,
                    const SpecUtils::Measurement &displayed_foreground, const Options &options )
{
  // The arithmetic below deliberately matches what the Flux Tool did before this function existed,
  //  operation for operation, so with no options turned on, results are identical to the last bit
  //  (test_FluxTool checks this).  Corrections are branched around, never multiplied in as 1.0.
  const float live_time = displayed_foreground.live_time();
  if( !(live_time > 0.0f) || std::isinf(live_time) )
    throw runtime_error( "FluxToolCalc::peak_flux: invalid live time" );

  const bool fixed_geom = det.isFixedGeometry();
  const double distance = options.distance;
  if( !fixed_geom && (!(distance >= 0.0) || std::isinf(distance)) )
    throw runtime_error( "FluxToolCalc::peak_flux: invalid distance" );

  PeakFlux answer;

  const double energy = peak.mean();
  const double amp = peak.peakArea();
  const double ampUncert = peak.peakAreaUncert();

  double cps = amp / live_time;
  double cpsUncert = ampUncert / live_time;

  if( options.scaled_background_peaks && !options.scaled_background_peaks->empty() )
  {
    const GammaInteractionCalc::BackgroundOverlap back
                     = GammaInteractionCalc::overlapping_background( peak, *options.scaled_background_peaks );
    if( back.counts > 0.0 )
    {
      // Background areas were scaled to the foreground live time, so these are foreground CPS.
      answer.background_subtracted = true;
      answer.background_cps = back.counts / live_time;
      answer.background_cps_stat = std::sqrt( back.uncert2 ) / live_time;
      cps = (amp - back.counts) / live_time;
      cpsUncert = std::sqrt( ampUncert*ampUncert + back.uncert2 ) / live_time;
    }//if( back.counts > 0.0 )
  }//if( subtract background )

  const double geomEff = geometric_efficiency( det, distance );
  const double totaleff = det.efficiency( energy, distance );  //distance ignored for fixed geometry

  // The intrinsic efficiency at this distance - per gamma crossing the crystal face - so the intrinsic
  //  and geometric columns multiply to the efficiency used.  Exactly the far-field curve unless the
  //  DRF has a geometry-modeled response (which accounts for where in the crystal gammas interact).
  const double intrinsic = (geomEff > 0.0) ? (totaleff / geomEff) : 0.0;

  // The efficiency uncertainty the DRF states, plus any geometry-model part (never an assumed one).
  const DrfModifyCalc::UncertSummary drf_uncert = DrfModifyCalc::uncertSummary( det, energy, distance );
  const double eff_uncert = drf_uncert.usedFrac( false );

  answer.energy = energy;
  answer.cps = cps;
  answer.cps_stat = cpsUncert;
  answer.geom_eff = geomEff;
  answer.intrinsic_eff = intrinsic;
  answer.intrinsic_syst = intrinsic * eff_uncert;
  answer.eff_uncert_frac = eff_uncert;
  answer.eff_uncert_states_own = drf_uncert.statesOwn;
  answer.outside_drf_energy_range = outside_energy_range( det, energy );

  if( totaleff <= 0.0 || intrinsic <= 0.0 )
  {
    answer.flux_on_det  = std::numeric_limits<double>::infinity();
    answer.flux_per_cm2 = std::numeric_limits<double>::infinity();
    answer.gammas_4pi   = std::numeric_limits<double>::infinity();
    return answer;
  }//if( totaleff <= 0.0 || intrinsic <= 0.0 )

  //gammas into 4pi
  const double gammaInto4pi = cps / totaleff;
  answer.gammas_4pi = gammaInto4pi;
  answer.gammas_4pi_stat = cpsUncert / totaleff;
  answer.gammas_4pi_syst = std::fabs(gammaInto4pi) * eff_uncert;  //fabs: net CPS can be negative

  // A fixed-geometry efficiency already describes the whole source-detector setup, so there is no
  //  distance, flux at the detector, or air path to speak of.
  if( fixed_geom )
    return answer;

  const double fluxOnDet = cps / intrinsic;

  //Flux in g/cm2/s, at the detector face; computed from the gammas into 4pi *before* any air
  //  correction, as it is what actually reaches the detector.
  const double distance_cm = distance / PhysicalUnits::cm;
  const double area_4pi = 4*M_PI*distance_cm*distance_cm;
  const double flux = gammaInto4pi / area_4pi;

  answer.flux_on_det = fluxOnDet;
  answer.flux_on_det_stat = cpsUncert / intrinsic;
  answer.flux_on_det_syst = std::fabs(fluxOnDet) * eff_uncert;

  answer.flux_per_cm2 = flux;
  answer.flux_per_cm2_stat = cpsUncert / totaleff / area_4pi;
  answer.flux_per_cm2_syst = std::fabs(flux) * eff_uncert;

  if( options.correct_for_air )
  {
    // Air path is source to detector face (not on to the crystal), as for every other air
    //  correction in InterSpec, e.g. DetectorPeakResponse::absoluteToIntrinsicMultiple.  The
    //  efficiency itself includes no air: far-field curves have had it removed, and geometry-modeled
    //  responses are computed in vacuum.
    const double mu = GammaInteractionCalc::transmission_coefficient_air( energy, distance );
    const double transmission = std::exp( -mu );
    if( (transmission > 0.0) && !std::isinf(transmission) && !std::isnan(transmission) )
    {
      answer.air_transmission = transmission;
      answer.gammas_4pi /= transmission;
      answer.gammas_4pi_stat /= transmission;
      answer.gammas_4pi_syst /= transmission;
    }
  }//if( options.correct_for_air )

  return answer;
}//peak_flux(...)

}//namespace FluxToolCalc
