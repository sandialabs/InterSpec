#ifndef FluxToolCalc_h
#define FluxToolCalc_h
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

#include <vector>

class PeakDef;
class DetectorPeakResponse;

namespace SpecUtils
{
  class Measurement;
}

/** The per-peak math of the Flux Tool (FluxToolWidget), kept free of the GUI so it can be unit tested.

 All distances are in PhysicalUnits, and measured from the center of the source to the face of the
 detector (the endcap front - the surface a user can put a ruler against).  The detector crystal sits
 some distance behind that face: the DRFs setback, or for a DRF with a geometry-modeled (CeeLo)
 response, the responses endcap offset; see #crystal_setback.
 */
namespace FluxToolCalc
{
  struct Options
  {
    /** Source center to detector face; not used for fixed-geometry DRFs. */
    double distance = 0.0;

    /** Whether to correct the gammas-into-4pi rate for attenuation in the air between the source and
     detector face.  Not used for fixed-geometry DRFs.
     */
    bool correct_for_air = false;

    /** Background peaks to subtract, already put on the foregrounds footing (see
     GammaInteractionCalc::scale_background_peaks).  Null, or empty, for no subtraction.
     */
    const std::vector<PeakDef> *scaled_background_peaks = nullptr;
  };//struct Options


  /** The Flux Tool results for one peak.

   Where both are given, the `_stat` member is the counting-statistics (peak area, and any
   subtracted background area) uncertainty, and the `_syst` member the detector-efficiency part.
   */
  struct PeakFlux
  {
    double energy = 0.0;

    /** The peak CPS; net of any background subtracted. */
    double cps = 0.0, cps_stat = 0.0;

    /** Whether a background peak was subtracted. */
    bool background_subtracted = false;

    /** The CPS of the background peak(s) subtracted - zero if none. */
    double background_cps = 0.0, background_cps_stat = 0.0;

    /** The fraction of the solid angle the crystal face takes up; 1 for fixed geometry. */
    double geom_eff = 0.0;

    /** Efficiency per gamma crossing the crystal face (the efficiency used divided by #geom_eff). */
    double intrinsic_eff = 0.0, intrinsic_syst = 0.0;

    /** Gammas per second crossing the crystal face; 0 for fixed geometry. */
    double flux_on_det = 0.0, flux_on_det_stat = 0.0, flux_on_det_syst = 0.0;

    /** Gammas per cm2 per second at the detector face (i.e., at #Options::distance); 0 for fixed
     geometry.  Never corrected for air - it is what reaches the detector.
     */
    double flux_per_cm2 = 0.0, flux_per_cm2_stat = 0.0, flux_per_cm2_syst = 0.0;

    /** Gammas per second the source emits into 4pi (corrected for air, if asked to).  For a
     fixed-geometry DRF, the sources emission rate - per cm2, m2, or gram, for those geometry types.
     */
    double gammas_4pi = 0.0, gammas_4pi_stat = 0.0, gammas_4pi_syst = 0.0;

    /** The fraction of gammas transmitted through the air, that #gammas_4pi was divided by; 1 when not
     corrected.
     */
    double air_transmission = 1.0;

    /** The fractional detector-efficiency uncertainty used (see DrfModifyCalc::UncertSummary), and
     whether the DRF states one of its own.
     */
    double eff_uncert_frac = 0.0;
    bool eff_uncert_states_own = false;

    /** Whether the peak is outside the energy range the DRF states it is valid for. */
    bool outside_drf_energy_range = false;
  };//struct PeakFlux


  /** Computes the Flux Tool results for a peak.

   @param peak The foreground peak.
   @param det A valid detector response.
   @param displayed_foreground The displayed foreground histogram, whose live time gives the CPS.  This
          is deliberately a histogram, not a live time: a spectrum file (SpecMeas) is not a
          SpecUtils::Measurement, so the live time of the whole file can not be used by mistake.
   @param options The distance, and optional corrections.

   Throws std::exception if the live time is not positive, or (for far-field DRFs) the distance is
   negative or not finite.
   */
  PeakFlux peak_flux( const PeakDef &peak, const DetectorPeakResponse &det,
                      const SpecUtils::Measurement &displayed_foreground, const Options &options );

  /** How far the crystal face is behind the detector face: for a DRF with a geometry-modeled (CeeLo)
   response, the responses endcap offset; otherwise the DRFs setback.  Zero for fixed geometry.
   */
  double crystal_setback( const DetectorPeakResponse &det );

  /** The crystal diameter the geometric efficiency is computed for: for a DRF with a geometry-modeled
   response, the responses crystal; otherwise the DRFs stated diameter.  Zero for fixed geometry.
   */
  double crystal_diameter( const DetectorPeakResponse &det );

  /** The fraction of the solid angle the crystal face takes up, for a point source `distance` from the
   detector face (the crystal is #crystal_setback further away).  1 for fixed-geometry DRFs.
   */
  double geometric_efficiency( const DetectorPeakResponse &det, const double distance );

  /** Whether `energy` is outside the energy range the DRF states it is valid for.

   A DRF that does not state a range (both zero), or states one that is not a range (upper <= lower,
   as some GADRAS detectors do), is never out of range.
   */
  bool outside_energy_range( const DetectorPeakResponse &det, const double energy );
}//namespace FluxToolCalc

#endif //FluxToolCalc_h
