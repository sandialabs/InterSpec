#ifndef MakeDrfCalc_h
#define MakeDrfCalc_h
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

#include <memory>
#include <string>
#include <vector>

#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/DetectorPeakResponse.h"

class PeakFitDetPrefs;
class MeasuredDrfPoints;
namespace ceelo{ struct GeometryDescriptor; class DetectorResponse; }

/** The non-GUI core of the "Create Detector Response Function" tool: turning measured peak
 efficiency points into the points the efficiency equation is fit to, and assembling the DRF from
 the fit results.  Kept free of Wt so it can be exercised (and validated against Monte Carlo) in
 unit tests.
 */
namespace MakeDrfCalc
{
  /** How the user described the detector. */
  struct GeometryChoice
  {
    /** The full detector geometry when the user entered one; null for the flat-disk description
     (just a diameter, and optionally a setback). */
    std::shared_ptr<const ceelo::GeometryDescriptor> geometry;

    /** Detector (crystal) diameter, PhysicalUnits.  Used by the flat-disk model, and recorded on
     the DRF either way. */
    double diameter = 0.0;

    /** Distance from the detector face to the crystal, PhysicalUnits; flat-disk model only (with a
     full geometry it follows from the layers). */
    double setback = 0.0;

    /** True for the fixed-geometry (e.g. activity per area) efficiency types: the measured points
     already are what the equation describes, no source geometry is divided out. */
    bool fixedGeometry = false;
  };//struct GeometryChoice


  /** The far-field intrinsic efficiency each absolute measured point corresponds to, ready to be
   fit: absolute efficiency divided by the geometric factor of its own source distance - the
   flat-disk solid angle fraction, or with a full geometry the kernel factor of
   CeeLoUtils::GeometryKernel::intrinsicFactor (so points taken at different distances land on one
   curve, near-field effects included).  A point's distance uncertainty becomes a fractional
   efficiency uncertainty through the factor's distance slope, correlated within its source.

   Energies are keV, or MeV when `inMeV`.  Points without a distance are skipped unless
   `geom.fixedGeometry`.  Throws std::runtime_error when the geometry is unusable.
   */
  std::vector<MakeDrfFit::EffFitPoint> intrinsicFitPoints( const MeasuredDrfPoints &points,
                                                           const GeometryChoice &geom,
                                                           const bool inMeV );

  /** Everything #assembleDrf needs besides the points and geometry. */
  struct FitResults
  {
    /** Efficiency equation fit (energy units per `effInMeV`). */
    MakeDrfFit::EffFitResult eff;
    bool effInMeV = true;

    /** Energy range (keV) the equation is valid over - normally the lowest/highest peak used. */
    float lowerEnergy = 0.0f, upperEnergy = 0.0f;

    /** FWHM fit; `fwhm.coefs` empty when there is none. */
    MakeDrfFit::FwhmFitResult fwhm;
    DetectorPeakResponse::ResolutionFnctForm fwhmForm = DetectorPeakResponse::kNumResolutionFnctForm;

    DetectorPeakResponse::EffGeometryType geometryType = DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;

    std::shared_ptr<const PeakFitDetPrefs> peakFitPrefs;
  };//struct FitResults


  /** Builds the DRF: the fitted equation with its full coefficient covariance, the FWHM equation
   with its uncertainties, the raw measured points (and source table) for provenance, the geometry,
   and a CeeLo response - `mcResponse` grounded to the raw points when one was generated, else
   (with a full geometry) the instant curve-transfer response anchored on the fitted curve with its
   covariance (CeeLoUtils::curveAnchorWithCovarianceForDrf).  A flat-disk description gets no
   response.

   Throws std::runtime_error when the fit results cannot make a valid DRF.
   */
  std::shared_ptr<DetectorPeakResponse> assembleDrf( const std::string &name,
                                                     const std::string &description,
                                                     const MeasuredDrfPoints &points,
                                                     const GeometryChoice &geom,
                                                     const FitResults &fit,
                                                     std::shared_ptr<ceelo::DetectorResponse> mcResponse );


  /** The geometry `drf` describes itself with, ready to hand to #intrinsicFitPoints or
   #refitEfficiencyFromPoints: its geometry descriptor (when it has one and is far-field), diameter,
   setback, and whether it is a fixed-geometry type.
   */
  GeometryChoice geometryChoiceForDrf( const DetectorPeakResponse &drf );


  /** Re-fits `drf`'s efficiency equation to `points` and installs the result: the equation, its
   coefficient covariance (the authoritative uncertainty of an equation curve), the block covariance
   the points imply (kept as the fallback for other representations), the points themselves, and the
   energy range they span.

   This is the ONLY way the measured points and the curve are allowed to change together outside the
   Create-DRF tool - editing the points and re-interpolating them as the curve would silently swap an
   intrinsic curve for absolute measurements, and editing them without re-fitting would leave the
   refit source describing a curve nobody would get back.

   Everything else about the DRF is left exactly as it was: geometry, diameter, setback, geometry
   type, flags, FWHM, DRF source, timestamps, and any attached CeeLo response (which the caller is
   responsible for regenerating or detaching - it answers every query while attached).

   @param nterms Number of equation terms; <= 0 keeps however many the current curve has.
   @param warnings Receives the fit's non-fatal warnings, if any.
   Throws std::runtime_error (leaving `drf` untouched) when the points cannot be fit.
   */
  void refitEfficiencyFromPoints( DetectorPeakResponse &drf,
                                  const MeasuredDrfPoints &points,
                                  const int nterms,
                                  std::string &warnings );
}//namespace MakeDrfCalc

#endif //MakeDrfCalc_h
