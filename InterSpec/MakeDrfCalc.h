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
#include "InterSpec/DetectorEfficiency.h"
#include "InterSpec/DetectorPeakResponse.h"

class PeakDef;
struct Material;
class PeakFitDetPrefs;
class MeasuredDrfPoints;
namespace SandiaDecay{ struct Nuclide; }
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

    /** The diameter the efficiency curve is quoted against, PhysicalUnits - the DRF's stated
     detector diameter: the crystal's for a plain flat disk, the whole object's
     (2*transverse_half_extent) for Make Detector Response's geometry mode.  Used by the flat-disk
     model, as the far-field footing of the kernel factor (see #intrinsicFitPoints), and recorded on
     the DRF either way. */
    double diameter = 0.0;

    /** Distance from the detector face to the crystal, PhysicalUnits; flat-disk model only (with a
     full geometry it follows from the layers). */
    double setback = 0.0;

    /** True for the fixed-geometry (e.g. activity per area) efficiency types: the measured points
     already are what the equation describes, no source geometry is divided out. */
    bool fixedGeometry = false;
  };//struct GeometryChoice


  /** A gamma line counts toward the nearest fit peak in its spectrum when within this many of the
   peak's sigma - the rule the Activity/Shielding fit uses. */
  constexpr double sm_photopeakClusterSigma = 1.25;

  /** A source's decay during the measurement is corrected for when the measurement's real time
   exceeds this fraction of its half-life - the Activity/Shielding fit's default rule. */
  constexpr double sm_decayDuringMeasFraction = 0.005;


  /** A calibration source in one source spectrum, as the user described it. */
  struct CalSource
  {
    /** Null for a source given only by its activity: each of its peaks then carries the branching
     ratio the user entered (#SamplePeak::userBr). */
    const SandiaDecay::Nuclide *nuclide = nullptr;

    /** Activity at the start of the measurement (PhysicalUnits), and its fractional 1-sigma
     uncertainty. */
    double activity = 0.0;
    double fracActivityUncert = 0.0;

    /** Age at the start of the measurement, PhysicalUnits; only used with a nuclide. */
    double age = 0.0;

    /** Source to detector face, and its 1-sigma uncertainty, PhysicalUnits; air attenuation uses
     the distance.  Negative for the fixed-geometry efficiency types. */
    double distance = -1.0;
    double distanceUncert = 0.0;

    /** Shielding around the source: a named material and its thickness (PhysicalUnits) or, with no
     material, an atomic number and areal density (PhysicalUnits); an atomic number below 1 or a
     zero areal density is no shielding. */
    std::shared_ptr<const Material> shieldMaterial;
    double shieldThickness = 0.0;
    double shieldAtomicNumber = 0.0;
    double shieldArealDensity = 0.0;

    /** False when the inputs have a problem.  Then the source's own peaks are left out, and so are
     peaks of the other sources its lines would add more than a trace to (that can't be corrected for
     without knowing its activity). */
    bool usable = true;
  };//struct CalSource


  /** A fit peak in a source spectrum. */
  struct SamplePeak
  {
    std::shared_ptr<const PeakDef> peak;

    /** Index into #SampleInput::sources of the source this peak's efficiency is for; negative when
     the peak is not used for the efficiency.  Unused peaks still matter: each gamma line counts
     toward the nearest of ALL the spectrum's peaks. */
    int source = -1;

    /** For a peak of a source without a nuclide: the branching ratio the user entered; negative if
     none (or it isn't valid). */
    double userBr = -1.0;
  };//struct SamplePeak


  /** A fit peak in a background spectrum. */
  struct BackgroundPeak
  {
    std::shared_ptr<const PeakDef> peak;

    /** Live time of the background spectrum the peak is in, seconds. */
    double liveTime = 0.0;

    /** Identifies that background spectrum: its live time counts once however many of its peaks
     match a source peak (e.g., a doublet). */
    size_t spectrum = 0;
  };//struct BackgroundPeak


  /** A source spectrum (one or more samples of a file, summed), with its peaks and sources. */
  struct SampleInput
  {
    std::vector<CalSource> sources;
    std::vector<SamplePeak> peaks;

    /** Live and real time, seconds; decay during the measurement is computed over the real time. */
    double liveTime = 0.0;
    double realTime = 0.0;
  };//struct SampleInput


  /** What one of the peaks used for the efficiency gave. */
  struct PeakPoint
  {
    enum class Status
    {
      /** #measured holds the efficiency point. */
      Used,

      /** The peak's source has no nuclide, and no branching ratio was entered for the peak. */
      NoUserBranchRatio,

      /** An unusable source (#interferer) would add more than a trace of gammas to the peak. */
      NeedsUnusableSource,

      /** The sources put no gammas into the peak (e.g., a short-lived source decayed away by the
       dates entered, or no line of the nuclide is near the peak). */
      NoSourceGammas,

      /** Nothing is left of the peak after subtracting the background. */
      NoNetArea
    };//enum class Status

    /** Index of the peak in #SampleInput::peaks. */
    size_t peak = 0;

    Status status = Status::Used;

    /** For NeedsUnusableSource: index of the unusable source in #SampleInput::sources. */
    int interferer = -1;

    /** keV: the energy of the gamma assigned to the peak, else (no source assigned) its mean. */
    double energy = 0.0;

    /** Gammas per second from all the spectrum's sources that land in the peak, after their
     shielding and the air, but before the source geometry - and its fractional 1-sigma uncertainty
     from the source activities. */
    double sourceRate = 0.0;
    double fracSourceUncert = 0.0;

    /** Peak area after any background subtraction, and its 1-sigma uncertainty. */
    double peakArea = 0.0;
    double peakAreaUncert = 0.0;

    /** Summed area of the matching background peaks, and the summed live time of the background
     spectra they are in; zero when there is no background to subtract. */
    double bkgArea = 0.0;
    double bkgLiveTime = 0.0;

    /** Whether the FWHM fit uses the peak: not a data-defined peak, which has no width, nor 511 keV
     annihilation photons, which are Doppler broadened whatever made the positron. */
    bool useForFwhm = false;

    /** For a Used peak, the absolute-efficiency point with its provenance, but for the caller's
     identity strings (sourceKey, fileName, sampleNumbers). */
    MeasuredEffPoint measured;
  };//struct PeakPoint


  /** The efficiency the peaks of one source spectrum give. */
  struct SampleResult
  {
    /** One entry per peak with a source, other than those of unusable sources, in peak order. */
    std::vector<PeakPoint> points;

    /** Sources found unusable while computing (index into #SampleInput::sources, and why); they act
     as if marked unusable. */
    std::vector<std::pair<size_t,std::string>> sourceErrors;
  };//struct SampleResult


  /** The absolute efficiency each of a source spectrum's peaks gives: its net area over the gammas
   the spectrum's sources put into it during the measurement.

   - Each gamma line of a source counts toward the nearest of the spectrum's peaks, if within
     #sm_photopeakClusterSigma of its sigma (as the Activity/Shielding fit does), so a line between
     two peaks counts once, and every source's lines in a peak are counted.
   - The decay during the measurement is corrected for when the real time exceeds
     #sm_decayDuringMeasFraction of the source's half-life.
   - The source's shielding, and the air (if `airAttenuation`) over its distance, are divided out.
   - Background peaks within 1.5 sigma (or assigned the same gamma) are subtracted, as the rate of
     their summed area over the summed live time of the spectra they are in.
   */
  SampleResult samplePoints( const SampleInput &sample,
                             const std::vector<BackgroundPeak> &background,
                             const bool airAttenuation );


  /** The number of efficiency or FWHM equation terms offered by default for `npoints` points: two
   fewer than the points (one, with three points), so chi2/dof says something about the curve; at
   least min(npoints,2), and at most `maxTerms`.  Zero with no points. */
  int defaultEquationTerms( const int npoints, const int maxTerms );


  /** The far-field intrinsic efficiency each absolute measured point corresponds to, ready to be
   fit: absolute efficiency divided by the geometric factor of its own source distance - the
   flat-disk solid angle fraction, or with a full geometry the kernel factor of
   CeeLoUtils::GeometryKernel::intrinsicFactor quoted against `geom.diameter` (so points taken at
   different distances land on one curve, near-field effects included, that the DRF's flat-disk
   model reads back in the far field).  A point's distance uncertainty becomes a fractional
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
