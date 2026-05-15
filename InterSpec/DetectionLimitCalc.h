#ifndef DetectionLimitCalc_h
#define DetectionLimitCalc_h
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
#include <string>
#include <memory>
#include <vector>
#include <ostream>
#include <utility>

#include "InterSpec/PeakDef.h"

// Forward declarations
struct Material;
struct DetectorPeakResponse;

namespace SpecUtils
{
class Measurement;
class SpecFile;
}

namespace SandiaDecay
{
  class Nuclide;
}

namespace SpecUtils
{
class Measurement;
}//namespace SpecUtils


/** This header/src file holds the non-gui aware code for calculating minimum detectable activity (MDA), and detection confidence
 intervals.
 
 Notes on naming and methodology:
 - All calculations labeled as "currie" actually follow ISO 11929 methodologies for the simple gross-counts style calculation.
  E.g., A peak region is defined, and then a given number of channels are used on either side of the peak region to estimate the
  number of continuum counts inside the peak region, and variance, and then this is used to set limit from the total number or counts
  in the peak region.
  - Only single peak calculations, with no interferences are implemented.
 - All calculations labeled as "decon" use a more sophisticated  "de-convoluted" method that takes into account the shape of the peak
  and better takes into account all information provided, as well as using multiple peaks of an isotope to derive limits.  This methodology
  seems to follow the intent of Annex B of ISO 11929-3:2019, but instead these calculations form a large chi2/likelihood calculation to
  co-compute everything and hopefully do a better job.
 - Note that "Currie" is the detection-limit name, while "Curie" is the activity unit.
  
 Note: As of 20210724, these calculations have only had cursory checks performed, and have not been verified and validated to a level
  appropriate to use them for anything of importance.
 
 References consulted in developing calculations
 - ISO 11929-1:2019, ISO 11929-3:2019, ISO 11929-4:2020
 - IAEA /AQ /48
  Determination and Interpretation of Characteristic Limits for Radioactivity Measurements,
  Decision Threshold, Detection Limit and Limits of the Confidence Interval
  https://www-pub.iaea.org/MTCD/Publications/PDF/AQ-48_web.pdf
 - CURRIE, L.A., Limits for qualitative detection and quantitative determination.
  Application to radiochemistry, Anal. Chem., 40(3) (1968) 586.
  https://pubs.acs.org/doi/10.1021/ac60259a007
 - P.A. Zyla et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2020, 083C01 (2020) and 2021
  https://pdg.lbl.gov
 
 General TODO items:
 - The Currie-style limit assumes 100% of the peaks counts are within the ROI - perhaps allow a correction for this?
 - The deconvolution-style limit should allow uncertainty on the FWHM and energy calibration, to better get limits
 - More test cases need implementing
 */

namespace DetectionLimitCalc
{
#if( PERFORM_DEVELOPER_CHECKS )
void batch_test();
#endif
  
/** The input to a simple "Currie" style minimum detectable activity and detection confidence interval calculation.
 
 */
struct CurrieMdaInput
{
  /** The spectrum the calculations will be performed on.  */
  std::shared_ptr<const SpecUtils::Measurement> spectrum;
  
  /** The energy (in keV) of the photopeak limit is being derived for.
   
   This value doesnt enter into the calculation, other than it is the reference energy used for the continuum equation
   (see #CurrieMdaResult::continuum_eqn), and also its assumed this is the energy #detection_probability is derived from.
   
   Required to be between #roi_lower_energy and #roi_upper_energy, or else exception will be thrown.
   */
  float gamma_energy;
  
  /** The lower energy (in keV) of region to check for excess in.
   
   Must be within range of #spectrum.
   
   The actual energy used will be rounded to the nearest channel boundary; see #CurrieMdaResult::first_peak_region_channel.
   */
  float roi_lower_energy;
  
  /** The upper energy (in keV) of region to check for excess in.
   
   Must be greater than #roi_lower_energy and within range of #spectrum.
   
   The actual energy used will be rounded to the nearest channel boundary; see #CurrieMdaResult::last_peak_region_channel.
   */
  float roi_upper_energy;
  
  /** The number of channels below #roi_lower_energy to use to estimate the continuum. */
  size_t num_lower_side_channels;
  
  /** The number of channels above #roi_upper_energy to use to estimate the continuum. */
  size_t num_upper_side_channels;
  
  /** A value less than 1.0.  Currently used for k_alpha and k_beta.  Typically will be 0.95.
   
   In the future may split this up as two separate values, namely:
   - alpha: probability of false-positive decision (saying the signal is there, when it isnt)
   - beta: probability of false-negative decision (saying no signal, when there is signal there)
   */
  double detection_probability;
  
  /** The additional uncertainty to include when calculating limits.
   This would probably be from DRF and distance uncertainties.
   
   In AQ-48 this roughly corresponds to u_{rel}(w).
   */
  float additional_uncertainty;
  
  
  /** Default constructor that just zeros everything out. */
  CurrieMdaInput();
};//struct CurrieMdaInput


/** The results of simple "Currie" style (e.g., ) calculations. */
struct CurrieMdaResult
{
  CurrieMdaInput input;
  
  size_t first_lower_continuum_channel;
  size_t last_lower_continuum_channel;
  float lower_continuum_counts_sum;
  
  size_t first_upper_continuum_channel;
  size_t last_upper_continuum_channel;
  float upper_continuum_counts_sum;
  
  size_t first_peak_region_channel;
  size_t last_peak_region_channel;
  float peak_region_counts_sum;
  
  
  /** The Ax+b equation of the continuum, where x == (energy -  input.gamma_energy).
   E.g., the same form of polynomial
   */
  double continuum_eqn[2];
  
  float estimated_peak_continuum_counts;
  float estimated_peak_continuum_uncert;
  
  
  /** This is the number of counts in the peak region, _above_ the predicted continuum number of counts, at which point we will
   consider signal to be present.
   
   This quantity is also known as A*.
   See eqn 128 (pg47) of AQ-48, but note that equation is in activity, where this variable is in counts.
   Note that our derivation of this quantity is also slightly different to account for non-uniform energy binning, but in the limit of uniform
   binning, gives the same answer.
   
   Note: I believe this quantity corresponds to Currie's "critical level" ( L_c ),  the net signal level (instrument response) above which an
   observed signal may be reliably recognized as "detected".
   */
  float decision_threshold;
  
  
  /** Estimate of the true number of signal counts, where we would reliably detect a signal above the decision threshold.
   
   This quantity is also known as A#.
   See eqn 130 (pg 47) of AQ-48, but note that equation is in activity, where this variable is in counts.
   
   Note: I believe this quantity corresponds to Currie's "detection limit" (L_d) that is the “true” net signal level which may be a priori
   expected to lead to detection.
   */
  float detection_limit;
  
  /** This gives the nominal (best estimate) signal counts observed in the peak region.
   E.g., the amount excess above the estimated continuum.
   
   This value will be negative if a deficit of counts is observed.
   */
  float source_counts;
  
  /** This gives the lower limit, in terms of counts, on the true number of counts from signal.
   E.g. corresponds to the number of expected signal counts that we can be 95% certain the true signal is greater than.
   
   Note: may be less than zero.
   
   \sa CurrieMdaInput::detection_probability
   */
  float lower_limit;
  
  
  /** This gives the upper limit, in terms of counts, on the true number of counts from signal.
   E.g. corresponds to the number of expected signal counts that we can be 95% certain the true signal is less than.
   
   \sa CurrieMdaInput::detection_probability
   */
  float upper_limit;
  
  
  /** Default constructor that just zeros everything out. */
  CurrieMdaResult();
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  /** Checks that the \p test and \p expected results are the same.
   
   Code for checking the continuum equation is currently commented out (if there is something wrong with the continuum computation
   it will be caught by one of the other values anyway).
   */
  static void equal_enough( const CurrieMdaResult &test, const CurrieMdaResult &expected );
#endif
  
  //std::ostream &operator<<( std::ostream &strm )
  //{
  //  return strm;
  //}
};//struct CurrieMdaResult

/** Prints a summary of a CurrieMdaResult.
 
 @param strm The stream to print information out to.
 @param result The results to print out
 @param w The factor to convert from counts to source activity.  If zero or less, then only info about counts will be printed out.
  This value will usually be something like: `w = 1.0 / (shielding_transmission * drf_eff * gamma_br * live_time);`
 
 Primarily for development/debugging.
 */
std::ostream &print_summary( std::ostream &strm, const CurrieMdaResult &result, const float w );


/** For a given ROI lower and upper energy, returns the first to last channels (inclusive) that comprise the ROI, by rounding to nearest
 channels.
 
 Throws exception if energy ranges are outside spectrum, or invalid spectrum passed in.
*/
std::pair<size_t,size_t> round_roi_to_channels( std::shared_ptr<const SpecUtils::Measurement> spectrum,
                                    const float roi_lower_energy,
                                    const float roi_upper_energy );
  

/** Performs the simple gross-counts in regions style calculation.
 
 Will throw exception if input is invalid, or runs into any errors.
 */
CurrieMdaResult currie_mda_calc( const CurrieMdaInput &input );


/** How the continuum for peaks should be normalized.
 */
enum class DeconContinuumNorm : int
{
  /** The Gaussian (at a given activity) plus continuum are summed and compared to data (e.g., normal peak fitting). 
   
   I think this gives you what you usually want, e.g., the most likely peak with your current data, at the activity you
   are asserting.
   
   Note that with this option a large activity will cause the continuum to clearly be below the data, and above data for too small of activity;
   e.g., the continuum will help make up for the incorrectness of the Gaussian area.
   */
  Floating,
  
  /** A straight line is created from the edges past the ROI.
   
   I think this is sorta equivalent of saying the amplitude of the signal will not affect the
   area in the ROI, of your current data.
   
   With this option, then the continuum offset+linear coefficients will be fixed by the channels above and below the ROI, as
   specified by `DeconRoiInfo::num_lower_side_channels` and `DeconRoiInfo::num_upper_side_channels`.
   Currently when this is done, the statistical uncertainty of the above/below regions are not accounted for, and the continuum
   is just fixed - I need to put in more thought around properly handling this.
   */
  FixedByEdges,
  
  /** The continuum is fit, assuming a Gaussian area of zero, then the Gaussian component is added on top of this.
   
   This options is asserting there is no signal in your current data (most if you are doing an MDA).
   
   No uncertainty is used with the ROI (which is slightly kinda sorta, just a tiny bi defensible, if you are asserting there
   is no signal present in the data - the data variance will be used when evaluating the Chi2 anyway).
   */
  FixedByFullRange
};//enum class DeconContinuumNorm
  

/** Information about a single Region Of Interest (ROI) that is input to the deconvolution method of estimating peaks and chi2 for a
 given activity and distance.
 
 */
struct DeconRoiInfo
{
  /** The energy of the ROI start.
   
   Must be less than #roi_end, and less than #PeakInfo::energy of all #peak_infos.
   
   Will be rounded to nearest channel edge.
   */
  float roi_start;
  
  /** The energy of the ROI end.
 
   Must be greater than #roi_start, and greater than #PeakInfo::energy of all #peak_infos.
   
   Will be rounded to nearest channel edge.
   */
  float roi_end;
 
  /** The continuum type to use.
   
   Must be Linear or quadratic, for the moment...
   */
  PeakContinuum::OffsetType continuum_type;
  
  /** Whether to allow the continuum to float in the fit, or to fix the continuum using the peaks bordering the ROI, or
   use the whole ROI to determine the continuum with the assumption no signal is there.
   
   \sa num_lower_side_channels
   \sa num_upper_side_channels
   */
  DeconContinuumNorm cont_norm_method;
  
  /** The number of channels below #roi_lower_energy to use to estimate the continuum.
   
   Directly used for `cont_norm_method` of `DeconContinuumNorm::FixedByEdges`,
   and can not be zero for that case.  For the other continuum normalization methods, is used
   for the initial starting continuum equations parameters that will be fit for; in these cases the
   value will be clamped between 2 and 16 channels (incase z zero or garbage value is provided).
   */
  size_t num_lower_side_channels;
  
  /** The number of channels above #roi_upper_energy to use to estimate the continuum.
   
   Only used if `cont_norm_method` is `DeconContinuumNorm::FixedByEdges`.
   
   Note that for the other continuum normalizations, the initial (before fitting) continuum equation
   parameters are estimated using number of side channels, which is currently fixed at 4.
   */
  size_t num_upper_side_channels;
  
  
  /** Information about a photopeak-peak, for a given gamma line. */
  struct PeakInfo
  {
    /** Energy of the peak (gamma-line), in keV. */
    float energy;
    
    /** Full-width at half maximum, as given by the Detector Response Function, for the peak. */
    float fwhm;
    
    /** Counts from the source, per Bq, into 4 pi.
     
     If applicable, must have effects of shielding already accounted for.
     This number must not have effects of attenuation in air, or detector intrinsic efficiency accounted for; these will be applied during
     call to #decon_compute_peaks.
     
     For example, this value would be `branching_ratio * live_time * shielding_attenuation`,
     Where:
      - `branching_ratio` is number of this energy gamma, per decay of the parent nuclide.
      - `live_time` is the live time of the spectrum, and
      - `shielding_attenuation` is the fraction of gammas through the shielding without interaction,
      i.e. in range (0,1], where 1.0 is no shielding.
     */
    double counts_per_bq_into_4pi;
    
    /** Default zero constructor. */
    PeakInfo();
  };//struct PeakInfo
  
  /** There must be at least one peak information given, but there may be multiple.
   
   The peak means must all be between #roi_start and #roi_end, after those energies are rounded to the nearest channel edges.
   All returned fit peaks corresponding to this ROI will share a #PeakContinuum.
   */
  std::vector<DeconRoiInfo::PeakInfo> peak_infos;
  
  /** Default constructor that just zeros things out. */
  DeconRoiInfo();
};//struct DeconRoiInfo


struct DeconComputeInput
{
  double distance;
  double activity;
  
  /** Wether or not to include attenuation in air of #distance - #shielding_thickness */
  bool include_air_attenuation;
  
  /** The thickness of the shielding; zero if generic shielding or no shielding.  Subtracted from #distance if #include_air_attenuation
   is true in order to calculate attenuation in air.
   
   Must not be negative, inf, or NaN.
   */
  double shielding_thickness;
  
  std::shared_ptr<const DetectorPeakResponse> drf;
  std::shared_ptr<const SpecUtils::Measurement> measurement;
  
  std::vector<DeconRoiInfo> roi_info;
  
  /** Default constructor just zeros things out. */
  DeconComputeInput();
};//struct DeconComputeInput


struct DeconComputeResults
{
  DeconComputeInput input;
  
  double chi2;
  int num_degree_of_freedom;
  
  std::vector<PeakDef> fit_peaks;
};//struct DeconComputeResults

/** Computes the most consistent peaks, and their chi2, for a given input activity and distance.
 
 Throws exception if input is invalid, or error during calculation.
 */
DeconComputeResults decon_compute_peaks( const DeconComputeInput &input );

  
struct DeconActivityOrDistanceLimitResult
{
  /** The input variables to the calculations */
  bool isDistanceLimit;
  double confidenceLevel;
  double minSearchValue;
  double maxSearchValue;
  DeconComputeInput baseInput;
  
  // TODO: refactor generating this text to a separate function.
  std::string limitText;
  std::string quantityLimitStr;
  std::string bestCh2Text;
  
  double overallBestChi2;
  double overallBestQuantity;
  std::shared_ptr<const DeconComputeResults> overallBestResults;
  
  bool foundUpperCl;
  double upperLimit;
  double upperLimitChi2;
  std::shared_ptr<const DeconComputeResults> upperLimitResults;
  
  bool foundLowerCl;
  double lowerLimit;
  double lowerLimitChi2;
  std::shared_ptr<const DeconComputeResults> lowerLimitResults;
  
  bool foundUpperDisplay = false;
  double upperDisplayRange;
  bool foundLowerDisplay;
  double lowerDisplayRange;
  
  /** The Chi2s for a series of the quantity being found - for generating the displayed Chi2 chart from
   {quantity,Chi2}
   */
  std::vector<std::pair<double,double>> chi2s;
  
  DeconActivityOrDistanceLimitResult();
};//struct DeconActivityOrDistanceLimitResult
  
  
DeconActivityOrDistanceLimitResult get_activity_or_distance_limits( const float wantedCl,
                        const std::shared_ptr<const DeconComputeInput> base_input,
                        const bool is_dist_limit,
                        const double min_search_quantity,
                        const double max_search_quantity,
                        const bool useCurie );


// ---------------------------------------------------------------------------
// Dynamic MDA: Currie-style limit for a moving point source (or moving
// detector).  Source moves in a straight line past the detector at constant
// speed `v` with distance-of-closest-approach `d0`.  Counts are integrated
// over a short transit window `T_w` centered on closest approach, and the
// decision threshold is set against a "false-positives per hour" budget,
// since the system is implicitly running many windows.
// ---------------------------------------------------------------------------

/** Result of inverting the gross-counts Currie test with separate α / β. */
struct CurrieLcLd
{
  /** Decision threshold (counts above continuum). */
  double L_c;
  /** Detection limit (counts above continuum). */
  double L_d;
  /** Continuum uncertainty σ_B used in derivation. */
  double sigma_B;
};//struct CurrieLcLd

/** Solves the Currie L_c / L_d with separate k_α, k_β and (optional) relative
 uncertainty u_rel.  Iterative when u_rel > 0.

 @param B          Expected continuum counts under the peak.
 @param sigma_B    Uncertainty on B from side-channel extrapolation (already
                   including the Poisson term).
 @param k_alpha    Φ⁻¹(1 − α), where α is the per-trial false-positive prob.
 @param k_beta     Φ⁻¹(detection_probability).
 @param u_rel      Relative uncertainty from DRF / distance (0 if none).
 */
CurrieLcLd currie_lc_ld( const double B,
                         const double sigma_B,
                         const double k_alpha,
                         const double k_beta,
                         const double u_rel );


/** Inputs to the Dynamic MDA calculator.

 Exactly one of `speed_m_per_s`, `dca_m`, `activity_bq` must be left
 zero or negative — that quantity is the one to be solved for.  The
 other two are inputs.

 TODO:
 - move physical quantities to all use PhysicalUnits units
 */
struct DynamicMdaInput
{
  /** Background spectrum used to estimate B and σ_B in the ROI. */
  std::shared_ptr<const SpecUtils::Measurement> background_spectrum;

  /** Detector response (may be a fixed-geometry DRF, in which case the
   integral over r(t) degenerates to intrinsicEff · T_w). */
  std::shared_ptr<const DetectorPeakResponse> drf;

  /** Source nuclide (may be null; in that case branch_ratio is used as-is). */
  const SandiaDecay::Nuclide *nuclide;
  /** Age of the source, in seconds (used for branching ratio computation
   when nuclide is non-null and branch_ratio is non-positive). */
  double age;

  /** Photopeak energy, in keV. */
  float gamma_energy;

  /** Number of gammas emitted per source decay at `gamma_energy`.
   If non-positive and nuclide is non-null, this is computed from the
   nuclide / age at call time. */
  double branch_ratio;

  /** ROI lower/upper energies (keV). */
  float roi_lower_energy, roi_upper_energy;

  /** Number of channels above/below ROI for side-channel continuum.  Set
   both to zero to use the ROI itself as the continuum estimate. */
  size_t num_lower_side_channels, num_upper_side_channels;

  /** Probability of correctly detecting a real source at the limit (β). */
  double detection_probability;

  /** False-positive budget across the scenario, expressed as expected
   per-hour count (e.g., 1.0 means "at most one false alarm per hour"). */
  double max_false_positives_per_hour;

  /** Additional relative uncertainty (u_rel) on DRF / distance — folded
   into L_d iteratively. */
  double additional_uncertainty;

  /** Source speed (m/s).  Zero/negative ⇒ unknown to solve for. */
  double speed_m_per_s;
  /** Distance of closest approach (m).  Zero/negative ⇒ unknown to solve for. */
  double dca_m;
  /** Source activity (Bq).  Zero/negative ⇒ unknown to solve for. */
  double activity_bq;

  /** Window length (s).  If <= 0, auto = 2.784 · dca / speed. */
  double window_length_s;

  /** If true and `sample_period_s > 0`, T_w is rounded to a whole number of
   samples. */
  bool snap_window_to_samples;

  /** Median per-sample real time of the passthrough foreground; 0 if no
   passthrough loaded. */
  double sample_period_s;

  /** Stride between consecutive windows (s).  Required > 0.  Caller should
   pre-convert "samples per step" to seconds. */
  double window_stride_s;

  /** Optional shielding directly at the source.  When `shielding_generic` is
   true, the AN/AD pair defines the shielding; otherwise `shielding_material`
   (non-null) and `shielding_thickness` (PhysicalUnits internal length units)
   are used.  Transmission is treated as a constant multiplier at
   `gamma_energy`, independent of `r(t)`. */
  bool shielding_generic;
  double shielding_atomic_number;     // 1..100
  double shielding_areal_density;     // PhysicalUnits internal (mass per area)
  std::shared_ptr<const Material> shielding_material;
  double shielding_thickness;         // PhysicalUnits internal length units

  DynamicMdaInput();
};//struct DynamicMdaInput


/** Result of `compute_dynamic_mda`. */
struct DynamicMdaResult
{
  DynamicMdaInput input;

  /** All three filled in: the two passed in, and the one that was solved for. */
  double speed_m_per_s;
  double dca_m;
  double activity_bq;

  /** Final T_w (after the auto-or-snap step). */
  double window_length_s;
  double window_stride_s;

  /** N_per_hour = 3600 / T_step. */
  double trials_per_hour;

  /** α = F_per_hour / N_per_hour. */
  double per_trial_alpha;
  /** Φ⁻¹(1 − α). */
  double k_alpha;
  /** Φ⁻¹(detection_probability). */
  double k_beta;

  /** Expected continuum counts under the peak in T_w. */
  double background_counts_in_window;
  /** ∫ DRF::eff(E, r(t)) · exp(-μ_air · r(t)) dt over T_w, times branch ratio. */
  double signal_counts_per_bq;

  /** L_c / L_d / σ_B in counts. */
  CurrieLcLd thresholds;

  /** Diagnostic sweep for plotting (quantity, A·S − L_d).  May be empty
   if the unknown is `activity_bq`. */
  std::vector<std::pair<double,double>> sweep;

  /** Non-empty if a bracket / bisection issue prevented a solution. */
  std::string warning;

  DynamicMdaResult();
};//struct DynamicMdaResult


/** Compute the Dynamic MDA.  Throws on invalid inputs. */
DynamicMdaResult compute_dynamic_mda( const DynamicMdaInput &input );


/** One window-level "hit" from the search-mode sliding window. */
struct DynamicSearchHit
{
  int first_sample, last_sample;       // inclusive
  double window_start_s;               // seconds since first measurement
  double window_real_time_s;           // sum of real_time over window
  double counts_in_roi;
  double expected_background;
  double decision_threshold;           // L_c scaled to this window
  double excess_counts;                // counts_in_roi - expected_background
  double n_sigma_excess;               // excess_counts / sqrt(expected_background + sigmaB^2)
};//struct DynamicSearchHit

struct DynamicSearchResult
{
  std::vector<DynamicSearchHit> hits;
  size_t total_windows_evaluated;
  std::string warning;

  DynamicSearchResult();
};//struct DynamicSearchResult

/** Slide a window over the passthrough/foreground spectrum and flag windows
 whose ROI excess exceeds the decision threshold derived from
 `mda_result.thresholds.L_c` (scaled by window real time).

 @param file               The spec file to search.
 @param detectors          Detectors to sum (empty ⇒ all gamma detectors).
 @param background_samples Samples to sum into the background rate; if empty
                           and `mda_result.input.background_spectrum` is set,
                           that measurement is used instead.
 @param mda_result         Result of `compute_dynamic_mda`, which carries
                           ROI / window / threshold / DRF info.
 */
DynamicSearchResult search_passthrough_for_source(
    const std::shared_ptr<const SpecUtils::SpecFile> &file,
    const std::vector<std::string> &detectors,
    const std::set<int> &background_samples,
    const DynamicMdaResult &mda_result );

}//namespace DetectionLimitCalc


#endif //DetectionLimitCalc_h
