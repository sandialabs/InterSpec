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

#include <atomic>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <ostream>
#include <utility>

#include "InterSpec/PeakDef.h"

// Forward declarations
class DetectorPeakResponse;

namespace SpecUtils
{
class Measurement;
}

namespace SpecUtils
{
class Measurement;
}//namespace SpecUtils


/** This header/src file holds the non-gui aware code for calculating minimum detectable activity (MDA), and detection confidence
 intervals.
 
 Notes on naming and methodology:
 - All calculations labeled as "currie" follow the ISO 11929-1:2019 procedure for the simple
  gross-counts style calculation.  E.g., A peak region is defined, and then a given number of channels are used on either side of the
  peak region to estimate the number of continuum counts inside the peak region, and variance, and then this is used to set limit from
  the total number or counts in the peak region.
  - Only single peak calculations, with no interferences are implemented.

  What conforms, with the ISO 11929-1:2019 formula each implements:
   - The variance model, including the calibration term y^2*u_rel^2(w)                    Formula (25)
   - The uncertainty as a function of an assumed true value, u~(y~)                       Formulas (29)/(31)
   - Decision threshold, y* = k_{1-alpha} * u~(0)                                          Formula (32)
   - Detection limit, the smallest solution of y# = y* + k_{1-beta} * u~(y#)               Formula (34)
   - Its existence condition, k_{1-beta} * u_rel(w) < 1                                    Formula (35)
   - Limits of the probabilistically symmetric coverage interval, which apply the prior
     knowledge that the measurand is non-negative by truncating the Gaussian at zero and
     renormalizing (`omega`)                                                               Formulas (38)-(40)
   - The best estimate y^ and its associated standard uncertainty u(y^) - the mean and standard
     deviation of that same truncated, renormalized Gaussian - as
     `CurrieMdaResult::best_estimate` and `::best_estimate_uncert`.  Formula (46) (that y^ = y and
     u(y^) = u(y) are sufficient above 4*u(y)) is reproduced by the exact form rather than branched
     on, so there is no discontinuity at that threshold.                                   Formulas (44)-(46)
   - y^ - not the primary result - is what the GUI, the batch reports and the LLM JSON quote as the
     measured value, as clause 10 NOTE 1 has it.  The primary result is retained beside it (it is
     what the decision rule uses) rather than omitted.                                     Clause 10, NOTE 1
   - alpha, beta, and the coverage probability kept as three separate quantities            Clauses 8.1, 9

  What does NOT conform, and must not be claimed:
   - Only the probabilistically symmetric coverage interval is offered; the shortest coverage
     interval is not.  ISO leaves the choice to the regulator, and the two differ by ~10% at low
     counts.                                                                               Formulas (42)/(43)
   - `CurrieMdaInput::detection_probability` is a ONE-SIDED confidence level, so selecting "95%"
     gives ISO's gamma = 0.10 interval, not gamma = 0.05.  To reproduce an ISO two-sided 95%
     interval, pass 0.975.
   - The model is the ROI specialization of Y = (X1 - X2)*W; the further input quantities X3 and X4
     (blank, shielding factor, ...) are not carried separately - all non-counting uncertainty is
     collapsed into #CurrieMdaInput::additional_uncertainty as u_rel(w).
   - This is the ISO 11929-1 (Gaussian) route, not the ISO 11929-2 (Bayesian/Monte Carlo) one.  The
     two differ at very low counts; for the clause 7 example ISO publishes a# of 3,90E-2 by Part 1
     against 3,58E-2 by Part 2.
   - The rounding and documentation stipulations of clauses 5.3 and 11 are not implemented.
 - All calculations labeled as "decon" use a more sophisticated  "de-convoluted" method that takes into account the shape of the peak
  and better takes into account all information provided, as well as using multiple peaks of an isotope to derive limits.  This methodology
  seems to follow the intent of Annex B of ISO 11929-3:2019, but instead these calculations form a large chi2/likelihood calculation to
  co-compute everything and hopefully do a better job.
  - It finds its limit where the profiled statistic rises by a threshold, and only ever evaluates
   non-negative activities, so it does not use the `omega` construction above and makes no ISO
   conformance claim.  Its measured coverage runs conservative - see `DeconCoverageStudy`.
 - Note that "Currie" is the detection-limit name, while "Curie" is the activity unit.

 Verification status:
 - The Currie path is checked against published worked examples: ISO 11929-4:2020 clause 6
  (`Iso11929Part4Example1`) and clause 7 (`Iso11929Part4Example2LowCounts`), and IAEA AQ-48 Table 16
  (`Table16OfAQ48`).  Every published value is reproduced to within the resolution ISO prints it at;
  the tests emit the measured agreement so it is recorded rather than hidden behind a tolerance.
  `CurrieCoverageIntervalIsCalibrated` additionally verifies, by numerical integration rather than by
  the analytic form, that the coverage interval covers what it claims, and
  `CurrieBestEstimateMatchesTruncatedGaussianMean` does the same for the best estimate: it integrates
  the truncated Gaussian's mean and variance directly, so the check cannot agree with Formulas
  (44)/(45) by sharing a mistake.  `CurrieBestEstimateConvergesToPrimaryResult` sweeps y/u(y) and
  checks the relations ISO states, including that Formula (46) is reached without being branched on.
  - Note AQ-48 Table 16's printed a# = 0.471705 omits the u_rel(w) = 0.0478 that the same table lists
   as an input; the error is the document's, and is not asserted against.
 - The decon path has NOT been checked against any published example.  As of 20210724 these
  calculations had only had cursory checks performed; that remains true for everything outside the
  Currie path above, which has not been verified and validated to a level appropriate to use it for
  anything of importance.

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
  
  /** The confidence level of the reported interval; e.g., 0.95.  Must be greater than 0.05 and less
   than 1.0.

   This is the coverage of #CurrieMdaResult::lower_limit / #CurrieMdaResult::upper_limit - the
   "less than X at Y CL" the tools report.

   It is also the fallback for the two decision rates: when #alpha or #beta is left unspecified, that
   quantile is taken from this value, which reproduces exactly the historical behavior of using one
   number for all three roles.
   \sa alpha, beta
   */
  double detection_probability;

  /** Probability of deciding a signal is present when there is none - the false-positive rate that
   sets #CurrieMdaResult::decision_threshold,  L_c = k_{1-alpha} * sigma_0.

   Must be less than 0.5, and at least 1.0E-9 when given.  Zero or negative means "not specified", in
   which case `1 - detection_probability` is used - so a default-constructed #CurrieMdaInput, and
   every caller written before this field existed, behaves exactly as before.
   \sa currie_k_alpha_beta
   */
  double alpha;

  /** Probability of failing to detect a signal whose true size is the detection limit - the
   false-negative rate that sets #CurrieMdaResult::detection_limit, L_d.

   Must be less than 0.5, and at least 1.0E-9 when given.  Zero or negative means "not specified", in
   which case `1 - detection_probability` is used.
   \sa currie_k_alpha_beta
   */
  double beta;

  /** The relative systematic uncertainty to include when calculating limits: everything that scales
   the expected counts and is not counting statistics - detector efficiency, gamma branching ratio,
   and the source-to-detector distance (which enters squared, a point sources efficiency going as
   1/r^2).

   In AQ-48 this roughly corresponds to u_{rel}(w).  Must be in [0, 1).
   \sa combine_systematic_uncertainty
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
  
  /** This gives the nominal signal counts observed in the peak region.
   E.g., the amount excess above the estimated continuum.

   This value will be negative if a deficit of counts is observed.

   This is ISO 11929-1:2019's *primary result*, y - the unconstrained net signal, and is what the
   decision rule compares against #decision_threshold.  It is not what ISO says to quote as the
   measured value; that is #best_estimate.
   */
  float source_counts;

  /** The best estimate of the signal counts, ISO 11929-1:2019 Formula (44) - the mean of the
   Gaussian after truncation at zero and renormalization.

   This, not #source_counts, is what ISO says to quote as the measured value.  The two converge as
   the signal strengthens - they agree to about 3 parts in 1E5 at #source_counts equal to four times
   its uncertainty, which is where ISO's Formula (46) says treating them as equal is sufficient, and
   to better than a part in 1E9 by seven times - and diverge as it weakens.  #source_counts remains
   the primary result and is what the decision rule compares against #decision_threshold.
   */
  float best_estimate;

  /** Standard uncertainty of #best_estimate, ISO 11929-1:2019 Formula (45).

   Never larger than the uncertainty of the primary result: truncation removes the part of the
   distribution that sat on impossible negative activities.
   */
  float best_estimate_uncert;

  /** This gives the lower limit, in terms of counts, on the true number of counts from signal.
   E.g. corresponds to the number of expected signal counts that we can be 95% certain the true signal is greater than.

   Non-negative: the coverage interval applies the prior knowledge that the measurand cannot be
   negative, per ISO 11929-1:2019 clause 9, by truncating the Gaussian at zero and renormalizing.
   (Before that treatment was added this could come out below zero, which is not a possible activity.
   #source_counts is the unconstrained estimate and may still be negative.)

   \sa CurrieMdaInput::detection_probability, source_counts
   */
  float lower_limit;


  /** This gives the upper limit, in terms of counts, on the true number of counts from signal.
   E.g. corresponds to the number of expected signal counts that we can be 95% certain the true signal is less than.

   Also from the zero-truncated, renormalized Gaussian; see #lower_limit.  This makes the limit
   *larger* than the plain `source_counts + k*sigma` - by about 19% at 95% and a factor of two at the
   68.27% setting when no signal is present - because the probability that would have sat on negative
   activities has to go somewhere, and it goes onto the positive side.

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


/** The normal quantiles #currie_mda_calc will actually use, with the #CurrieMdaInput::alpha and
 #CurrieMdaInput::beta "not specified" sentinels resolved against
 #CurrieMdaInput::detection_probability.

 This is the one place the sentinel is interpreted, so that a GUI showing the effective error rates
 and the calculation using them cannot drift apart.

 @returns The pair {k_alpha, k_beta}.

 Throws exception if any of the three probabilities is out of range; these are exactly the checks
 #currie_mda_calc makes on them.
 */
std::pair<double,double> currie_k_alpha_beta( const CurrieMdaInput &input );


/** Combines a distance uncertainty and a relative efficiency uncertainty into the single relative
 systematic uncertainty #CurrieMdaInput::additional_uncertainty takes.

 A point sources detection efficiency goes as 1/r^2, so a relative distance uncertainty reaches the
 expected counts with a factor of two:

     u_rel = sqrt( (2*sigma_d/d)^2 + u_eff^2 )

 @param distance Source to detector distance.  Only consulted when \p distance_uncertainty is
        positive, and must then be positive.
 @param distance_uncertainty The 1-sigma uncertainty of \p distance, in the same units.  Zero or
        negative means none.
 @param efficiency_rel_uncertainty The 1-sigma _relative_ uncertainty of the counts expected per unit
        activity - detector efficiency and gamma branching ratio together, which enter identically.
        Zero or negative means none.
 @param inverse_square Whether the efficiency actually follows 1/r^2.  Pass false for a fixed-geometry
        detector response (\sa DetectorPeakResponse::isFixedGeometry), where the distance does not
        enter the efficiency at all, and the distance term must be dropped rather than scaled.
 @returns The combined relative uncertainty; always non-negative.  May be greater than or equal to
        one, which #currie_mda_calc rejects: that is a user-input problem, so the caller decides how
        to report it rather than a clamp silently changing the number the user entered.

 Throws exception if an argument is not finite, or if \p distance_uncertainty is positive while
 \p distance is not.
 */
double combine_systematic_uncertainty( const double distance,
                                       const double distance_uncertainty,
                                       const double efficiency_rel_uncertainty,
                                       const bool inverse_square );


/** Options for the Currie-style detection limit check made for a single peak.

 \sa currie_check_for_peak
 */
struct PeakCurrieCheckOptions
{
  /** Confidence level of the limits; e.g., 0.95 for 95%.  Must be greater than 0.5 and less than 1. */
  double confidence_level = 0.95;

  /** Number of channels on each side of the peak region used to estimate the continuum.

   A value of zero puts `currie_mda_calc(...)` into its "the spectrum is asserted to be background"
   mode, which is a different calculation than this check performs, so zero is treated as one.
   */
  size_t num_side_channels = 4;

  /** The width of the peak region, in multiples of the peak FWHM.
   The default of 2.5 (i.e., mean +-1.25 FWHM) is what ISO 11929:2010 recommends.
   */
  double roi_num_fwhm = 2.5;
};//struct PeakCurrieCheckOptions


/** The outcome of a Currie-style detection limit check made for a single peak.

 This is the "is there a detectable signal here, and how much could be hiding" question, asked of a
 peak region of the spectrum.  It is meaningful both for a peak that could not be fit (where the
 limit is the answer), and for a peak that was fit (where it is a quality check - a peak can be fit,
 yet still sit below the level at which a signal can be reliably claimed).

 \sa currie_check_for_peak
 */
struct PeakCurrieCheck
{
  /** How the observed counts in the peak region compare to what the continuum predicts. */
  enum class ResultType : int
  {
    /** Signal is below the decision threshold (L_c); only an upper limit can be quoted. */
    NotDetected,

    /** Signal is above the decision threshold (L_c). */
    Detected,

    /** Fewer counts observed than the neighboring continuum predicts; even the upper limit is
     less than zero counts.  Usually means the continuum estimate is poor.
     */
    Deficit,

    /** The check could not be performed; see `error_message`. */
    Error
  };//enum class ResultType

  /** Whether `result` holds a valid calculation; when false see `error_message`. */
  bool computed = false;

  /** Why the check could not be made; only non-empty when `computed` is false. */
  std::string error_message;

  /** The underlying calculation.  When `computed` is false, only `result.input` is meaningful, and
   only to report the options that were attempted.
   */
  CurrieMdaResult result;

  ResultType result_type = ResultType::Error;

  /** True when the peak region holds essentially no counts, so that the Gaussian statistics the
   calculation uses have broken down (the decision threshold and upper limit both collapse to
   zero, which would make a single stray count read as a detection).
   */
  bool region_is_empty = false;

  /** The fraction of a Gaussian peak that falls inside the peak region; see
   `gaussian_fraction_in_roi(...)`.  The limit is on the counts within the region, and is not
   corrected for the signal outside it, so a value meaningfully below 1 makes the limit optimistic.
   */
  double signal_fraction_in_roi = 1.0;

  /** A table-cell sized summary; e.g., "Less than Lc". */
  std::string short_description;

  /** A sentence or two describing the outcome, with the numbers in it. */
  std::string result_summary;
};//struct PeakCurrieCheck


/** Returns a short, stable string for the result type; e.g. "NotDetected". */
const char *to_str( const PeakCurrieCheck::ResultType type );


/** Performs a Currie-style detection limit check for a single peak.

 The peak region is `peak.mean()` +- half of `options.roi_num_fwhm` FWHM, with the continuum
 estimated from `options.num_side_channels` channels on each side.  Does not throw; failures are
 reported through `PeakCurrieCheck::computed` and `error_message`.

 @param peak The peak to evaluate; its mean and width define the region.  For "data defined" peaks,
        see `peak_width_for_currie_check(...)` for how the width is determined.
 @param spectrum The spectrum to evaluate the peak region in.
 @param options See `PeakCurrieCheckOptions`.
 @param peak_was_fit Whether this peak was actually fit in `spectrum`.  Only affects the wording of
        `short_description` and `result_summary` - the numbers are the same either way.
 */
PeakCurrieCheck currie_check_for_peak( const PeakDef &peak,
                              const std::shared_ptr<const SpecUtils::Measurement> &spectrum,
                              const PeakCurrieCheckOptions &options,
                              const bool peak_was_fit );


/** Returns the FWHM to use for a peaks detection limit check; for peaks that arent Gaussian (i.e.,
 "data defined" peaks), the ROI is assumed to be the usual 2.5 FWHM wide, so that the default
 `PeakCurrieCheckOptions::roi_num_fwhm` reproduces the peaks own ROI width.

 Returns zero if a width could not be determined.
 */
double peak_width_for_currie_check( const PeakDef &peak );


/** Returns the fraction of a Gaussian peak that lies within +-(num_fwhm/2) FWHM of its mean.

 E.g., 2.5 gives 0.9968, and 1.0 gives 0.7610.  Returns zero for non-positive input.
 */
double gaussian_fraction_in_roi( const double num_fwhm );


/** Formats a confidence level for display; e.g., 0.95 becomes "95%", and 0.9999 becomes "1-1E-04".

 The complement form for very high levels is a published batch-report contract
 (`MdaConfidenceLevelStr` / `ConfidenceLevelStr`); use `confidence_level_pct_str` for anything
 user-facing that has to line up with the confidence-level selectors.
 */
std::string confidence_level_str( const double confidence_level );


/** Formats a confidence level as a plain percentage, with enough decimals that the complement keeps
 two significant figures; e.g. 0.95 becomes "95%" and 0.999999426696856 becomes "99.999943%".

 Exists because the 4-sigma and 5-sigma selections differ from 100% only in the 5th and 7th decimal
 place, so any fixed precision renders both as "100%" and makes two different confidence levels
 indistinguishable - which is what the deconvolution result text used to do.  Returns "?" for a
 value that is not strictly between 0 and 1, including NaN.

 \sa confidence_level_str, which keeps the "1-5.7E-07" complement form the batch reports document.
 */
std::string confidence_level_pct_str( const double confidence_level );


/** Returns the energy range a Currie-style limit used, including the side channels used to estimate
 the continuum.  Useful for detecting that a neighboring peak may contaminate the continuum.
 */
std::pair<double,double> currie_check_energy_range( const CurrieMdaResult &result );


/** The Poisson deviance ("Cash" statistic) contribution of a single channel:

 \code
   2*( expected - observed + observed*log(observed/expected) )
 \endcode

 with the logarithmic term defined to be zero when \p observed is zero.  Summed over independent
 channels this is the likelihood-ratio statistic of Cash (1979), and differences of it between
 nested models follow the same asymptotic chi-square law that a Gaussian chi-square difference
 would - but without the modified-Neyman variance bias, and without needing a variance floor for
 empty channels.

 \p expected must be strictly positive; \p observed must be non-negative.  Throws if either is
 NaN/infinite, or if \p expected is not positive.
 */
double poisson_deviance( const double observed, const double expected );


/** The smallest expected counts a channel is allowed, given the mean observed counts over the
 channels being scored.

 A *relative* floor, unlike the modified-Neyman variance floor of 1.0 count this calculation used
 to rely on: it keeps `log(observed/expected)` finite without imposing an absolute scale on spectra
 that have been scaled, background subtracted, or taken for very short dwells.

 Exported because the coverage study has to apply the identical floor - a study that floors
 differently is characterizing a different estimator than the one that ships - and because the two
 production call sites otherwise each re-derive it and can drift apart.

 \p mean_observed_counts is the mean over the scored channels; values below one are treated as one,
 so an almost-empty region does not drive the floor to zero.
 */
double min_expected_channel_counts( const double mean_observed_counts );


/** Result of fitting polynomial continuum coefficients to Poisson channel data.

 \sa fit_continuum_poisson
 */
struct PoissonContinuumFit
{
  /** The fit coefficients, defined relative to the reference energy passed in, and using the same
   channel-integrated polynomial basis as `PeakContinuum::offset_eqn_integral(...)`.
   */
  std::vector<double> coefficients;

  /** Row-major `num_coefficients*num_coefficients` covariance of #coefficients, from the inverse
   of the expected (Fisher) information.  Empty if the information matrix could not be inverted.
   */
  std::vector<double> covariance;

  /** Row-major `num_coefficients*num_coefficients` expected (Fisher) information of #coefficients,
   in the same units as #covariance.

   Provided alongside #covariance because a caller that wants to turn this fit into a Gaussian
   constraint on a later fit needs the *precision* matrix, and taking it from here avoids inverting
   a covariance that was itself produced by an inversion.

   #covariance is the Moore-Penrose pseudo-inverse of this, so the two are true mutual inverses
   only at full rank; at less than full rank the pseudo-inverse silently drops the deficient
   directions, and this matrix is the more faithful of the pair.  Empty if the information could
   not be formed - which is not the same condition under which #covariance is empty, so test each
   before using it.  \sa fit_continuum_poisson
   */
  std::vector<double> information;

  /** The Poisson deviance at #coefficients, including any constraint penalty. */
  double statistic = 0.0;

  /** Total number of optimizer iterations used, summed over all attempts. */
  size_t num_iterations = 0;

  /** How far down the fallback ladder the fit had to go; zero means plain damped Newton from the
   initial estimate succeeded.  Useful as a numerical-health diagnostic.
   */
  size_t num_restarts = 0;

  bool converged = false;

  /** Description of why the fit failed; non-empty if, and only if, #converged is false. */
  std::string error;
};//struct PoissonContinuumFit


/** One independent Poisson observation that constrains the continuum coefficients.

 Listing the channels explicitly, rather than as a contiguous run of a spectrum, is what lets a
 single optimizer serve the two cases that are not a simple energy range:

 - `DeconContinuumNorm::FixedByEdges` fits the channels *beside* an ROI, which are two blocks with
   the ROI as a gap between them.
 - `DeconMeasurementModel::BackgroundReference` fits one continuum to two different observations of
   the same energies - the reference spectrum and the projected sample measurement - which have
   different exposures and so different #continuum_scale.

 Channels are treated as statistically independent, so the same energy range may legitimately
 appear more than once when it was measured more than once.
 */
struct PoissonChannel
{
  /** Lower and upper energy of the channel; #upper_energy must exceed #lower_energy. */
  double lower_energy = 0.0;
  double upper_energy = 0.0;

  /** Observed counts in this channel; must be non-negative. */
  double observed = 0.0;

  /** Counts contributed by the fixed trial signal; must be non-negative. */
  double fixed_signal = 0.0;

  /** Multiplies the polynomial's contribution to this channel's expected counts, so that

   \code
     expected = continuum_scale*continuum(coefficients) + fixed_signal
   \endcode

   Use 1 for a channel of the spectrum the continuum is defined on, and
   `sample_exposure/reference_exposure` for a channel of a measurement being projected from it.
   Must be strictly positive.
   */
  double continuum_scale = 1.0;
};//struct PoissonChannel


/** Finds the polynomial continuum coefficients that minimize the Poisson deviance of

 \code
   expected_i = continuum_i(coefficients) + fixed_signal_i
 \endcode

 against \p observed, subject to every \p expected_i remaining strictly positive.

 This is the nuisance-parameter profiling step of the deconvolution detection limit: the trial
 signal is held fixed and the continuum is re-solved, so that the statistic reported afterwards is
 evaluated at the continuum that actually minimizes it.  Note the contrast with
 `PeakFit::fit_amp_and_offset_imp(...)`, which minimizes a modified-Neyman chi-square of
 `observed - fixed_signal` after clipping that difference at zero, and so does not minimize the
 statistic subsequently reported.

 @param channel_edges Lower energies of each channel; must have at least \p nbin + 1 entries.
 @param observed The observed counts of each channel; must have at least \p nbin entries, all
        non-negative.
 @param fixed_signal Counts in each channel from the fixed trial signal; may be nullptr, which is
        equivalent to all zeros.  Must be non-negative.
 @param nbin The number of channels.
 @param num_coefficients The number of polynomial terms; 1 is constant, 2 linear, 3 quadratic.
 @param reference_energy The energy the polynomial is defined relative to.
 @param initial_coefficients Starting estimate; may be empty, in which case a weighted
        least-squares solve provides the seed.
 @param constraint_center Optional Gaussian-constraint center for the coefficients; when non-empty
        must have \p num_coefficients entries, and \p constraint_inverse_covariance must also be
        supplied.  Adds `(c - center)^T InvCov (c - center)` to the objective, which is how
        `DeconContinuumNorm::FixedByEdges` carries the statistical uncertainty of its side
        channels.  This is a full multivariate penalty: for a linear continuum it constrains the
        offset and the slope jointly, including their correlation.
 @param constraint_inverse_covariance Optional row-major
        `num_coefficients*num_coefficients` inverse covariance for \p constraint_center.

 Never throws for a merely difficult fit; instead returns a result with `converged == false` and a
 populated `error`.  Throws only for structurally invalid input (null pointers, zero channels,
 mismatched constraint sizes, non-finite data).
 */
PoissonContinuumFit fit_continuum_poisson( const float * const channel_edges,
                                          const double * const observed,
                                          const double * const fixed_signal,
                                          const size_t nbin,
                                          const size_t num_coefficients,
                                          const double reference_energy,
                                          const std::vector<double> &initial_coefficients,
                                          const std::vector<double> &constraint_center = {},
                                          const std::vector<double> &constraint_inverse_covariance = {} );


/** As above, but over an explicit list of independent Poisson observations rather than a
 contiguous run of channels.

 This is the general form; the array-based overload builds \p channels from its arguments and
 calls straight through, so the two cannot give different answers for the same problem.  Use this
 one when the channels are not a single contiguous energy range, or when they do not all share one
 exposure.  \sa PoissonChannel

 @param channels The observations; must have at least \p num_channels entries.
 @param num_channels The number of observations; must be at least \p num_coefficients.

 Throws for structurally invalid input, including a non-positive `continuum_scale` or a channel
 whose upper energy does not exceed its lower energy.
 */
PoissonContinuumFit fit_continuum_poisson( const PoissonChannel * const channels,
                                          const size_t num_channels,
                                          const size_t num_coefficients,
                                          const double reference_energy,
                                          const std::vector<double> &initial_coefficients,
                                          const std::vector<double> &constraint_center = {},
                                          const std::vector<double> &constraint_inverse_covariance = {} );


/** Which quantity the deconvolution profile scan should report.

 These are genuinely different products, and they use different thresholds on the profiled
 statistic.  Reporting one while describing the other is the single easiest way to misstate a
 detection limit, so the choice is explicit rather than inferred.
 */
enum class DeconLimitType : int
{
  /** A one-sided upper bound: "the activity is less than U at this confidence".

   Threshold `quantile( chi_squared(1), 2*CL - 1 )`; at 95% that is 2.705543.  A chi-square(1)
   variate is the square of a standard normal, so its distribution has already folded both normal
   tails together - a one-sided normal confidence level `CL` therefore corresponds to chi-square
   probability `2*CL - 1`, not `CL`.

   This is the default, and what every existing caller wants.
   */
  OneSidedUpperLimit = 0,

  /** A central two-sided interval: "the activity is between L and U at this confidence".

   Threshold `quantile( chi_squared(1), CL )`; at 95% that is 3.841459.

   Note that taking the two roots found at the *one-sided* threshold does not give this: those are
   two individually-one-sided bounds, and together they cover only about `2*CL - 1` centrally, so
   two 95% one-sided roots form roughly a 90% central interval.
   */
  CentralInterval = 1
};//enum class DeconLimitType


/** The change in the profiled statistic that defines \p limit_type at confidence level
 \p confidence_level.

 Throws if \p confidence_level is not strictly between 0.5 and 1, including for NaN.
 */
double decon_limit_delta( const double confidence_level, const DeconLimitType limit_type );


/** How the continuum under a peak is determined.

 This says *how the background beneath the peak is estimated*, and nothing about whose measurement
 the limit describes - that is `DeconMeasurementModel`, and the two are chosen independently.
 Conflating them is what produced the retired #FixedByFullRange option below.

 Every value here propagates the uncertainty of whatever data determines the continuum; an option
 that froze the continuum and propagated nothing is what #FixedByFullRange was, and it under-covered
 badly.
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
  
  /** The channels beside the ROI are background-control observations, and the continuum is fit to
   them and the ROI together.

   The side channels named by `DeconRoiInfo::num_lower_side_channels` and
   `DeconRoiInfo::num_upper_side_channels` enter the *same* Poisson likelihood as the ROI's own
   channels, carrying no trial signal.  The continuum coefficients are then profiled over the
   union, so the ROI may pull the continuum but only as far as the side counts allow, and the side
   channels' own counting uncertainty propagates into the limit.  "Anchored by the channels beside
   the ROI", not "fixed".

   Two earlier constructions were tried and are worth not repeating.  Freezing a line through the
   side channels propagates none of their uncertainty and under-covers badly.  Fitting them
   separately and carrying their covariance in as a Gaussian penalty is only a second-order
   approximation of the joint likelihood used here, and it fails where this option matters most:
   the expected information weights each channel by `1/E`, so an empty side channel driven onto
   the positivity floor contributes an enormous `1/E` and over-constrains the continuum, biasing
   limits low below roughly three counts per channel.  Summing the deviance directly has no
   information matrix, no `1/E`, and no quadratic approximation.

   Side channels that fall inside another ROI are skipped, since they are already contributing to
   the summed statistic as that ROI's data.
   */
  FixedByEdges,

  /** DEPRECATED.  Retained only so that saved states and URLs that selected it still decode.

   This asserted the spectrum contains no signal, and used it to predict a future measurement's
   sensitivity - which is a different *measurement model*, not a continuum treatment, and is now
   `DeconMeasurementModel::BackgroundReference` combined with #Floating.  As a continuum treatment
   it was circular: any real signal was absorbed into the zero-signal continuum estimate and the
   uncertainty of that estimate was then ignored, which measured 40.05% coverage where 95% was
   claimed.

   Never silently reinterpret a stored selection of this value; migrate it visibly.
   */
  FixedByFullRange
};//enum class DeconContinuumNorm


/** Which measurement the limit is a statement about.

 Distinct from `DeconContinuumNorm`, which says how the continuum under a peak is determined.  This
 says whose spectrum is being described: the one in hand, or one that has not been taken yet.
 */
enum class DeconMeasurementModel : int
{
  /** A limit on how much signal could be present in the spectrum that was loaded.

   The ordinary case, and what every caller wanted before this enum existed.
   */
  CurrentSpectrum = 0,

  /** The loaded spectrum is asserted to be signal-free background, and is used to predict the
   sensitivity of a *future* measurement of `DeconComputeInput::sample_exposure`.

   The reference counts and the projected sample are scored jointly, sharing one continuum, so the
   reference's finite counting statistics propagate into the answer; as the reference exposure
   grows the result converges to the ideal known-background case.  Because there is no future
   sample to observe, the sample is taken at the counts it expects, with no measurement noise, making the result a
   *predicted median sensitivity* rather than an observed bound.

   The result is therefore not an upper limit on signal in the loaded spectrum, and must not be
   presented as one - `DeconActivityOrDistanceLimitResult::is_predicted_sensitivity` records this.

   **How good the expected-counts step is, measured.**  Setting #observed_sample re-scores this same two-block
   likelihood with a real sample in place of the expected-counts block, which is the estimand the prediction is
   supposed to be the median of.  `PredictedSampleBiasParameterSweep` compares the two over 64 cells - four
   sample-to-reference exposure ratios `k`, four continuum levels and four peak widths, spanning
   region background counts from 0.5 to 4000.

   The expected-counts step holds up: **no trend in `k` at all** (medians 1.02 / 1.02 / 0.98 / 0.99 at
   `k` = 0.25 / 1 / 4 / 16), and every cell whose region holds more than ~5 counts lands within
   0.70-1.35.  Below that the median of a limit is not a stable quantity and the ratio scatters from
   0.21 to 1.51, which is small-number noise rather than bias.

   What does *not* hold up is the continuum model, and it is not specific to this measurement model:
   see `ContinuumMisspecificationCollapsesLimit`, where a `Linear` continuum fitted to a truth it cannot
   represent collapses the limit by up to 7x - on the ordinary `CurrentSpectrum` path just as much as
   here.  A `BackgroundReference` prediction that disagrees with a real analysis of the same region is
   far more likely to be telling you the region's continuum is not linear than that the expected-counts step
   failed.
   */
  BackgroundReference = 1
};//enum class DeconMeasurementModel


/** Combo-box index to enum, and back, for the three deconvolution selectors.

 The combo order is not the enum order for every widget - Simple MDA's continuum selector had
 index 1 meaning `FixedByFullRange` and index 2 meaning `FixedByEdges`, the reverse of the enum,
 with nothing guarding it - and a wrong index silently selects a different statistical treatment
 rather than failing.  Both widgets and both URL decoders go through these, so the order lives in
 exactly one place.

 The `*_from_index` functions throw for an index outside the offered range.
 */
DeconContinuumNorm continuum_norm_from_index( const int index );

/** Throws for `DeconContinuumNorm::FixedByFullRange`, which is deprecated and is not offered as a
 selectable choice; a stored state carrying it must be migrated, never displayed.
 */
int index_from_continuum_norm( const DeconContinuumNorm norm );

/** How many continuum treatments a user may pick from; `FixedByFullRange` is excluded. */
int num_selectable_continuum_norms();

DeconMeasurementModel measurement_model_from_index( const int index );
int index_from_measurement_model( const DeconMeasurementModel model );

DeconLimitType limit_type_from_index( const int index );
int index_from_limit_type( const DeconLimitType type );


/** Decodes a stored continuum-normalization URL token, migrating the deprecated ones.

 Both tools persist their state in the database and in QR codes, so tokens written by any past
 version have to keep decoding forever.  The vocabularies are disjoint, so this accepts all of them
 regardless of the URI's version number - a hand-edited or re-saved `VER=2` URI can still carry a
 legacy token, and migrating on the token rather than the version is what makes that safe.

 | token | result |
 |---|---|
 | `UNKNOWN` (Simple, v1) / `FLOAT` (v2) | #DeconContinuumNorm::Floating, #DeconMeasurementModel::CurrentSpectrum |
 | `FIXED` (Simple, v1) / `EDGES` (v2)   | #DeconContinuumNorm::FixedByEdges, #DeconMeasurementModel::CurrentSpectrum |
 | `NOSIG` (Simple, v1) / `FULL` (Tool, v1) | Floating + BackgroundReference, and \p migrated set true |

 @param token The token text; matched case-insensitively, by prefix.
 @param norm Set on success; untouched otherwise.
 @param model Set on success; untouched otherwise.
 @param migrated Set true only when a deprecated token was reinterpreted.  The caller MUST surface a
        visible notice when it is - the retired option measured about 40% coverage where 95% was
        claimed, so silently loading it as something else would change what a saved number means
        without telling anyone.

 @returns false for an unrecognized token, leaving all outputs untouched, so the caller can log and
        keep its default rather than guess.
 */
bool decode_continuum_norm_token( const std::string &token,
                                  DeconContinuumNorm &norm,
                                  DeconMeasurementModel &model,
                                  bool &migrated );

/** The token `encodeStateToUrl` writes for \p norm; the inverse of
 `decode_continuum_norm_token` for the two non-deprecated values.
 Throws for `DeconContinuumNorm::FixedByFullRange`, which must be migrated, never re-emitted.
 */
std::string continuum_norm_token( const DeconContinuumNorm norm );
  

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

  /** Whether the limit describes #measurement, or a future measurement predicted from it.
   \sa DeconMeasurementModel
   */
  DeconMeasurementModel measurement_model = DeconMeasurementModel::CurrentSpectrum;

  /** Live time, in seconds, of the future measurement that
   `DeconMeasurementModel::BackgroundReference` predicts.

   Ignored by `DeconMeasurementModel::CurrentSpectrum`, and superseded by #observed_sample when that
   is set.  Zero or negative means "the same exposure as #measurement", which is the natural default
   and makes the reference-to-sample ratio one.
   */
  double sample_exposure = 0.0;

  /** The sample measurement, once it has actually been taken.

   Only meaningful under `DeconMeasurementModel::BackgroundReference`, and ignored otherwise.  Null -
   the default - leaves that model predicting: the sample block of the two-block likelihood is the
   counts it expects under the null, with no measurement noise, and the answer is a *predicted* sensitivity.

   Setting it replaces that expected-counts block with the sample's real counts, so #measurement and this are
   scored as a genuine joint likelihood against one continuum.  The answer is then a bound on the
   signal in this sample, informed by the reference, rather than a prediction.  #sample_exposure is
   ignored: the sample's own live time is the exposure.

   Must have the same number of channels and the same energy calibration as #measurement, and a
   positive live time; `decon_compute_peaks` throws otherwise.  The two spectra's channels are
   paired by index.
   */
  std::shared_ptr<const SpecUtils::Measurement> observed_sample;

  /** Default constructor just zeros things out. */
  DeconComputeInput();
};//struct DeconComputeInput


struct DeconComputeResults
{
  DeconComputeInput input;

  /** The value of the test statistic named by #statistic_name.

   This is a Poisson deviance (Cash statistic), not a chi-square: see #poisson_deviance.  It is
   summed once per independent channel, over regions of interest that have been merged where they
   overlapped, so no channel contributes twice.  Only *differences* of it between trial
   activities/distances are meaningful; its absolute value has no goodness-of-fit interpretation
   the way a chi-square would.

   Specifically `E[chi2] != nbin - npar`, so #num_degree_of_freedom exists for display continuity
   and `chi2/DOF` is *not* a fit-quality measure - at low counts per channel it is not even
   approximately one for a correct model.
   */
  double chi2;

  /** Names the statistic in #chi2, for display and reports; currently always "Cash". */
  std::string statistic_name;

  /** Channels minus the nuisance parameters estimated from those same channels.

   Display-only metadata.  It does *not* select the confidence threshold: that is set by the single
   activity or distance parameter being profiled, via `decon_limit_delta`.  And because #chi2 is a
   Poisson deviance rather than a chi-square, `chi2/DOF` is not a goodness-of-fit ratio; see #chi2.
   */
  int num_degree_of_freedom;

  /** Total continuum-optimizer iterations used, summed over regions of interest.  A numerical
   health diagnostic.
   */
  size_t num_continuum_iterations;

  /** How often the continuum optimizer had to fall back to an alternate starting estimate or to
   its derivative-free method.  Should normally be zero.
   */
  size_t num_continuum_restarts;

  /** Human-readable notes about anything the calculation had to change about the input; most
   commonly that overlapping regions of interest were combined.  Empty in the ordinary case.
   */
  std::vector<std::string> warnings;

  /** Future-measurement exposure divided by the reference's, or 1 for every other measurement
   model.

   \sa DeconMeasurementModel::BackgroundReference, fit_peaks
   */
  double exposure_ratio = 1.0;

  /** The trial peaks, in the exposure of #DeconComputeInput::measurement - that is, of the spectrum
   these will be drawn over - so a peak's Gaussian and its continuum are always in the same units.

   Under `DeconMeasurementModel::BackgroundReference` that means they show what the trial activity
   would have looked like *in the reference measurement*, not in the future one being predicted.
   Multiply an amplitude by #exposure_ratio for counts in the predicted measurement; do not scale
   the continuum with it, which belongs to the reference.
   */
  std::vector<PeakDef> fit_peaks;
};//struct DeconComputeResults

/** Computes the most consistent peaks, and their chi2, for a given input activity and distance.
 
 Throws exception if input is invalid, or error during calculation.
 */
DeconComputeResults decon_compute_peaks( const DeconComputeInput &input );


/** Which sentence describes a completed profile scan.

 `DeconActivityOrDistanceLimitResult::limitText` is hard-coded English because this file has no Wt
 dependency and batch mode never loads message bundles - but the GUI has to say the same thing in
 the user's language.  Rather than have each display re-derive the case from the flags (three
 places, three chances to drift), the scan records which sentence it chose and every consumer
 switches on this.

 \sa decon_limit_text_kind
 */
enum class DeconLimitTextKind : int
{
  /** Nothing could be stated; see `DeconActivityOrDistanceLimitResult::errorMessage`. */
  None = 0,

  /** "Less than U at CL (one-sided upper limit)" - an upper bound on the loaded spectrum. */
  OneSidedUpper,

  /** "Predicted sensitivity: less than U at CL for a T s measurement."

   Only when `DeconActivityOrDistanceLimitResult::is_predicted_sensitivity`; this must never be
   worded as an upper limit on the loaded spectrum.  \sa DeconMeasurementModel
   */
  PredictedSensitivity,

  /** "Between L and U at CL central" - only when `limitType == DeconLimitType::CentralInterval`. */
  CentralInterval,

  /** "Between 0 and U at CL central" - a central interval whose lower endpoint ran into the
   physical boundary at zero, so only the upper crossing exists.

   This is the ORDINARY outcome of asking for a central interval on a spectrum with no signal, and
   it must not be worded as a one-sided upper limit: the root was found at
   `quantile(chi2(1), CL)`, so quoting it as a one-sided bound at CL understates its coverage
   (a 95% central threshold gives a 97.5% one-sided bound).  Reporting it as the interval [0, U]
   is conservative - a parameter bounded below at zero gives coverage of at least CL.
   */
  CentralIntervalUpperBound,

  /** "Two one-sided CL bounds: L to U (~(2*CL-1) central coverage)."

   Two roots of a *one-sided* threshold are not a central interval; at 95% one-sided the pair covers
   only about 90% centrally.  This is the case the pre-Increment-C text mislabelled as
   "Between L and U at 95% CL".
   */
  TwoOneSidedBounds,

  /** "Distance >= L at CL" - a distance scan that bracketed only a lower crossing. */
  DistanceLowerBound
};//enum class DeconLimitTextKind


/** Classifies a completed scan into the sentence that describes it.

 Called once by `get_activity_or_distance_limits(...)` when it builds its English text, and
 available to display code that has a result and wants the localized equivalent.  Keeping the
 decision in one function is what stops the English and translated wordings from describing
 different things.

 Note the two scan types are not symmetric.  An activity scan reports from its *upper* crossing, so
 an upper crossing alone is a complete answer.  A distance scan reports from its *lower* crossing -
 "you could be at least this far away and still see it" - so for it an upper crossing alone states
 nothing, and returns #None.

 Consequence worth knowing: `DeconMeasurementModel::BackgroundReference` forces `found_lower_cl`
 false (the null expectation has no excess), so a *distance* limit under that model can never be
 reported.  The combination is therefore not offered in the UI.
 */
DeconLimitTextKind decon_limit_text_kind( const bool is_dist_limit,
                                          const bool found_lower_cl,
                                          const bool found_upper_cl,
                                          const bool is_predicted_sensitivity,
                                          const DeconLimitType limit_type );


struct DeconActivityOrDistanceLimitResult
{
  /** The input variables to the calculations */
  bool isDistanceLimit;
  double confidenceLevel;

  /** Which quantity was asked for; determines the threshold used, and how the result must be
   described.  \sa DeconLimitType
   */
  DeconLimitType limitType = DeconLimitType::OneSidedUpperLimit;

  /** The change in the profiled statistic that defined #upperLimit / #lowerLimit. */
  double limitDelta = 0.0;

  /** True when the result predicts a future measurement's sensitivity rather than bounding signal
   in the spectrum that was loaded, i.e. when the input used
   `DeconMeasurementModel::BackgroundReference`.

   Exists so that no display or report can call this an upper limit on the current spectrum.  When
   true, #foundLowerCl is always false: the null expectation has no excess to exclude zero with.
   */
  bool is_predicted_sensitivity = false;

  /** LIVE time, in seconds, of the future measurement #is_predicted_sensitivity refers to; zero
   otherwise.  This is what the likelihood actually used - see `DeconComputeInput::sample_exposure`.

   Do not quote this to a user; quote #sampleRealTime.  \sa plan_measurement
   */
  double sampleExposure = 0.0;

  /** REAL time, in seconds, of that same future measurement - the number the user entered, and the
   one every result string quotes.

   Derived from #sampleExposure using the reference spectrum's dead-time fraction, so the two
   describe one measurement.  A detector reporting only one of the two times is taken to have zero
   dead time, making them equal.
   */
  double sampleRealTime = 0.0;

  double minSearchValue;
  double maxSearchValue;
  DeconComputeInput baseInput;

  // TODO: refactor generating this text to a separate function.
  std::string limitText;
  std::string quantityLimitStr;
  std::string bestCh2Text;

  /** Which statement #limitText makes, so a display can render the same statement localized rather
   than re-deriving the case from the flags.  \sa DeconLimitTextKind
   */
  DeconLimitTextKind limitTextKind = DeconLimitTextKind::None;

  /** Non-empty when the scan completed but could not bracket the requested limit. */
  std::string errorMessage;

  /** Notes about anything that qualifies the result but did not prevent it: overlapping regions of
   interest being combined, or a profile that crosses the threshold more than twice.  These should
   be surfaced to the user alongside the limit, not dropped.
   */
  std::vector<std::string> warnings;
  
  /** The smallest value of the profiled statistic found by the scan.

   Named "chi2" for continuity with the display and the report keys, but it is a Poisson deviance
   (Cash statistic), not a chi-square - see `DeconComputeResults::chi2`.  Only *differences* between
   trial values are meaningful; the difference is what `decon_limit_delta` thresholds.  Its absolute
   value has no goodness-of-fit meaning: `E[chi2] != nbin - npar`, so `chi2/DOF` is not a fit-quality
   measure.
   */
  double overallBestChi2;
  double overallBestQuantity;
  std::shared_ptr<const DeconComputeResults> overallBestResults;

  bool foundUpperCl;
  double upperLimit;
  /** The profiled statistic at #upperLimit; a Cash statistic, see #overallBestChi2.
   Meaningful only when #foundUpperCl - otherwise it is a stale or default value.
   */
  double upperLimitChi2;
  std::shared_ptr<const DeconComputeResults> upperLimitResults;

  bool foundLowerCl;
  double lowerLimit;
  /** The profiled statistic at #lowerLimit; a Cash statistic, see #overallBestChi2.
   Meaningful only when #foundLowerCl.
   */
  double lowerLimitChi2;
  std::shared_ptr<const DeconComputeResults> lowerLimitResults;
  
  bool foundUpperDisplay = false;
  double upperDisplayRange;
  bool foundLowerDisplay;
  double lowerDisplayRange;
  
  /** The profiled statistic at a series of trial values, as {quantity, chi2}, for the displayed
   chart.  Each value is a Cash statistic, see #overallBestChi2 - which is why the chart is only
   ever read as "how far above the minimum", never as an absolute fit quality.
   */
  std::vector<std::pair<double,double>> chi2s;
  
  DeconActivityOrDistanceLimitResult();
};//struct DeconActivityOrDistanceLimitResult
  
  
/** Profiles activity (or distance) and returns where the statistic rises by the threshold that
 \p limit_type defines at \p wantedCl.

 @param wantedCl The confidence level; must be strictly between 0.5 and 1.
 @param limit_type Which quantity to report.  Defaults to a one-sided upper limit, which is what
        every caller wanted before this parameter existed, so the default does not change any
        previously reported number.
 */
/** How much of a scan's output the caller actually needs.

 A scan does considerably more than find the limit: it also fills a ~26-point chi2 display grid and
 three more full peak fits (best / lower / upper), each of which is a complete `decon_compute_peaks`
 call.  That is well over half the work, and a caller that only reads `upperLimit` pays for all of
 it.  \sa get_activity_or_distance_limits
 */
enum class DeconLimitDetail
{
  /** The limit, the chi2 display grid, and the best/lower/upper peak fits - everything a GUI needs
   to draw its chart.  The default, and what every interactive caller wants.
   */
  Full,

  /** The limit and nothing else: no display grid, no result peak fits.

   Also runs its internal searches serially rather than on a `SpecUtilsAsync::ThreadPool`.  The only
   caller is the projection Monte Carlo, which runs hundreds of these and supplies its own
   parallelism a level up; nesting a second pool there both wastes threads and, on builds where
   `SpecUtilsAsync` submits into the Wt io service, competes with request handling.
   */
  LimitOnly
};//enum class DeconLimitDetail


DeconActivityOrDistanceLimitResult get_activity_or_distance_limits( const double wantedCl,
                        const std::shared_ptr<const DeconComputeInput> base_input,
                        const bool is_dist_limit,
                        const double min_search_quantity,
                        const double max_search_quantity,
                        const bool useCurie,
                        const DeconLimitType limit_type = DeconLimitType::OneSidedUpperLimit,
                        const DeconLimitDetail detail = DeconLimitDetail::Full );


/** Returns a copy of `input` whose channel counts, live_time, and real_time are linearly
 scaled by `new_real_time / input->real_time()`, so the Measurement appears to have been
 taken for `new_real_time` seconds at the same per-channel rate and the same dead-time
 fraction.

 Used by the MDA / detection-confidence tools to answer "what would the limit be if my
 spectrum had been taken for a different dwell time?".  Counts are scaled to the expected
 value at the new dwell; the calc downstream then treats them as Poisson with variance equal
 to the (scaled) count.

 That gets the *median* right and the *spread* wrong, and the two should not be confused.  Variance
 equal to the scaled count asserts \p input's own rate is known exactly, when in fact the projection
 inherits \p input's counting noise: the projected limit's relative scatter goes as one over the
 root of the counts \p input holds, while the future limit's goes as one over the root of the counts
 the future measurement will hold.  With `k = new_real_time / input->real_time()` those differ by
 about `sqrt(k)`, so the honest predictive spread is about `sqrt(1 + k)` times what the projection
 implies.

 Measured by `RealSpectrumDwellProjectionAccuracy` on an 11.5 hour real background, on the
 detection limit - the quantity the tools quote when a measurement time is entered: the median
 projected Currie limit was within 4% of what the dwell actually delivered over k from 0.25 to 16
 (ratios 1.037 / 1.010 / 0.995 / 1.011), while the measured spread understatement ran
 1.1 / 1.6 / 1.7 / 4.0 over the same k.  Immaterial up to about `k = 1`; a planning estimate rather
 than a bound past it.  This is what `dl-plan-time-tt` and the two help pages tell the user, in the
 same terms.

 Some detectors report only one of the two times, so a non-positive `real_time()` falls back to
 `live_time()` - treating that as zero dead time rather than refusing to project.

 Throws if `input` is null, `new_real_time <= 0`, or both of the input's times are non-positive
 (NaN included: the checks are written `!(x > 0)` so NaN is rejected rather than let through).
 */
std::shared_ptr<const SpecUtils::Measurement>
scale_spectrum_for_dwell( const std::shared_ptr<const SpecUtils::Measurement> &input,
                          const float new_real_time );


/** The prior strength added to each channel's observed count before its rate is drawn, in
 `draw_projected_measurement`.  Jeffreys, and not an arbitrary choice - see that function.
 */
double projected_limit_prior_strength();


/** One Monte Carlo realisation of what a measurement of \p new_real_time *would* record, given that
 \p reference recorded what it did.

 Where `scale_spectrum_for_dwell` multiplies the counts by `k` and hands back the expectation, this
 draws an actual future measurement, per channel, in two stages:

 \code
   rate  ~ Gamma( n + alpha, 1 )    // the channel's true rate, given it recorded n
   counts ~ Poisson( k * rate )     // what a dwell of k times the reference then records
 \endcode

 The two stages are the two things a projection has to account for and a plain scaling accounts for
 neither: we do not know the rate, and the future measurement will fluctuate around it.  Together
 they give `Var = k(n+alpha)(1+k)` against the `k*n` a scaling asserts - the `sqrt(1+k)` the
 user-facing text describes.

 \p alpha is Jeffreys (0.5), and the value matters more than it looks.  At `alpha = 0` this is a
 plain Poisson-of-Poisson with exactly the right mean, but a channel that happened to record zero
 is then locked at zero in every realisation, which is most channels in the low-count regime a
 detection limit lives in.  At `alpha = 1` (flat prior) nothing locks, but a whole count is added to
 *every* channel, and a region holding a couple of counts per channel is inflated by tens of
 percent, which the limit inherits as its square root.  Measured over 400 reference/future pairs on
 a real background, the median predicted Currie limit ran 0.93 -> 1.18 across `k` at `alpha = 0` and
 1.04 -> 1.25 at `alpha = 1`, against 1.06 / 1.01 / 1.09 / 1.09 at Jeffreys.
 \sa ProjectedMeasurementMonteCarlo

 Only channels in `[first_channel, last_channel]` are drawn; every other channel is given its
 expectation `k*n`.  Callers must therefore pass a window covering every channel the limit will
 read - the regions of interest and their side channels.  Drawing all of a 16k-channel spectrum
 costs about 200x what the limit itself does, so this is what makes the Monte Carlo affordable at
 all rather than an optimisation.

 \p generator is advanced; callers wanting a reproducible limit must seed it deterministically from
 the input (a reported MDA that changes when nothing changed is worse than one that is slightly
 wrong).

 Throws on the same conditions as `scale_spectrum_for_dwell`, plus an out-of-range channel window.
 */
std::shared_ptr<const SpecUtils::Measurement>
draw_projected_measurement( const std::shared_ptr<const SpecUtils::Measurement> &reference,
                            const float new_real_time,
                            const size_t first_channel,
                            const size_t last_channel,
                            std::mt19937 &generator );


/** A limit predicted for a measurement that has not been taken, as a distribution rather than a
 point.

 \sa currie_projected_limit, decon_projected_limit
 */
struct ProjectedLimit
{
  /** Whether #median and the band below are usable at all. */
  bool valid = false;

  /** The median over Monte Carlo realisations.

   Both tools use this as the *scale* for the band rather than as the reported limit: they show
   `lower/median` and `upper/median` as multipliers on the limit they already computed directly.
   Quoting #median itself would make the displayed number jitter with the draw seed, for a
   difference the measurement puts at a few percent - see `RealSpectrumDwellProjectionAccuracy`.
   The consequence, which the help text states, is that the displayed band is centred on the
   directly computed limit and not exactly on this median.
   */
  double median = 0.0;

  /** The 16th and 84th percentiles over realisations: what the dwell could plausibly deliver.
   Measured to cover near their nominal rate on a real background - 0.65-0.74 against 0.68, and
   0.88-0.94 against 0.90 for the 5th-95th - for both methods and over a 64-fold span of dwell
   ratios.  \sa ProjectedMeasurementMonteCarlo
   */
  double lower = 0.0;
  double upper = 0.0;

  /** Realisations that produced a limit, and those attempted.  A cell where many trials failed is
   reported rather than silently averaged over the survivors.
   */
  size_t num_used = 0;
  size_t num_attempted = 0;
};//struct ProjectedLimit


/** The Currie limit for a measurement of \p planned_real_time, as a predicted distribution.

 \p input describes the region on the *reference* spectrum (`input.spectrum` is the reference); each
 realisation redraws the region and its side channels and recomputes `currie_mda_calc`.  Cheap - the
 calculation is closed-form arithmetic, so a few hundred realisations cost milliseconds.

 Scores `CurrieMdaResult::detection_limit`, because that is what the tools quote when a measurement
 time is entered.  Scoring `upper_limit` instead would describe a different and wider distribution,
 since that one also carries the observed peak-region excess of each realisation.

 Note this differs from `decon_projected_limit`, which scores the profile *upper* limit - and it
 differs for the same reason, that each band describes the number displayed beside it.  The visible
 consequence is that the two bands move in opposite directions as the dwell grows: this one narrows
 towards the floor set by the reference's own counting statistics, while the deconvolution band
 widens, because the excess it carries grows with the exposure being predicted.  That is not an
 inconsistency between the two; it is the difference between "how well do I know my sensitivity"
 and "how much could the measurement come back saying".

 Returns an invalid result rather than throwing if too few realisations produced a limit.
 */
ProjectedLimit currie_projected_limit( const CurrieMdaInput &input,
                                       const double planned_real_time,
                                       const size_t num_trials );


/** Which limit each Monte Carlo realisation is scored as.  These are different questions, and a
 prediction is only meaningful against the one the caller will actually go on to compute.
 \sa decon_projected_limit
 */
enum class ProjectedLimitScoring
{
  /** Each realisation analysed on its own, as an analyst does with a single loaded spectrum.
   Predicts what `DeconMeasurementModel::CurrentSpectrum` will return for that measurement.
   */
  SampleOnly,

  /** Each realisation scored jointly with the reference, sharing one continuum.  Predicts what
   `DeconMeasurementModel::BackgroundReference` describes - so this is the one to use when giving a
   spread for that mode, and it is legitimately tighter than #SampleOnly because the reference
   contributes information.
   */
  JointWithReference
};//enum class ProjectedLimitScoring


/** The deconvolution limit for a measurement of \p planned_real_time, as a predicted distribution.

 Each realisation is a full profile scan, so this is orders of magnitude more expensive than the
 Currie equivalent and is run across \p num_threads threads.  \p base_input carries the reference
 spectrum; whether each realisation is then scored alone or together with that reference is
 \p scoring, and the two answer different questions - see `ProjectedLimitScoring`.

 The realisations are drawn up front on the calling thread and written back by index, so the result
 does not depend on \p num_threads: one thread and eight give the same numbers.

 Set \p cancel to abandon the work: it is checked between realisations, so a caller that starts one
 of these on every input change does not leave superseded runs burning a thread to completion.  A
 cancelled run returns an invalid result.

 Returns an invalid result rather than throwing if too few realisations produced a limit.
 */
ProjectedLimit decon_projected_limit( const std::shared_ptr<const DeconComputeInput> &base_input,
                                      const double wantedCl,
                                      const double planned_real_time,
                                      const double max_search_quantity,
                                      const bool useCurie,
                                      const DeconLimitType limit_type,
                                      const ProjectedLimitScoring scoring,
                                      const size_t num_trials,
                                      const size_t num_threads,
                                      const std::shared_ptr<std::atomic_bool> &cancel = nullptr );


/** What one user-entered "planned measurement time" resolves to for each calculation.
 \sa plan_measurement
 */
struct PlannedMeasurement
{
  /** What Currie-style calculations use: the reference projected to the planned time, or the
   reference itself when no planned time was requested.
   */
  std::shared_ptr<const SpecUtils::Measurement> currie;

  /** What `DeconComputeInput::measurement` must be.  The same object as #currie, except under
   `DeconMeasurementModel::BackgroundReference`, where it is the *unscaled* reference - the
   reference counts have to stay real counts at their own exposure for their counting statistics to
   propagate.
   */
  std::shared_ptr<const SpecUtils::Measurement> decon;

  /** `DeconComputeInput::sample_exposure`, a LIVE time in seconds.  Zero unless the model is
   `BackgroundReference` and a planned time was given.
   */
  double sample_exposure = 0.0;

  /** The planned measurement time as the user entered it, in seconds of REAL time; zero when no
   projection was requested.  Users think in real time and results are reported in real time, so
   this - not #sample_exposure - is the number to quote back.
   */
  double planned_real_time = 0.0;

  /** `sample_exposure / decon->live_time()`, or 1.0.  The factor taking a per-Bq yield computed at
   the reference exposure to the projected measurement; display code that pairs a Currie number
   with a per-Bq yield needs it, because the two are then at different exposures.
   */
  double exposure_ratio = 1.0;
};//struct PlannedMeasurement


/** Resolves the single "planned measurement time" control into what each calculation needs.

 There is one control, meaning one thing - "the dwell I am asking about" - but it reaches the
 calculation two different ways, and mixing them silently changes every number:

 - Currie always, and the deconvolution method under `DeconMeasurementModel::CurrentSpectrum`:
   the spectrum itself is projected with `scale_spectrum_for_dwell`, and #sample_exposure is zero.
 - The deconvolution method under `DeconMeasurementModel::BackgroundReference`: the spectrum is left
   alone and the time is carried in `DeconComputeInput::sample_exposure` instead.  Scaling as well
   would apply the projection twice, since `decon_compute_peaks` already multiplies the trial signal
   by `sample_exposure/live_time`.

 Note `DeconComputeInput::sample_exposure` is compared against `Measurement::live_time()`, so it is
 a LIVE time, whereas \p planned_real_time_s is entered as a REAL time.  The conversion here
 preserves the reference's dead-time fraction, which makes

 \code
   result.sample_exposure == scale_spectrum_for_dwell( reference, planned_real_time_s )->live_time()
 \endcode

 exactly - and that identity is what makes the Currie and deconvolution paths describe the *same*
 planned measurement rather than two dwells differing by the dead-time fraction.

 A non-positive \p planned_real_time_s means "not requested", and returns the reference unchanged
 for both spectra.  Throws only what `scale_spectrum_for_dwell` throws.
 */
PlannedMeasurement plan_measurement( const std::shared_ptr<const SpecUtils::Measurement> &reference,
                                     const double planned_real_time_s,
                                     const DeconMeasurementModel model );


/** The two band endpoints of \p limit, as multiples of its median, ready to put on screen.

 Both tools show the band as a multiplier on the limit they display, so this is the one place that
 decides how it reads.  Fixed decimals rather than significant figures, because the two endpoints
 have to be on the same scale to be read as a pair: a +0.4% band formatted to two significant
 figures gives "0.99 to 1", which looks like a bound at exactly one rather than the tight band it
 is.  Two decimals gives "0.99 to 1.00", and stays legible for wide bands too ("0.61 to 1.72").

 Returns false, leaving the strings untouched, when \p limit is unusable or when both endpoints
 render identically - a band of "1.00 to 1.00" carries no information and only adds clutter.
 */
bool projected_band_endpoints( const ProjectedLimit &limit,
                               std::string &lower_multiple,
                               std::string &upper_multiple );

}//namespace DetectionLimitCalc


#endif //DetectionLimitCalc_h
