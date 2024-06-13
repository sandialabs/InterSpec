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

#include <memory>
#include <ostream>

#include "InterSpec/PeakDef.h"

// Forward declarations
struct DetectorPeakResponse;

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
 - All calculations labeled as "currie" actually follow ISO 11929 methodologies for the simple gross-counts style calculation.
  E.g., A peak region is defined, and then a given number of channels are used on either side of the peak region to estimate the
  number of continuum counts inside the peak region, and variance, and then this is used to set limit from the total number or counts
  in the peak region.
  - Only single peak calculations, with no interferences are implemented.
 - All calculations labeled as "decon" use a more sophisticated  "de-convoluted" method that takes into account the shape of the peak
  and better takes into account all information provided, as well as using multiple peaks of an isotope to derive limits.  This methodology
  seems to follow the intent of Annex B of ISO 11929-3:2019, but instead these calculations form a large chi2/likelihood calculation to
  co-compute everything and hopefully do a better job.
  
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
  
}//namespace DetectionLimitCalc


#endif //DetectionLimitCalc_h
