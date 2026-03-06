#ifndef PeakFitSpecImp_h
#define PeakFitSpecImp_h
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

// This header contains internal implementation details for peak finding and
// FWHM-based detector classification. Only PeakFitUtils.cpp and GA optimization
// code (target/peak_fit_improve/) should include this header. All other InterSpec
// code should use coarse_det_type() from PeakFitUtils.h.
//
// Stage 2 (NaI/CsI vs MedRes) was attempted but performed poorly:
//   39.6% correct, 1.7% incorrect, 58.6% unknown across 119,812 spectra
//   Optimized settings:
//     nai_fwhm_bias_mult=1.319, medres_fwhm_frac=0.511, hpge_fwhm_bias_mult=1.196
//     nai_distance_weight=2.082, medres_distance_weight=4.679
//     amp_clamp_denom=18.257, weight_clamp_max=6.976
//     min_peak_energy=63.985, max_peak_energy=3782.517
//     unknown_threshold=0.762, min_peaks_for_classify=3
//     min_peak_significance=7.310, sig_fraction_of_max=0.059

#include <memory>
#include <string>
#include <vector>

#include "InterSpec/PeakFitUtils.h"

class PeakDef;

namespace SpecUtils
{
  class Measurement;
}


namespace PeakFitSpec
{

/** Result of FWHM-based HPGe vs non-HPGe classification.
 Only 3 outcomes are possible from the spectral classifier.
 Converted to CoarseResolutionType at the coarse_det_type() boundary.
 */
enum class SpecClassType : int
{
  LowOrMedRes,  // Not HPGe (NaI, CsI, LaBr, CZT, etc.)
  High,         // HPGe
  Unknown
};


/** Settings for the candidate peak finder (Savitzky-Golay second-derivative method). */
struct FindCandidateSettings
{
  int num_smooth_side_channels = 9;
  int smooth_polynomial_order = 2;
  double threshold_FOM = 0.758621;
  double more_scrutiny_FOM_threshold = 1.598265;
  float pos_sum_threshold_sf = 0.119178f;

  /** For second-derivative, how many channels are required to be above threshold,
   in-order to signal a transition.
   */
  size_t num_chan_fluctuate = 1;

  float more_scrutiny_coarser_FOM = 3.001943f;

  /** The minimum Chi2 required, of at least one channel in ROI, to be above a straight
   line predicted by the channels on either side of the ROI.
   */
  float more_scrutiny_min_dev_from_line = 6.816465;

  float amp_to_apply_line_test_below = 6.000000;

  // Energy-adaptive smoothing: SG window scales with channel position to match peak width
  //  vs energy.  Below ref channel (smooth_ref_fraction * nchannel), uses num_smooth_side_channels.
  //  Above ref channel, side_bins = num_smooth_side_channels * (ch/ref_ch)^smooth_scale_power.
  //  Set smooth_scale_power = 0 to disable (use fixed num_smooth_side_channels everywhere).
  float smooth_ref_fraction = 0.0f;         // Reference channel as fraction of nchannel (e.g. 0.1)
  float smooth_scale_power = 0.0f;          // Scaling exponent (e.g. 0.5 = sqrt for NaI)

  // Minimum significance of 2nd derivative at peak center to accept the candidate.
  //  Significance = |second_deriv[min]| / sqrt(variance[min]).  Set to 0 to disable.
  float min_second_deriv_significance = 0.0f;

  // Compton backscatter asymmetry test: checks that positive 2nd-derivative regions flanking
  //  the peak are reasonably symmetric. Set compton_next_ratio_max >= 50 to disable.
  float compton_next_ratio_max = 100.0f;    // Old algo: 4.0
  float compton_prev_ratio_min = 0.0f;      // Old algo: 0.2
  float compton_total_ratio_min = 0.0f;     // Old algo: 0.3

  // Low-energy drop-off test: validates peaks below a threshold energy have a Gaussian-like
  //  shape, not a detector efficiency edge. Set low_energy_test_max_keV = 0 to disable.
  float low_energy_test_max_keV = 0.0f;     // Old algo: 130 keV
  float low_energy_drop_fraction = 0.45f;   // Minimum actual/expected Gaussian drop ratio

  // PCGAP-style ROI: walking background feature detection, blended with derivative-based ROI.
  //  pcgap_roi_blend_weight: 0 = only derivative-based ROI, 1 = only PCGAP ROI.
  //  Blending applies only when FOM < pcgap_fom_blend_threshold.
  float pcgap_feature_nsigma = 3.0f;        // Feature detection threshold (sigma)
  float pcgap_max_extent_nsigma = 7.5f;     // Maximum ROI walk distance from mean (in sigma)
  float pcgap_roi_blend_weight = 0.0f;      // 0 = disabled (current behavior)
  float pcgap_fom_blend_threshold = 5.0f;   // FOM below which PCGAP blending applies

  std::string print( const std::string &var_name ) const;
  std::string to_json() const;
};//struct FindCandidateSettings


/** Finds candidate peaks in a spectrum using second-derivative analysis.

 Uses a Savitzky-Golay smoothed second derivative to identify peak-like features,
 then estimates sigma via a two-pass log-parabola Gaussian fit.

 @param data The spectrum measurement
 @param start_channel First channel to search (0 for beginning)
 @param end_channel Last channel to search (0 for end)
 @param settings The candidate peak finder parameters
 @return Vector of candidate PeakDef objects with estimated mean, sigma, amplitude, and ROI
 */
std::vector<PeakDef> find_candidate_peaks(
  const std::shared_ptr<const SpecUtils::Measurement> &data,
  size_t start_channel,
  size_t end_channel,
  const FindCandidateSettings &settings );


/** Settings for FWHM-based HPGe vs non-HPGe classification, optimized via GA. */
struct LowHighResClassifySettings
{
  // Multiplicative biases on reference FWHM curves (nominal ~1.0, range [0.7, 1.3])
  double hpge_fwhm_bias_mult = 1.0;
  double nai_fwhm_bias_mult = 1.0;

  // Per-type distance weights (range [0.1, 5.0])
  double hpge_distance_weight = 1.0;
  double nai_distance_weight = 1.0;

  // Peak weighting: weight = min( sqrt(amplitude) / amp_clamp_denom, weight_clamp_max )
  double amp_clamp_denom = 10.0;         // range [1, 50]
  double weight_clamp_max = 10.0;        // range [1, 30]

  // Energy filtering
  double min_peak_energy = 50.0;         // range [30, 100] keV
  double max_peak_energy = 3000.0;       // range [2000, 5000] keV

  // Confidence threshold: if best type's vote fraction < this, return Unknown
  double unknown_threshold = 0.4;        // range [0.2, 0.8]

  // Min peaks required (after filters), else Unknown
  int min_peaks_for_classify = 2;        // range [1, 5]

  // Min peak significance (sqrt(peakArea), i.e., Poisson SNR) to include in voting
  double min_peak_significance = 2.0;    // range [0.5, 5.0]

  // Relative significance threshold: only consider peaks with significance
  // >= sig_fraction_of_max * max_peak_significance. Filters noise-floor peaks.
  double sig_fraction_of_max = 0.05;     // range [0.005, 0.3]

  // Asymmetric FWHM penalty: extra multiplier on HPGe distance when peak is
  // narrower than HPGe reference. Modulated by significance: high-significance
  // narrow peaks get less penalty (could be real), low-significance narrow peaks
  // get more penalty (likely noise).
  double narrow_penalty_mult = 2.0;      // range [0.5, 10.0]
  double narrow_sig_threshold = 20.0;    // range [5, 100] - significance above which no extra penalty

  // Best-peak CZT check: for the most significant peak, compare distance to CZT
  // reference vs HPGe reference. If closer to CZT, apply penalty to HPGe vote.
  // Only applied when HPGe vote fraction is below best_peak_conf_threshold.
  // Helps discriminate CZT_H3D (whose peaks have HPGe-like widths) from true HPGe.
  double best_peak_czt_penalty = 1.0;          // range [0, 5] - penalty weight
  double best_peak_conf_threshold = 0.85;      // range [0.5, 1.0] - only apply when HPGe conf below this

  std::string print( const std::string &var_name ) const;
  std::string to_json() const;
};//struct LowHighResClassifySettings


/** Classifies candidate peaks as HPGe (High) vs non-HPGe (LowOrMedRes) using
 2-way FWHM soft voting with outlier rejection and best-peak CZT check.

 @param candidate_peaks Peaks from find_candidate_peaks
 @param settings The classification parameters
 @param confidence_out If non-null, receives the winning type's vote fraction [0,1]
 @return LowOrMedRes, High, or Unknown
 */
SpecClassType classify_lowres_highres_from_peak_candidates(
  const std::vector<PeakDef> &candidate_peaks,
  const LowHighResClassifySettings &settings,
  double *confidence_out = nullptr );


/** Convenience function that classifies a spectrum as HPGe vs non-HPGe using
 hardcoded GA-optimized settings.

 Internally finds candidate peaks and runs classify_lowres_highres_from_peak_candidates.

 @param measurement The spectrum to classify
 @param confidence_out If non-null, receives the winning type's vote fraction [0,1]
 @return LowOrMedRes, High, or Unknown
 */
SpecClassType initial_lowres_highres_classify(
  const std::shared_ptr<const SpecUtils::Measurement> &measurement,
  double *confidence_out = nullptr );

}//namespace PeakFitSpec

#endif //PeakFitSpecImp_h
