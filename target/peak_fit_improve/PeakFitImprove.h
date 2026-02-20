#ifndef PeakFitImprove_h
#define PeakFitImprove_h
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


#include <memory>
#include <string>
#include <vector>

namespace SpecUtils
{
  class SpecFile;
  class Measurement;
}


#define WRITE_ALL_SPEC_TO_HTML 0

namespace PeakFitImprove
{
const bool debug_printout = true;

const double debug_lower_energy = 250;
const double debug_upper_energy = 240;

extern size_t sm_num_optimization_threads;
extern size_t sm_num_threads_per_individual;
extern size_t sm_ga_population;
extern size_t sm_ga_generation_max;
extern size_t sm_ga_best_stall_max;
extern size_t sm_ga_elite_count;
extern double sm_ga_crossover_fraction;
extern double sm_ga_mutation_rate;
extern double sm_ga_mutate_threshold;
extern double sm_ga_crossover_threshold;
}



struct PeakTruthInfo
{
  float energy;
  float cps;
  float area;
  float fwhm;

  /**
   Where:
   const double SC = 2.0 + 0.02*(lowSkew + highSkew)*std::pow(energy/661.0, resolutionPower );
   FullWidth = SC * GetFWHM(energy, resolutionOffset, resolution661, resolutionPower );
   */
  float full_width;

  /**
   "S.E.", "D.E.", "EscapeXRay", "Peak"
   */
  std::string label;

  PeakTruthInfo();

  explicit PeakTruthInfo( const std::string &line );
};//struct PeakTruthInfo


struct InjectSourceInfo
{
  std::string file_base_path;
  std::string src_name;

  std::vector<PeakTruthInfo> source_lines;
  std::vector<PeakTruthInfo> background_lines;

  /** The PCF file for this source. */
  std::shared_ptr<SpecUtils::SpecFile> spec_file;

  /** The source + background spectra, that are Poisson varied. */
  std::vector<std::shared_ptr<const SpecUtils::Measurement>> src_spectra;

  /** Background that is Poisson varied, and same duration as src spectra. */
  std::shared_ptr<const SpecUtils::Measurement> short_background;

  /** A longer background, with Poisson variation. */
  std::shared_ptr<const SpecUtils::Measurement> long_background;

  /** Source + background without Poisson variation. */
  std::shared_ptr<const SpecUtils::Measurement> src_no_poisson;

  /** Background, with no Poisson variation. */
  std::shared_ptr<const SpecUtils::Measurement> background_no_poisson;

  /** Whether this data was loaded from the compact data format (4 measurements per PCF,
   pre-computed truth info) vs the original format (16 measurements per PCF, CSV truth).
   */
  bool from_compact_data = false;
};//struct InjectSourceInfo


struct DetectorInjectSet
{
  std::string detector_name;
  std::string location_name;
  std::string live_time_name;

  std::vector<InjectSourceInfo> source_infos;
};//struct DetectorInjectSet



struct ExpectedPhotopeakInfo
{
  /** An approximate detection sigma, i.e., 2.33 times this value is kinds 95% DL */
  double nsigma_over_background;

  double roi_lower;
  double roi_upper;
  double effective_energy;
  double peak_area;
  double continuum_area;
  double effective_fwhm;

  std::vector<PeakTruthInfo> gamma_lines;
};//struct ExpectedPhotopeakInfo


struct DataSrcInfo
{
  std::string detector_name;
  std::string location_name;
  std::string live_time_name;

  InjectSourceInfo src_info;

  std::vector<ExpectedPhotopeakInfo> expected_photopeaks;
  std::vector<ExpectedPhotopeakInfo> expected_signal_photopeaks;
  std::vector<ExpectedPhotopeakInfo> expected_background_photopeaks;
};//struct DataSrcInfo



/** The weights and limits to apply to the various optimizations.

 These are all chosen using 'expert judgment', which means they may not be that well chosen.
 */
namespace JudgmentFactors
{
  // Lets not worry about super-small peaks, even where there is little to no background
  const double min_truth_peak_area = 5;

  //Photopeak clusters below this next number of sigma will be discarded, since we really
  //  shouldnt find these peaks.
  //  i.e., The threshold at which we will start punishing if a peak is not expected at this
  //  threshold; these source photopeaks have already been removed.
  const double min_truth_nsigma = 1.0;

  const double def_want_nsigma = 4;   // i.e., above 4 sigma, lets weight all peaks the same
  const double min_def_wanted_counts = 15; //i.e., if expected peak area is below 15 counts, we wont punish for not finding
  const double lower_want_nsigma = 2; // The number of sigma above which we will positively reward finding a peak
  // Between `def_want_nsigma` and `lower_want_nsigma` we will linearly weight for not finding a peak

  const double found_extra_punishment = 0.25; // 1/this-value gives the trade-off of finding extra peaks, verses not finding peaks

  // When a peak between lower_want_nsigma and def_want_nsigma is found, the minimum value we should assign
  const double min_initial_fit_maybe_want_score = 0.25;

  // Cost for fitting an extra peak after the initial proper fit
  const double initial_fit_extra_peak_punishment = 0.75;

  // Multiple for the fraction of additional candidates tried, that didnt stick around
  //  E.g., with a value of 1.0, if 100% of tried peaks failed, this would be equivalent of not
  //    finding one peaks
  const double extra_add_fits_punishment = 1.0;
}//namespace JudgmentFactors


struct FindCandidateSettings
{
  int num_smooth_side_channels = 9; // low res more
  int smooth_polynomial_order = 2;  // highres 3, lowres 2
  double threshold_FOM = 0.758621;  // accept peaks higher than this FOM
  double more_scrutiny_FOM_threshold = 1.598265; // Peaks bellow this get extra scrutiny
  float pos_sum_threshold_sf = 0.119178f;

  /** For second-derivative, how many channels are required to be above threshold, in-order to signal a transition */
  size_t num_chan_fluctuate = 1;

  //float min_counts_per_channel = 1.0f;
  float more_scrutiny_coarser_FOM = 3.001943f;

  /** The minimum Chi2 required, of at least one channel in ROI, to be above a straight
   line predicted by the channels on either side of the ROI.
   */
  float more_scrutiny_min_dev_from_line = 6.816465;

  float amp_to_apply_line_test_below = 6.000000;

  std::string print( const std::string &var_name ) const;

  std::string to_json() const;
};//struct FindCandidateSettings


struct InitialPeakFindSettings
{
  /** Filter significance, after initial fit, using `chi2_significance_test(...)` function

   These two parameters should be optimized separate from everything else

   `initial_stat_threshold`: This is how incompatible with background/continuum the data
   must be, before a peak is allowed to exist. A reasonable search range of values is maybe between 0 and 8.

   `initial_hypothesis_threshold`:  this specifies how well the peak must match in shape to a gaussian
   in order to keep the peak.  The higher this number, the more like a gaussian it must be. It is actually the ratio of
   the null hypothesis chi2 to the test hypothesis chi2.  A reasonable search range of values is maybe between -0.1 and 10.
   */
  double initial_stat_threshold = 1.951264;
  double initial_hypothesis_threshold = 0.673169;

  double initial_min_nsigma_roi = 2.246770;
  double initial_max_nsigma_roi = 6.378162;

  enum class FwhmFcnForm : int {
    Gadras,

    SqrtPolynomialTwoCoefs,  //FWHM = sqrt( Sum_i{A_i*pow(x/1000,i)} );
    SqrtPolynomialThreeCoefs,

    SqrtEnergyPlusInverse, // FWHM = `sqrt(A0 + A1*E + A2/E)`

    NumFwhmFcnForm
  };

  FwhmFcnForm fwhm_fcn_form = FwhmFcnForm::SqrtPolynomialTwoCoefs;


  /** Adds peaks based on now kinda known FWHM functional form.
   Lets just do something stupid, and slide along non-ROI areas, and have an ROI ~5 sigma wide,
   and draw a line using surrounding channels to fit a flat continuum, and if the deficit is
   some reasonable value (2.5 sigma? - need to do a time-tradeoff comparison to get something
   reasonable), then call `fit_amp_and_offset(...)` and `chi2_significance_test(...)`

   `search_roi_nsigma_deficit`: Reasonable search range of 1 to 10.
   `search_stat_threshold`: Reasonable search range of 0 and 8.
   `search_hypothesis_threshold`: Reasonable search range of 0 and 8. -0.1 and 10.
   `search_stat_significance`: Reasonable search range of 1 and 8.
   */
  double search_roi_nsigma_deficit = 4.241748;
  //double search_stat_threshold = 8.051485;
  double search_hypothesis_threshold = 3.342207;
  double search_stat_significance = 2.025582;

  /** Add peaks to ROIs.  WIP.
   `ROI_add_nsigma_required`: required previous and new peaks to be better than. Reasonable search range: 1 to 8
   `ROI_add_chi2dof_improve`: Reasonable search range: 0 to 8
   */
  double ROI_add_nsigma_required = 3.526017;
  double ROI_add_chi2dof_improve = 0.516229;

  std::string print( const std::string &var_name ) const;
};//struct InitialPeakFindSettings


struct FinalPeakFitSettings
{
  /** Combine ROIs, based on how near means are, and how much initial ROIs are overlapping.

   Reasonable search range of 1 to 10
   */
  double combine_nsigma_near_;

  /** Some very rough guesses as to what we should about expect

    < 0.8 FWHMMust combineSingle peak artifact
   0.8-2.0 FWHM - Should combineStrong interference
   2.0-3.0 FWHM - Prefer combineContinuum sharing
   3.0-4.0 FWHM - Either OKWeak interference
     > 4.0 FWHM - Prefer separateIndependent peaks
  */
  double require_combine_num_fwhm_near;
  double not_allow_combine_num_fwhm_near;


  /** Combine ROIs, based on how much initial ROIs are overlapping.
   
   Reasonable search range of -1 to 1
   */
  //double combine_ROI_overlap_frac;


  /** How many nsigma the peak needs to be to consider modifying the continuum type; e.g., how statistically significant
   the needs to be to have anything other than a linear continuum - applies to the most stat sig peak in a ROI

   Reasonable search range: 5 to 60. ???
   */
  double cont_type_peak_nsigma_threshold;

  /** How many stat higher the average channel is on the left than the right, where stat is based of side edge area, to consider a stepped continuum.

   Reasonable search range: 1 to 20. ???
   */
  double cont_type_left_right_nsigma;

  /** To see if flat, linear, or bilinear stepped continuum
   //Cont P1: -1.6985 on left, -0.3472 on right
   */
  //double cont_type_sum_slopes;

  /** Like above, but for stepped continua.
   If multiple peaks in ROI, then
   */
  //double stepped_roi_extent_lower_side_stat_multiple, stepped_roi_extent_upper_side_stat_multiple;

  /** How much of an improvement quadratic or cubic continuum type needs on the chi2/dof to be to actually use.

   Reasonable search range: 0 to 4. ???
   */
  double cont_poly_order_increase_chi2dof_required;

  /** How much of an improvement to the chi2/dof either moving to a stepped continuum type or increasing the
   * type of stepped continuum must give, in order to use the result

   Reasonable search range: 0 to 4. ???
   */
  double cont_step_type_increase_chi2dof_required;

  /** Check if we should add skew.
      Right now maybe just just ratio of area between continuum and data for -1 to -4 sigma from mean,
      compared to peak area, or maybe the ratio of sqrt of the areas.

   Reasonable search range: 0 to 10. ???
   */
  double skew_nsigma;

  /** After fitting a gaussian+continuum to data, for channels cooresponding between 1 and 1.5 FWHM, the difference between data and fit,
   divided by uncert, is summed, then divided by number of channels.  Large values indicate skew - I am guessing - but uncertain.

   This is an "or" to `skew_nsigma`.

   Reasonable range: ???
   */
  double left_residual_sum_min_to_try_skew, right_residual_sum_min_to_try_skew;

  /** How much adding skew needs to improve the ROI threshold in order to keep it.

   Reasonable search range: 0 to 4. ???
   */
  double skew_improve_chi2_dof_threshold;


  /** Determine ROI left and right base-widths for single peak ROIs.
   Width will be modified based on stat uncert of initial fit.

   Reasonable search range: 0.5 to 7 - nominal ~2.5 FWHM
  */
  double roi_extent_low_num_fwhm_base_highstat, roi_extent_high_num_fwhm_base_highstat;

  double roi_extent_low_num_fwhm_base_lowstat, roi_extent_high_num_fwhm_base_lowstat;

  double high_stat_threshold;

  /** The maximum amount extra FWHM beyond `roi_extent_low/high_num_fwhm_base` that the ROI can extend */
  double roi_extent_low_num_fwhm_extra;
  double roi_extent_high_num_fwhm_extra;

  /** We will walk away from the extent minimum, towards the extent maximum, but if we detect a second derivative
   that is above this statisticaly significant, we'll stop, and not include that channel.
   */
  double roi_end_second_deriv_thresh;

  /** If data is this number of statistical below away from fit continuum + peak, the multi-peak ROI will be broken up, but only if
   the resultant Chi2Dof improves by `break_multi_roi_up_required_chi2dof_improve`. */
  double break_multi_roi_up_continuum_away_sigma;
  double break_multi_roi_up_required_chi2dof_improve;

  std::string print( const std::string &var_name ) const;
};//struct FinalPeakFitSettings


#endif //PeakFitImprove_h
