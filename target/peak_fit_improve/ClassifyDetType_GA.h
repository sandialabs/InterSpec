#ifndef ClassifyDetType_GA_h
#define ClassifyDetType_GA_h
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

#include <string>
#include <vector>
#include <memory>
#include <functional>

#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitSpecImp.h"

struct DataSrcInfo;

class PeakDef;

namespace SpecUtils
{
  class Measurement;
}


/** Stage 1 settings: alias to the production LowHighResClassifySettings. */
using HighResClassifySettings = PeakFitSpec::LowHighResClassifySettings;


namespace ClassifyDetType_GA
{
  /** Returns the true detector type for a given detector name from the test data. */
  PeakFitUtils::CoarseResolutionType true_det_type_for_name( const std::string &detector_name );

  /** Returns human-readable string for a CoarseResolutionType. */
  const char *det_type_str( PeakFitUtils::CoarseResolutionType type );


  /** Stage 1 classification: High-resolution (HPGe) vs NaI.

   @param candidate_peaks Peaks from find_candidate_peaks
   @param settings The high-res classification parameters
   @param confidence_out If non-null, receives the winning type's vote fraction [0,1]
   @return High if HPGe, Low if NaI (to be refined by Stage 2), or Unknown
  */
  PeakFitUtils::CoarseResolutionType classify_highres_from_peaks(
    const std::vector<PeakDef> &candidate_peaks,
    const HighResClassifySettings &settings,
    double *confidence_out = nullptr );


  /** Precomputed candidate peaks for each DataSrcInfo, to avoid recomputing during GA. */
  struct PrecomputedCandidates
  {
    PeakFitUtils::CoarseResolutionType truth;
    std::vector<PeakDef> src_candidates;  // from source spectrum
    std::vector<PeakDef> bg_candidates;   // from short_background (empty if none)
    bool has_background = false;
  };

  /** Precompute candidate peaks for all input sources. */
  std::vector<PrecomputedCandidates> precompute_candidates(
    const std::vector<DataSrcInfo> &input_srcs,
    const FindCandidateSettings &candidate_settings );


  /** Evaluate Stage 1 (HPGe vs NaI) across all precomputed spectra. */
  double eval_highres_precomputed( const HighResClassifySettings &settings,
                                   const std::vector<PrecomputedCandidates> &precomputed );

  /** Print Stage 1 (HPGe vs NotHPGe) confusion matrix and per-detector results. */
  void print_stage1_results_precomputed(
    const std::vector<DataSrcInfo> &input_srcs,
    const std::vector<PrecomputedCandidates> &precomputed,
    const HighResClassifySettings &highres_settings );


  /** Print detailed per-detector diagnostic output for the classifier. */
  void print_diagnostic_output_precomputed(
    const std::vector<DataSrcInfo> &input_srcs,
    const std::vector<PrecomputedCandidates> &precomputed,
    const HighResClassifySettings &highres_settings );


  // --- Stage 1 GA infrastructure (HPGe vs NaI) ---

  struct HighResSolution
  {
    double hpge_fwhm_bias_mult;
    double nai_fwhm_bias_mult;
    double hpge_distance_weight;
    double nai_distance_weight;
    double amp_clamp_denom;
    double weight_clamp_max;
    double min_peak_energy;
    double max_peak_energy;
    double unknown_threshold;
    int min_peaks_for_classify;
    double min_peak_significance;
    double sig_fraction_of_max;
    double narrow_penalty_mult;
    double narrow_sig_threshold;
    double best_peak_czt_penalty;
    double best_peak_conf_threshold;

    std::string to_string( const std::string &separator ) const;
  };//struct HighResSolution

  struct HighResCost
  {
    double objective1;
  };//struct HighResCost


  void highres_init_genes( HighResSolution &p, const std::function<double(void)> &rnd01 );
  HighResClassifySettings highres_genes_to_settings( const HighResSolution &p );
  bool highres_eval_solution( const HighResSolution &p, HighResCost &c );

  HighResSolution highres_mutate( const HighResSolution &X_base,
                                  const std::function<double(void)> &rnd01,
                                  double shrink_scale );

  HighResSolution highres_crossover( const HighResSolution &X1,
                                     const HighResSolution &X2,
                                     const std::function<double(void)> &rnd01 );

  /** Run Stage 1 GA optimization (HPGe vs NaI). */
  HighResClassifySettings do_highres_ga_eval(
    std::function<double( const HighResClassifySettings & )> eval_fcn );

}//namespace ClassifyDetType_GA

#endif //ClassifyDetType_GA_h
