#ifndef CandidatePeak_GA_h
#define CandidatePeak_GA_h
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

#include <tuple>
#include <string>
#include <functional>

#define OPENGA_EXTERN_LOCAL_VARS 1

#include "openGA.hpp"

#include "PeakFitImprove.h"

namespace SpecUtils
{
  class Measurement;
}


namespace CandidatePeak_GA
{
  std::vector<PeakDef> find_candidate_peaks( const std::shared_ptr<const SpecUtils::Measurement> data,
                                            size_t start_channel,
                                            size_t end_channel,
                                            const FindCandidateSettings &settings );


  struct CandidatePeakScore
  {
    double score;
    size_t num_peaks_found;
    size_t num_def_wanted_not_found;
    size_t num_def_wanted_peaks_found;
    size_t num_possibly_accepted_peaks_not_found;
    size_t num_extra_peaks;
    
    // Peak tuples for visualization (only populated when needed for N42 output)
    std::vector<std::tuple<float,float,float>> detected_expected;    //{mean, sigma, amplitude}
    std::vector<std::tuple<float,float,float>> detected_not_expected; //{mean, sigma, amplitude}
    
    // Expected peak categorization for visualization (only populated when needed for N42 output)
    std::vector<ExpectedPhotopeakInfo> possibly_expected_but_not_detected;
    std::vector<ExpectedPhotopeakInfo> expected_and_was_detected;
    std::vector<ExpectedPhotopeakInfo> def_expected_but_not_detected;
  };//struct CandidatePeakScore

  /** Calculates candidate peak score for a single source.
   * 
   * This function compares detected peaks against expected peaks and calculates
   * a score based on how well the detected peaks match expectations.
   * 
   * @param detected_peaks The peaks that were detected/found
   * @param expected_photopeaks The expected peaks to match against
   * @return CandidatePeakScore with score and statistics for this source
   */
  CandidatePeakScore
  calculate_candidate_peak_score_for_source( const std::vector<PeakDef> &detected_peaks,
                                             const std::vector<ExpectedPhotopeakInfo> &expected_photopeaks );

  /** Corrects CandidatePeakScore fields to exclude escape peaks.
   *
   * Escape peaks (S.E. at 511 keV below parent, D.E. at 1022 keV below parent) are
   * typically not fit and should not penalize the score. This function identifies
   * escape peaks in def_expected_but_not_detected by checking if a much larger
   * parent peak exists at +511 or +1022 keV in the expected photopeaks list.
   *
   * @param score The CandidatePeakScore to correct (modified in place)
   * @param expected_photopeaks The expected peaks to check for parent peaks
   */
  void correct_score_for_escape_peaks( CandidatePeakScore &score,
                                       const std::vector<ExpectedPhotopeakInfo> &expected_photopeaks );

  CandidatePeakScore
  eval_candidate_settings( const FindCandidateSettings settings, const std::vector<DataSrcInfo> &input_srcs, const bool write_n42 );


  struct CandidatePeakSolution
  {
    int num_smooth_side_channels;
    int smooth_polynomial_order;
    double threshold_FOM;
    double more_scrutiny_FOM_threshold_delta;
    double pos_sum_threshold_sf;
    int num_chan_fluctuate;
    double more_scrutiny_coarser_FOM_delta;
    double more_scrutiny_min_dev_from_line;
    int amp_to_apply_line_test_below;

    std::string to_string( const std::string &separator ) const;
  };//struct CandidatePeakSolution

  struct CandidatePeakCost
  {
    // This is where the results of simulation
    // is stored but not yet finalized.
    double objective1;
  };//struct CandidatePeakCost

  typedef EA::Genetic<CandidatePeakSolution,CandidatePeakCost> GA_Type;
  typedef EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> Generation_Type;

  void init_genes(CandidatePeakSolution& p,const std::function<double(void)> &rnd01);

  FindCandidateSettings genes_to_settings( const CandidatePeakSolution &p );

  bool eval_solution( const CandidatePeakSolution &p, CandidatePeakCost &c );

  CandidatePeakSolution mutate( const CandidatePeakSolution &X_base,
                               const std::function<double(void)> &rnd01,
                               double shrink_scale );

  CandidatePeakSolution crossover( const CandidatePeakSolution& X1,
                                  const CandidatePeakSolution& X2,
                                  const std::function<double(void)> &rnd01);

  double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X);


  void SO_report_generation( int generation_number,
                            const EA::GenerationType<CandidatePeakSolution,CandidatePeakCost> &last_generation,
                            const CandidatePeakSolution& best_genes);

  /**  Performs the genetic optimization for finding candidate peaks.

   Currently, you can only call this function once per program execution.
  */
  FindCandidateSettings do_ga_eval( std::function<double(const FindCandidateSettings &)> ga_eval_fcn );
}//namespace CandidatePeak_GA

#endif //CandidatePeak_GA_h
