#ifndef FinalFit_GA_h
#define FinalFit_GA_h
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
#include <functional>

#define OPENGA_EXTERN_LOCAL_VARS 1

#include "openGA.hpp"

#include "InterSpec/PeakDef.h"

#include "PeakFitImprove.h"

namespace FinalFit_GA
{

/** Takes in peaks in a single ROI, and fits them.

Note: this function assumes you have already delt with `FinalPeakFitSettings::combine_nsigma_near`
 and `FinalPeakFitSettings::combine_ROI_overlap_frac`.
 */
std::vector<PeakDef> final_peak_fit_for_roi( const std::vector<PeakDef> &pre_fit_peaks,
                                            const FinalPeakFitSettings &final_fit_settings,
                                            const DataSrcInfo &src_info );

std::vector<PeakDef> final_peak_fit( const std::vector<PeakDef> &pre_fit_peaks,
                                    const FinalPeakFitSettings &final_fit_settings,
                                    const DataSrcInfo &src_info );

struct FinalFitScore
{
  /** The number of (truth-level) sigma away from expected area.  Scaled from 0 (perfect) to 1.0 (15 sigma off (arbitrary)) */
  double area_score = 0.0;
  /** The fractional difference from expected. Scaled from 0 (perfect) to 1.0 (100% or more off) `*/
  double width_score = 0.0;
  /** The number of (truth-level) peak-sigma away from expected, scaled from 0 (perfect) to 1.0 (1.5 sigma off (arbitrary)) */
  double position_score = 0.0;

  /** The number of unexpected peaks that were ignored. */
  size_t ignored_unexpected_peaks = 0u;
  double unexpected_peaks_sum_significance = 0.0;

  /** The score for the quality of the fit, for optimization purposes.
   Currently just the area score.
  */
  double total_weight = 0.0;

  size_t num_peaks_used = 0;

  std::string print( const std::string &varname ) const;
};//struct FinalFitScore


FinalFitScore eval_final_peak_fit( const FinalPeakFitSettings &final_fit_settings,
                                  const DataSrcInfo &src_info,
                                  const std::vector<PeakDef> &intial_peaks,
                                  const bool write_n42 );



void do_final_peak_fit_ga_optimization( const FindCandidateSettings &candidate_settings,
                                       const InitialPeakFindSettings &initial_fit_settings,
                                       const std::vector<DataSrcInfo> &input_srcs );



struct FinalFitSolution
{
  double combine_nsigma_near;
  double combine_ROI_overlap_frac;
  double cont_type_peak_nsigma_threshold;
  double cont_type_left_right_nsigma;
  //double stepped_roi_extent_lower_side_stat_multiple;
  //double stepped_roi_extent_upper_side_stat_multiple;
  double cont_poly_order_increase_chi2dof_required;
  double cont_step_type_increase_chi2dof_required;
  double skew_nsigma;
  double left_residual_sum_min_to_try_skew;
  double right_residual_sum_min_to_try_skew;
  double skew_improve_chi2_dof_threshold;
  double roi_extent_low_num_fwhm_base;
  double roi_extent_high_num_fwhm_base;
  int roi_extent_mult_type;
  double roi_extent_lower_side_stat_multiple;
  double roi_extent_upper_side_stat_multiple;
  double multi_roi_extent_lower_side_fwhm_mult;
  double multi_roi_extent_upper_side_fwhm_mult;
  
  std::string to_string( const std::string &separator ) const;
};//struct FinalFitSolution

struct FinalFitCost
{
  // This is where the results of simulation
  // is stored but not yet finalized.
  double objective1;
};

typedef EA::Genetic<FinalFitSolution,FinalFitCost> GA_Type;
typedef EA::GenerationType<FinalFitSolution,FinalFitCost> Generation_Type;

void init_genes(FinalFitSolution& p,const std::function<double(void)> &rnd01);

FinalPeakFitSettings genes_to_settings( const FinalFitSolution &solution );

bool eval_solution( const FinalFitSolution &p, FinalFitCost &c );

FinalFitSolution mutate(
                        const FinalFitSolution& X_base,
                        const std::function<double(void)> &rnd01,
                        double shrink_scale);

FinalFitSolution crossover(
                           const FinalFitSolution& X1,
                           const FinalFitSolution& X2,
                           const std::function<double(void)> &rnd01);

double calculate_SO_total_fitness(const GA_Type::thisChromosomeType &X);

void SO_report_generation(
                          int generation_number,
                          const EA::GenerationType<FinalFitSolution,FinalFitCost> &last_generation,
                          const FinalFitSolution& best_genes);

FinalPeakFitSettings do_ga_eval( std::function<double(const FinalPeakFitSettings &)> ga_eval_fcn );
};//namespace FinalFit_GA

#endif //FinalFit_GA_h
