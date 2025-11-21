#ifndef InitialFit_GA_GA_h
#define InitialFit_GA_GA_h
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
#include <functional>

#define OPENGA_EXTERN_LOCAL_VARS 1

#include "openGA.hpp"

#include "InterSpec/PeakDef.h"

#include "PeakFitImprove.h"

namespace SpecUtils
{
  class Measurement;
}



namespace InitialFit_GA
{
std::vector<PeakDef> initial_peak_find_and_fit( const InitialPeakFindSettings &fit_settings,
                              const FindCandidateSettings &candidate_settings,
                              const std::shared_ptr<const SpecUtils::Measurement> &data,
                                               const bool multithread,
                              size_t &num_add_candidates_fit_for,  //Only for eval purposes
                              size_t &num_add_candidates_accepted //Only for eval purposes
);


struct PeakFindAndFitWeights
{
  double find_weight = std::numeric_limits<double>::max();

  // For def wanted peaks, this is sum of fabs(fit_area - actual_area)/max(sqrt(actual_area),0.01*actual_area),
  //  divided by number of peaks. (if actual_area is really large, then sigma is capped as though the uncert of
  //  10k counts)
  double def_wanted_area_weight = std::numeric_limits<double>::max();
  double maybe_wanted_area_weight = std::numeric_limits<double>::max();

  /** The sum of sqrt of peak areas, that are less than 10 sigma significant (maybe we messed up and fifnt expect an actual real peak) */
  double not_wanted_area_weight = std::numeric_limits<double>::max();

  double def_wanted_area_median_weight = std::numeric_limits<double>::max();
  double maybe_wanted_area_median_weight = std::numeric_limits<double>::max();
  double not_wanted_area_median_weight = std::numeric_limits<double>::max();
};//struct PeakFindAndFitWeights


PeakFindAndFitWeights eval_initial_peak_find_and_fit( const InitialPeakFindSettings &fit_settings,
                                                     const FindCandidateSettings &candidate_settings,
                                                     const DataSrcInfo &src_info,
                                                     const bool multithread );

  // I think we could actually just use `InitialPeakFindSettings` instead of defining this struck - but I'll wait on that until we get it up and going a bit
  struct InitialFitSolution
  {
    double initial_stat_threshold;
    double initial_hypothesis_threshold;
    double initial_min_nsigma_roi;
    double initial_max_nsigma_roi;
    int fwhm_fcn_form;
    double search_roi_nsigma_deficit;
    //double search_stat_threshold;
    double search_hypothesis_threshold;
    double search_stat_significance;
    double ROI_add_nsigma_required;
    double ROI_add_chi2dof_improve;

    std::string to_string( const std::string &separator ) const;
  };//struct InitialFitSolution

  InitialPeakFindSettings genes_to_settings( const InitialFitSolution &solution );

  struct InitialFitCost
  {
    // This is where the results of simulation
    // is stored but not yet finalized.
    double objective1;
  };

  typedef EA::Genetic<InitialFitSolution,InitialFitCost> GA_Type;
  typedef EA::GenerationType<InitialFitSolution,InitialFitCost> Generation_Type;

  void init_genes( InitialFitSolution& p,const std::function<double(void)> &rnd01 );

  bool eval_solution( const InitialFitSolution &p, InitialFitCost &c );

  InitialFitSolution mutate( const InitialFitSolution& X_base,
                            const std::function<double(void)> &rnd01,
                            double shrink_scale);

  InitialFitSolution crossover( const InitialFitSolution& X1,
                               const InitialFitSolution& X2,
                               const std::function<double(void)> &rnd01 );

  double calculate_SO_total_fitness( const GA_Type::thisChromosomeType &X);


  void SO_report_generation( int generation_number,
                            const EA::GenerationType<InitialFitSolution,InitialFitCost> &last_generation,
                            const InitialFitSolution& best_genes );

  InitialPeakFindSettings do_ga_eval( std::function<double( const InitialPeakFindSettings &)> ga_eval_fcn );
}//namespace InitialFit_GA

#endif //InitialFit_GA_GA_h
