#ifndef FitPeaksForNuclideDev_h
#define FitPeaksForNuclideDev_h
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

#include <vector>

#include "InterSpec/PeakDef.h"

#include "PeakFitImprove.h"
#include "InitialFit_GA.h"
#include "FinalFit_GA.h"
#include "CandidatePeak_GA.h"

namespace SpecUtils
{
  class Measurement;
}


/** Combined score for evaluating peak fitting quality across multiple metrics.
 Used by both FitPeaksForNuclideDev and NuclideConfig_GA to score GA solutions.

 Combines three scoring components:
 - FinalFitScore: Peak area, width, position accuracy
 - PeakFindAndFitWeights: Peak detection completeness and area accuracy
 - CandidatePeakScore: Candidate peak detection rates
 - final_weight: Combined score (simple sum of find_weight + total_weight + score)
 */
struct CombinedPeakFitScore
{
  FinalFit_GA::FinalFitScore final_fit_score;
  InitialFit_GA::PeakFindAndFitWeights initial_fit_weights;
  CandidatePeak_GA::CandidatePeakScore candidate_peak_score;

  /** Combined weight score.
   * find_weight: higher is better (more correct peaks found)
   * total_weight: lower is better (more accurate areas)
   */
  double final_weight = 0.0;

  std::string print( const std::string &varname ) const
  {
    std::string answer;

    answer += "=== " + varname + " - Final Fit Metrics ===\n";
    answer += final_fit_score.print( varname + ".final_fit_score" );

    answer += "\n=== " + varname + " - Peak Finding Metrics ===\n";
    answer += varname + ".initial_fit_weights.find_weight =                     " + std::to_string( initial_fit_weights.find_weight ) + "\n";
    answer += varname + ".initial_fit_weights.def_wanted_area_weight =          " + std::to_string( initial_fit_weights.def_wanted_area_weight ) + "\n";
    answer += varname + ".initial_fit_weights.maybe_wanted_area_weight =        " + std::to_string( initial_fit_weights.maybe_wanted_area_weight ) + "\n";
    answer += varname + ".initial_fit_weights.not_wanted_area_weight =          " + std::to_string( initial_fit_weights.not_wanted_area_weight ) + "\n";
    answer += varname + ".initial_fit_weights.def_wanted_area_median_weight =   " + std::to_string( initial_fit_weights.def_wanted_area_median_weight ) + "\n";
    answer += varname + ".initial_fit_weights.maybe_wanted_area_median_weight = " + std::to_string( initial_fit_weights.maybe_wanted_area_median_weight ) + "\n";
    answer += varname + ".initial_fit_weights.not_wanted_area_median_weight =   " + std::to_string( initial_fit_weights.not_wanted_area_median_weight ) + "\n";

    answer += "\n=== " + varname + " - Candidate Peak Metrics ===\n";
    answer += varname + ".candidate_peak_score.score =                           " + std::to_string( candidate_peak_score.score ) + "\n";
    answer += varname + ".candidate_peak_score.num_peaks_found =                 " + std::to_string( candidate_peak_score.num_peaks_found ) + "\n";
    answer += varname + ".candidate_peak_score.num_def_wanted_not_found =         " + std::to_string( candidate_peak_score.num_def_wanted_not_found ) + "\n";
    answer += varname + ".candidate_peak_score.num_def_wanted_peaks_found =      " + std::to_string( candidate_peak_score.num_def_wanted_peaks_found ) + "\n";
    answer += varname + ".candidate_peak_score.num_possibly_accepted_peaks_not_found = " + std::to_string( candidate_peak_score.num_possibly_accepted_peaks_not_found ) + "\n";
    answer += varname + ".candidate_peak_score.num_extra_peaks =                " + std::to_string( candidate_peak_score.num_extra_peaks ) + "\n";

    answer += "\n=== " + varname + " - Combined Metrics ===\n";
    answer += varname + ".final_weight = " + std::to_string( final_weight ) + "\n";

    return answer;
  }//print(...)
};//struct CombinedPeakFitScore


namespace FitPeaksForNuclideDev
{
  void eval_peaks_for_nuclide( const std::vector<DataSrcInfo> &srcs_info );
}//namespace FitPeaksForNuclideDev

#endif //FitPeaksForNuclideDev_h
