#ifndef CeresUtils_hpp
#define CeresUtils_hpp
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

#include <limits>
#include <algorithm>

#include "ceres/solver.h"

namespace CeresUtils
{
  /** The Ceres cost (i.e., 0.5*sum of squared residuals) at the parameters `ceres::Solve(...)` left
   in the problem.

   Use this instead of `Solver::Summary::final_cost`, which is the minimum `cost` over *all*
   iterations - including rejected trial steps, for which Ceres records the candidate's cost - while
   the parameters returned are only ever the best *accepted* point.  A step that lowered the cost by
   less than `min_relative_decrease` of the predicted decrease is rejected (bounds-constrained
   problems make this common: the projected line search shortens a step without rescaling its
   predicted decrease), so `final_cost` can be below the cost of the solution actually returned.

   When the solution is not usable, Ceres restores the starting parameters, whose cost is
   `initial_cost` - except that `initial_cost` is for the start projected onto the parameter
   bounds, so it is not exact for a start outside them.  Returns +infinity if Ceres never
   evaluated the cost.
   */
  inline double cost_at_solution( const ceres::Solver::Summary &summary )
  {
    if( summary.initial_cost < 0.0 ) // still its default of -1: never evaluated
      return std::numeric_limits<double>::infinity();

    // Seeding with `initial_cost` covers the starting point for both minimizer types: the
    //  trust-region minimizer records iteration zero as a successful step, the line-search
    //  minimizer does not.
    double cost = summary.initial_cost;
    if( summary.IsSolutionUsable() )
    {
      for( const ceres::IterationSummary &iter : summary.iterations )
      {
        if( iter.step_is_successful )
          cost = (std::min)( cost, iter.cost );
      }
    }//if( summary.IsSolutionUsable() )

    return cost;
  }//double cost_at_solution( const ceres::Solver::Summary &summary )
}//namespace CeresUtils

#endif //CeresUtils_hpp
