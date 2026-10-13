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
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef RelActCalcAuto_ActiveSet_hpp
#define RelActCalcAuto_ActiveSet_hpp

#include <cmath>
#include <atomic>
#include <limits>
#include <memory>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <functional>

#include "Eigen/Dense"

#include "ceres/ceres.h"

#include "InterSpec/RelActCalcAuto_ProfileFit.h"

/** Active-set handling of parameters that sit on a bound, for Ceres' bound-constrained L-M.

 Why this exists: Ceres handles box constraints by computing the Levenberg-Marquardt step as if
 unconstrained, predicting the cost reduction of THAT step, and then clamping the trial point onto
 the box (`ParameterBlock::Plus`).  A step is accepted only if the reduction actually achieved by the
 clamped point is at least `min_relative_decrease` of the reduction predicted for the unclamped one.
 When a free parameter sits on its bound and its gradient pushes further outside, the part of the
 predicted reduction carried by that parameter can never be realized.  If that phantom part
 dominates (in Jacobi-scaled units, the ratio of realizable to predicted reduction stays below
 `min_relative_decrease` at every trust-region radius), every step is rejected, the radius collapses,
 and `FunctionToleranceReached()` - which Ceres evaluates on the REJECTED trial - reports
 CONVERGENCE at a point that is not a minimum.  Warm restarts cannot escape, since the ratio does not
 depend on the radius.

 This bites whenever a parameter is pinned against a bound while other parameters still have work to
 do; e.g., adding a nuclide whose true activity is zero (a nested model) drives its activity onto
 the zero bound early, and the subsequent solve then stalls with the efficiency correction never
 leaving its starting value - a far higher chi2 than the same fit without that nuclide.

 The remedy is the textbook active-set step: hold the blocked parameters constant (so the solver's
 step and its predicted reduction live entirely in the subspace that can move), re-solve, and then let
 any held parameter whose gradient has turned inward move again.
 */
namespace RelActCalcAutoImp
{
namespace ActiveSet
{

/** Ambient indices of the free parameters whose descent is blocked by a bound.

 A parameter is blocked when its gradient points out of the box at a real bound, and the room left
 to that bound is worth no more than `max_cost_fraction` of the cost to first order - so a parameter
 sitting a hair above its bound, which stalls the solver the same way, is included too.  Ceres clamps
 onto bounds exactly, so the usual case is zero room.

 @param x                 Full (ambient) parameter vector.
 @param tangent_gradient  Gradient as `ceres::Problem::Evaluate` returns it under a `SubsetManifold`:
                          one entry per free parameter, in increasing ambient order.
 @param free_indices      Ambient index of each tangent entry; see
                          `ProfileLikelihood::tangent_to_ambient_indices`.
 @param box               Ambient-length (lower, upper) bounds, using Ceres' +-max sentinels for "none".
 @param cost              The cost the gradient is of (any positive scale; only used for the tolerance).
 @param max_cost_fraction Room tolerance, as a fraction of `cost`.
 */
inline std::vector<int> bound_blocked_parameters( const std::vector<double> &x,
                                                  const std::vector<double> &tangent_gradient,
                                                  const std::vector<std::size_t> &free_indices,
                                                  const std::vector<std::pair<double,double>> &box,
                                                  const double cost,
                                                  const double max_cost_fraction )
{
  std::vector<int> blocked;
  if( tangent_gradient.size() != free_indices.size() )
    return blocked;

  const double max_first_order_change = (std::max)( 0.0, max_cost_fraction*std::fabs(cost) );

  for( std::size_t j = 0; j < free_indices.size(); ++j )
  {
    const std::size_t i = free_indices[j];
    const double g = tangent_gradient[j];
    if( (i >= x.size()) || (i >= box.size()) || !std::isfinite(g) || (g == 0.0) )
      continue;

    // Descent moves along -g: a positive gradient pushes toward the lower bound, a negative one
    //  toward the upper.
    double room = std::numeric_limits<double>::infinity();
    if( (g > 0.0) && ProfileLikelihood::has_real_bound(box[i].first) )
      room = x[i] - box[i].first;
    else if( (g < 0.0) && ProfileLikelihood::has_real_bound(box[i].second) )
      room = box[i].second - x[i];

    if( !std::isfinite(room) )
      continue;

    if( (std::max)(0.0, room)*std::fabs(g) <= max_first_order_change )
      blocked.push_back( static_cast<int>(i) );
  }//for( loop over free parameters )

  return blocked;
}//bound_blocked_parameters(...)


struct PolishSettings
{
  /** Maximum number of solves (each either holding the currently-blocked set, or plain). */
  std::size_t max_rounds = 6;
  /** See `bound_blocked_parameters`; 10x the main solve's function tolerance. */
  double max_cost_fraction = 1.0e-8;
  /** A round is kept only if it lowers the objective by more than both of these, so a fit that was
   already at its constrained minimum comes back bit-identical rather than nudged by round-off. */
  double min_relative_gain = 1.0e-7;
  double min_absolute_gain = 1.0e-2;
};//struct PolishSettings


struct PolishResult
{
  /** Number of `ceres::Solve` calls made. */
  std::size_t rounds = 0;
  /** Number of those whose result was kept. */
  std::size_t rounds_adopted = 0;
  /** Total solver iterations over all rounds. */
  std::size_t iterations = 0;
  double objective_before = std::numeric_limits<double>::quiet_NaN();
  double objective_after = std::numeric_limits<double>::quiet_NaN();
  /** Parameters found blocked at entry (ambient indices). */
  std::vector<int> first_blocked;
  /** How the last round's solve ended, and the steps it accepted (Ceres counts the initial
   evaluation as one) and rejected; used to tell a confirmed constrained optimum from a stall. */
  ceres::TerminationType last_termination = ceres::CONVERGENCE;
  int last_successful_steps = 0;
  int last_unsuccessful_steps = 0;
  /** A kept round ran out of its iteration/time budget: the point is better, but not shown stationary. */
  bool budget_exhausted = false;
  /** Set when the run was cancelled; `x` is then left as it was before the cancelled round, and
   `cancel_summary` is that round's summary with a USER_FAILURE termination. */
  bool cancelled = false;
  ceres::Solver::Summary cancel_summary;

  /** True when the polish held the blocked parameters, re-solved the rest properly, and found nothing
   more to gain: the starting point is then a constrained optimum, even if Ceres rejected steps
   against the bound.  A held round that failed, or itself only rejected steps, confirms nothing. */
  bool confirmed_constrained_optimum() const
  {
    return (rounds > 0) && (rounds_adopted == 0) && !first_blocked.empty() && !cancelled
           && ((last_termination == ceres::CONVERGENCE) || (last_termination == ceres::USER_SUCCESS))
           && ((last_successful_steps > 1) || (last_unsuccessful_steps == 0));
  }
};//struct PolishResult


/** Re-solves `problem` with its bound-blocked parameters held, until no blocked parameter is left or
 nothing more is gained.

 On entry the problem's manifold must be exactly `base_constants` (or none), and `x.data()` must be
 the problem's one registered parameter block; `x` is only ever written in place.  The manifold is
 restored to `base_constants` on every exit path.  If nothing is blocked at entry, no solve is made
 and `x` is untouched.  A round that does not lower `objective` by more than the settings' gain floor
 is undone and ends the polish, so the result is never worse than the input.

 @param objective  The full objective at a parameter vector (e.g., the cost functor's chi2); used to
                   rank rounds, independent of what Ceres reports for its iterations.
 */
inline PolishResult polish_bound_blocked( ceres::Problem &problem,
                                          std::vector<double> &x,
                                          std::vector<int> base_constants,
                                          const std::vector<std::pair<double,double>> &box,
                                          const ceres::Solver::Options &options,
                                          const std::function<double(const std::vector<double> &)> &objective,
                                          const std::shared_ptr<std::atomic_bool> &cancel,
                                          const PolishSettings &settings = PolishSettings() )
{
  PolishResult result;

  const std::size_t num_pars = x.size();
  double * const block = x.data();
  if( !num_pars || (box.size() != num_pars) || !problem.HasParameterBlock(block) )
    return result;

  std::sort( begin(base_constants), end(base_constants) );
  base_constants.erase( std::unique(begin(base_constants), end(base_constants)), end(base_constants) );

  const std::vector<std::size_t> free_indices
                        = ProfileLikelihood::tangent_to_ambient_indices( num_pars, base_constants );

  // Every round, and every way out of this function, leaves the problem with the caller's manifold.
  struct ManifoldRestorer
  {
    ceres::Problem &problem;
    double * const block;
    const std::vector<int> &base_constants;
    bool changed = false;

    void restore()
    {
      if( !changed )
        return;
      const int n = static_cast<int>( problem.ParameterBlockSize(block) );
      problem.SetManifold( block, base_constants.empty() ? nullptr
                                  : new ceres::SubsetManifold(n, base_constants) );
      changed = false;
    }

    ~ManifoldRestorer(){ restore(); }
  };//struct ManifoldRestorer

  ManifoldRestorer restorer{ problem, block, base_constants };

  const auto safe_objective = [&objective]( const std::vector<double> &pars ) -> double {
    try
    {
      const double value = objective( pars );
      return std::isfinite(value) ? value : std::numeric_limits<double>::infinity();
    }catch( std::exception & )
    {
    }
    return std::numeric_limits<double>::infinity();
  };

  double best = safe_objective( x );
  result.objective_before = best;
  result.objective_after = best;
  if( !std::isfinite(best) )
    return result;

  ceres::Problem::EvaluateOptions evaluate_options;
  evaluate_options.num_threads = (std::max)( 1, options.num_threads );

  bool last_round_was_plain = false;
  std::vector<int> last_held_blocked;
  for( std::size_t round = 0; round < settings.max_rounds; ++round )
  {
    if( cancel && cancel->load() )
      break;

    // Gradient under the caller's manifold, i.e. over exactly the parameters the solve may move.
    double cost = 0.0;
    std::vector<double> gradient;
    if( !problem.Evaluate( evaluate_options, &cost, nullptr, &gradient, nullptr ) )
      break;

    const std::vector<int> blocked = bound_blocked_parameters( x, gradient, free_indices, box, cost,
                                                               settings.max_cost_fraction );
    if( round == 0 )
      result.first_blocked = blocked;

    // Nothing blocked: on entry that means Ceres' answer stands as-is.  Later, it means the held
    //  parameters' gradients have all turned inward (or off their bound), so one plain solve lets them
    //  move - unless the previous round already was that plain solve.
    if( blocked.empty() && ((round == 0) || last_round_was_plain) )
      break;

    // Holding the same set again would just re-confirm the previous round's converged point.
    if( (round > 0) && !blocked.empty() && (blocked == last_held_blocked) )
      break;

    std::vector<int> held = base_constants;
    held.insert( end(held), begin(blocked), end(blocked) );
    std::sort( begin(held), end(held) );
    held.erase( std::unique(begin(held), end(held)), end(held) );
    if( held.size() >= num_pars )
      break; //A SubsetManifold may not fix every coordinate, and there would be nothing to solve

    const std::vector<double> x_before = x;
    if( !blocked.empty() )
    {
      problem.SetManifold( block, new ceres::SubsetManifold(static_cast<int>(num_pars), held) );
      restorer.changed = true;
    }

    ceres::Solver::Summary summary;
    ceres::Solve( options, &problem, &summary );
    restorer.restore();

    result.rounds += 1;
    result.iterations += summary.iterations.size();
    result.last_termination = summary.termination_type;
    result.last_successful_steps = summary.num_successful_steps;
    result.last_unsuccessful_steps = summary.num_unsuccessful_steps;
    last_round_was_plain = blocked.empty();
    last_held_blocked = blocked;

    if( (summary.termination_type == ceres::USER_FAILURE) || (cancel && cancel->load()) )
    {
      // An interrupted round is not checked, so its point is not kept; the caller must see a cancel.
      std::copy( begin(x_before), end(x_before), begin(x) );
      result.cancelled = true;
      result.cancel_summary = summary;
      result.cancel_summary.termination_type = ceres::USER_FAILURE;
      break;
    }

    const bool usable = (summary.termination_type == ceres::CONVERGENCE)
                        || (summary.termination_type == ceres::USER_SUCCESS)
                        || (summary.termination_type == ceres::NO_CONVERGENCE);
    const double trial = safe_objective( x );
    const double gain_floor = (std::max)( settings.min_absolute_gain, settings.min_relative_gain*best );
    if( !usable || !((best - trial) > gain_floor) )
    {
      std::copy( begin(x_before), end(x_before), begin(x) );
      break;
    }

    best = trial;
    result.objective_after = best;
    result.rounds_adopted += 1;

    // A budget-exhausted round is kept (it is the best point reached) but not built on further, the
    //  same rule the pre-fit stages use; the caller reports it.
    if( summary.termination_type == ceres::NO_CONVERGENCE )
    {
      result.budget_exhausted = true;
      break;
    }
  }//for( loop over rounds )

  return result;
}//polish_bound_blocked(...)

}//namespace ActiveSet


/** How much a parameter sitting at a limit matters to a derived quantity's uncertainty, and folding
 parameters held out of a fit's rank analysis back into its reported covariance.
 */
namespace BoundEffects
{

/** Quadratic form `g^T C g`; NaN if the sizes do not match. */
inline double quadratic_form( const std::vector<std::vector<double>> &cov, const std::vector<double> &g )
{
  const std::size_t n = g.size();
  if( cov.size() != n )
    return std::numeric_limits<double>::quiet_NaN();

  double answer = 0.0;
  for( std::size_t i = 0; i < n; ++i )
  {
    if( cov[i].size() != n )
      return std::numeric_limits<double>::quiet_NaN();
    if( g[i] == 0.0 )
      continue;
    for( std::size_t j = 0; j < n; ++j )
      answer += g[i]*cov[i][j]*g[j];
  }
  return answer;
}//quadratic_form(...)


/** For parameter `p` inside the covariance: the fractional increase in the quantity's uncertainty that
 letting `p` vary causes, `sqrt(var/var_held) - 1`, where `var_held` is the variance with `p` held at its
 value (a Schur complement).  NaN if it cannot be determined (no variance, or `p` has none); +infinity
 if `p` carries essentially all of the quantity's variance. */
inline double uncertainty_added_by_freedom( const std::vector<std::vector<double>> &cov,
                                            const std::vector<double> &g,
                                            const std::size_t p )
{
  const double var = quadratic_form( cov, g );
  if( (p >= g.size()) || !(var > 0.0) || !std::isfinite(var) || !(cov[p][p] > 0.0) )
    return std::numeric_limits<double>::quiet_NaN();

  double cross = 0.0;
  for( std::size_t k = 0; k < g.size(); ++k )
    cross += cov[p][k]*g[k];

  const double var_held = var - cross*cross/cov[p][p];
  if( !(var_held > 1.0e-12*var) )
    return std::numeric_limits<double>::infinity();

  return std::sqrt( var/var_held ) - 1.0;
}//uncertainty_added_by_freedom(...)


/** The covariance with held-out Jacobian columns appended: the block inverse of the information matrix
 `[[J^T J, J^T A], [A^T J, A^T A]]`, given `cov = (J^T J)^+` already.

 With `W = cov J^T A` and the Schur complement `S = A^T A - A^T J W`, the result is
 `cov11 = cov + W S^+ W^T`, `cov12 = -W S^+`, `cov22 = S^+`.  `S` is inverted in the basis equilibrated
 by the held columns' norms (each diagonal entry is then the fraction of that column the kept ones
 cannot explain), dropping directions below `min_unexplained` - a held column the kept ones fully
 explain gets no variance, exactly as the main covariance treats a degenerate direction.
 */
struct FoldedCovariance
{
  Eigen::MatrixXd cov11, cov12, cov22;
};

inline FoldedCovariance fold_in_held_columns( const Eigen::MatrixXd &cov, const Eigen::MatrixXd &J,
                                              const Eigen::MatrixXd &A, const double min_unexplained = 1.0e-12 )
{
  const Eigen::Index nk = A.cols();
  const Eigen::MatrixXd jt_A = J.transpose() * A;
  const Eigen::MatrixXd W = cov * jt_A;
  const Eigen::MatrixXd S = A.transpose() * A - jt_A.transpose() * W;

  Eigen::VectorXd d( nk );
  for( Eigen::Index k = 0; k < nk; ++k )
  {
    const double norm = A.col(k).norm();
    d(k) = (norm > 0.0) ? (1.0 / norm) : 0.0;
  }

  const Eigen::MatrixXd S_eq = d.asDiagonal() * S * d.asDiagonal();
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig( 0.5*(S_eq + S_eq.transpose()) );
  Eigen::VectorXd inv_lambda( nk );
  for( Eigen::Index k = 0; k < nk; ++k )
  {
    const double lambda = eig.eigenvalues()(k);
    inv_lambda(k) = (std::isfinite(lambda) && (lambda > min_unexplained)) ? (1.0 / lambda) : 0.0;
  }

  const Eigen::MatrixXd S_pinv = d.asDiagonal() * eig.eigenvectors() * inv_lambda.asDiagonal()
                                 * eig.eigenvectors().transpose() * d.asDiagonal();

  FoldedCovariance answer;
  answer.cov12 = -W * S_pinv;
  answer.cov11 = cov - answer.cov12 * W.transpose();
  answer.cov22 = S_pinv;
  return answer;
}//fold_in_held_columns(...)

}//namespace BoundEffects
}//namespace RelActCalcAutoImp

#endif //RelActCalcAuto_ActiveSet_hpp
