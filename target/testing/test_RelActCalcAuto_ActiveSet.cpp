/* Tests for the active-set handling of bound-pinned parameters used by RelActCalcAuto's solve. */
#include "InterSpec_config.h"

#include <set>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <iostream>

#define BOOST_TEST_MODULE RelActCalcAuto_ActiveSet_suite
#include <boost/test/included/unit_test.hpp>

#include "Eigen/Dense"

#include "ceres/ceres.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/RelActCalcAuto_ActiveSet.hpp"

using namespace std;
using namespace boost::unit_test;
using namespace RelActCalcAutoImp::ActiveSet;

namespace
{
const double sm_max = std::numeric_limits<double>::max();
const double sm_lowest = std::numeric_limits<double>::lowest();

/** p0 >= 0 is pushed hard against its bound (it wants to be -0.1), while (p1,p2) is a Rosenbrock
 valley started away from its minimum at (1,1).  The constrained optimum is (0,1,1), with a sum of
 squared residuals of exactly 1e6.  The two parts do not interact, so the only thing linking them is
 Ceres' step-acceptance test: p0 carries nearly all the predicted reduction, none of which a clamped
 step can realize. */
struct BoundDominatedResidual
{
  template<typename T>
  bool operator()( const T * const p, T *r ) const
  {
    r[0] = T(1.0e4) * (p[0] + T(0.1));
    r[1] = T(100.0) * (p[2] - p[1]*p[1]);
    r[2] = T(10.0) * (T(1.0) - p[1]);
    return true;
  }
};//struct BoundDominatedResidual


double sum_sq( const std::vector<double> &p )
{
  const double r0 = 1.0e4 * (p[0] + 0.1);
  const double r1 = 100.0 * (p[2] - p[1]*p[1]);
  const double r2 = 10.0 * (1.0 - p[1]);
  return r0*r0 + r1*r1 + r2*r2;
}


/** The options `RelActAutoCostFcn::solve_ceres` uses for its main solve. */
ceres::Solver::Options relact_like_options()
{
  ceres::Solver::Options options;
  options.minimizer_type = ceres::TRUST_REGION;
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  options.use_nonmonotonic_steps = true;
  options.max_consecutive_nonmonotonic_steps = 5;
  options.max_num_consecutive_invalid_steps = 10;
  options.initial_trust_region_radius = 1.0e4;
  options.max_trust_region_radius = 1.0e16;
  options.min_trust_region_radius = 1.0e-32;
  options.min_relative_decrease = 1.0e-3;
  options.linear_solver_type = ceres::DENSE_QR;
  options.max_num_iterations = 50000;
  options.function_tolerance = 1.0e-9;
  options.gradient_tolerance = 1.0e-13;
  options.parameter_tolerance = 1.0e-11;
  options.logging_type = ceres::SILENT;
  options.minimizer_progress_to_stdout = false;
  return options;
}


void setup_problem( ceres::Problem &problem, std::vector<double> &x )
{
  problem.AddResidualBlock( new ceres::AutoDiffCostFunction<BoundDominatedResidual,3,3>(
                                                                 new BoundDominatedResidual() ),
                            nullptr, x.data() );
  problem.SetParameterLowerBound( x.data(), 0, 0.0 );
}


const std::vector<std::pair<double,double>> sm_box = { {0.0, sm_max}, {sm_lowest, sm_max}, {sm_lowest, sm_max} };
}//namespace


BOOST_AUTO_TEST_CASE( blocked_parameter_classification )
{
  const std::vector<std::size_t> free_indices = { 0, 1, 2, 3 };
  const std::vector<std::pair<double,double>> box = { {0.0, 10.0}, {0.0, 10.0}, {sm_lowest, sm_max}, {-1.0, 1.0} };

  // On the lower bound pushing out; on the lower bound pulling in; unbounded; on the upper bound pushing out.
  {
    const std::vector<double> x = { 0.0, 0.0, 5.0, 1.0 };
    const std::vector<double> g = { 3.0, -3.0, 100.0, -2.0 };
    const std::vector<int> blocked = bound_blocked_parameters( x, g, free_indices, box, 1.0, 1.0e-8 );
    BOOST_CHECK( blocked == std::vector<int>({0, 3}) );
  }

  // A hair above the bound counts only while the room left is worth less than the cost tolerance.
  {
    const std::vector<double> x = { 1.0e-12, 1.0e-3, 5.0, 0.5 };
    const std::vector<double> g = { 3.0, 3.0, 0.0, -2.0 };
    const std::vector<int> blocked = bound_blocked_parameters( x, g, free_indices, box, 100.0, 1.0e-8 );
    BOOST_CHECK( blocked == std::vector<int>({0}) );
  }

  // Ceres' "no bound" sentinels are never treated as a bound, however close a value is to them.
  {
    const std::vector<double> x = { 5.0, 5.0, sm_max, 0.0 };
    const std::vector<double> g = { 0.0, 0.0, -1.0, 0.0 };
    BOOST_CHECK( bound_blocked_parameters( x, g, free_indices, box, 1.0, 1.0e-8 ).empty() );
  }

  // Tangent entries map through `free_indices` (parameter 1 held constant).
  {
    const std::vector<std::size_t> tangent_free = { 0, 2, 3 };
    const std::vector<double> x = { 3.0, 0.0, 5.0, -1.0 };
    const std::vector<double> g = { 1.0, 1.0, 2.0 };
    const std::vector<int> blocked = bound_blocked_parameters( x, g, tangent_free, box, 1.0, 1.0e-8 );
    BOOST_CHECK( blocked == std::vector<int>({3}) );
  }

  // A parameter closing on its bound under a large gradient (as Ceres' projected line search leaves it)
  //  is worth more than the cost tolerance, but is within the range fraction of its bound.
  {
    const std::vector<std::size_t> free = { 0, 1, 2 };
    const std::vector<std::pair<double,double>> range_box = { {-0.3, 0.3}, {0.0, sm_max}, {-0.3, 0.3} };
    const std::vector<double> x = { -0.3 + 2.3e-10, 2.3e-10, -0.3 + 1.0e-3 };
    const std::vector<double> g = { 5.3e5, 5.3e5, 5.3e5 };
    BOOST_CHECK( bound_blocked_parameters( x, g, free, range_box, 1578.0, 1.0e-8 ).empty() );
    // Only the parameter with both bounds, and within the fraction of its range, counts.
    const std::vector<int> blocked = bound_blocked_parameters( x, g, free, range_box, 1578.0, 1.0e-8, 1.0e-6 );
    BOOST_CHECK( blocked == std::vector<int>({0}) );
    // Pulled away from the bound, it is not blocked however close it is.
    const std::vector<double> inward = { -5.3e5, 5.3e5, 5.3e5 };
    BOOST_CHECK( bound_blocked_parameters( x, inward, free, range_box, 1578.0, 1.0e-8, 1.0e-6 ).empty() );
  }

  // A loose limit sets no scale: an activity-like parameter allowed up to 1E9 is pressed against its
  //  lower bound of 1 only within 1E-6 of it - not anywhere within 1E-6 of its range - but against its
  //  upper bound within 1E-6 of that bound's value.
  {
    const std::vector<std::size_t> free = { 0, 1, 2 };
    const std::vector<std::pair<double,double>> loose_box = { {1.0, 1.0e9}, {1.0, 1.0e9}, {1.0, 1.0e9} };
    const std::vector<double> x = { 100.0, 1.0 + 5.0e-7, 1.0e9 - 500.0 };
    const std::vector<double> g = { 1.0, 1.0, -1.0 };
    const std::vector<int> blocked = bound_blocked_parameters( x, g, free, loose_box, 1.0, 1.0e-8, 1.0e-6 );
    BOOST_CHECK( blocked == std::vector<int>({1, 2}) );
  }
}//BOOST_AUTO_TEST_CASE( blocked_parameter_classification )


BOOST_AUTO_TEST_CASE( pressed_parameters_snap_onto_their_bound )
{
  const std::vector<std::size_t> free = { 0, 1, 2, 3 };
  const std::vector<std::pair<double,double>> box = { {-0.3, 0.3}, {1.0, 1.0e9}, {1.0, 1.0e9}, {0.0, sm_max} };
  const std::vector<double> start = { -0.3 + 1.0e-8, 1.0 + 5.0e-7, 100.0, 1.0e-12 };

  // Pressed two-sided parameters land exactly on their bound; one far off a loose bound, and one with
  //  only one bound, do not move.
  std::vector<double> x = start;
  snap_pressed_to_bounds( x, { 1.0, 1.0, 1.0, 1.0 }, free, box, 1.0e-6 );
  BOOST_CHECK_EQUAL( x[0], -0.3 );
  BOOST_CHECK_EQUAL( x[1], 1.0 );
  BOOST_CHECK_EQUAL( x[2], 100.0 );
  BOOST_CHECK_EQUAL( x[3], 1.0e-12 );

  // Pulled away from the bound, nothing moves.
  x = start;
  snap_pressed_to_bounds( x, { -1.0, -1.0, -1.0, -1.0 }, free, box, 1.0e-6 );
  BOOST_CHECK( x == start );
}//BOOST_AUTO_TEST_CASE( pressed_parameters_snap_onto_their_bound )


BOOST_AUTO_TEST_CASE( bound_press_detector )
{
  // p0 two-sided, p1 one-sided (never watched), p2 two-sided but not free.
  const std::vector<std::pair<double,double>> box = { {0.0, 1.0}, {0.0, sm_max}, {0.0, 1.0} };
  const std::vector<std::size_t> free = { 0, 1 };

  // What Ceres hands an iteration callback: `objective` is twice its cost.
  const auto step = []( const int iteration, const double objective, const bool accepted ) {
    ceres::IterationSummary summary;
    summary.iteration = iteration;
    summary.cost = 0.5*objective;
    summary.step_is_successful = accepted;
    return summary;
  };

  // Fires on the `window`-th consecutive pressed iteration, when the objective barely moved over it.
  {
    BoundPressDetector detector( free, box, 4, 1.0e-6, 1.0 );
    BOOST_CHECK( detector.watching_anything() );
    const double x[] = { 1.0 - 1.0e-7, 0.5, 0.5 };
    BOOST_CHECK( !detector.iteration( x, step(0, 100.1, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(1, 100.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(2, 99.9, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(3, 99.8, true) ) );
    BOOST_CHECK( detector.iteration( x, step(4, 99.7, true) ) );
  }

  // Real progress while pressed restarts the watch; it fires once that progress dies down.
  {
    BoundPressDetector detector( free, box, 3, 1.0e-6, 1.0 );
    const double x[] = { 0.0, 0.5, 0.5 };
    BOOST_CHECK( !detector.iteration( x, step(0, 110.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(1, 100.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(2, 90.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(3, 80.0, true) ) );   // dropped 20 over the window: keep going
    BOOST_CHECK( !detector.iteration( x, step(4, 79.9, true) ) );
    BOOST_CHECK( detector.iteration( x, step(5, 79.8, true) ) );    // 0.2 over the window since 80
  }

  // The objective is the best point's, which `x` is: a rejected trial (reported at its own cost, up to
  //  Ceres' max-double sentinel) or a non-monotonic accepted step (above the best) ending a window does
  //  not look like a stall - nor, at its start, like progress.
  {
    BoundPressDetector detector( free, box, 3, 1.0e-6, 1.0 );
    const double x[] = { 0.0, 0.5, 0.5 };
    BOOST_CHECK( !detector.iteration( x, step(0, 100.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(1, 90.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(2, 80.0, true) ) );
    BOOST_CHECK( !detector.iteration( x, step(3, std::numeric_limits<double>::max(), false) ) ); // 90 -> 80
    BOOST_CHECK( !detector.iteration( x, step(4, 85.0, true) ) );    // non-monotonic: the best stays 80
    BOOST_CHECK( !detector.iteration( x, step(5, 70.0, true) ) );    // 80 -> 70 over this window
    BOOST_CHECK( !detector.iteration( x, step(6, 1.0e6, false) ) );
    BOOST_CHECK( detector.iteration( x, step(7, 75.0, true) ) );     // best still 70: stalled
  }

  // Leaving the bound resets the count; one-sided and held parameters are never watched.
  {
    BoundPressDetector detector( free, box, 2, 1.0e-6, 1.0 );
    const double pressed[] = { 1.0, 0.5, 0.5 }, away[] = { 0.5, 0.0, 0.0 };
    BOOST_CHECK( !detector.iteration( pressed, step(0, 10.0, true) ) );
    BOOST_CHECK( !detector.iteration( pressed, step(1, 10.0, true) ) );
    BOOST_CHECK( !detector.iteration( away, step(2, 10.0, true) ) );
    BOOST_CHECK( !detector.iteration( pressed, step(3, 10.0, true) ) );
    BOOST_CHECK( detector.iteration( pressed, step(4, 10.0, true) ) );

    BoundPressDetector one_sided( {1}, box, 1, 1.0e-6, 1.0 );
    BOOST_CHECK( !one_sided.watching_anything() );
    BOOST_CHECK( !one_sided.iteration( away, step(1, 10.0, true) ) );
  }
}//BOOST_AUTO_TEST_CASE( bound_press_detector )


BOOST_AUTO_TEST_CASE( polish_escapes_bound_dominated_stall )
{
  std::vector<double> x = { 0.0, -1.0, 1.0 };
  ceres::Problem problem;
  setup_problem( problem, x );

  const ceres::Solver::Options options = relact_like_options();
  ceres::Solver::Summary summary;
  ceres::Solve( options, &problem, &summary );

  // The precondition: plain Ceres reports CONVERGENCE well above the constrained optimum.  If a Ceres
  //  upgrade ever makes this fail, the stall this machinery works around may be gone - re-measure.
  const double stalled = sum_sq( x );
  BOOST_REQUIRE( summary.termination_type == ceres::CONVERGENCE );
  BOOST_REQUIRE_GT( stalled, 1.0e6 + 100.0 );

  const PolishResult result = polish_bound_blocked( problem, x, {}, sm_box, options, sum_sq, nullptr );

  BOOST_CHECK( result.first_blocked == std::vector<int>({0}) );
  BOOST_CHECK_GE( result.rounds_adopted, 1u );
  BOOST_CHECK_CLOSE_FRACTION( result.objective_before, stalled, 1.0e-12 );
  BOOST_CHECK_CLOSE_FRACTION( result.objective_after, sum_sq(x), 1.0e-12 );
  BOOST_CHECK_LT( sum_sq(x) - 1.0e6, 1.0e-3 );
  BOOST_CHECK_EQUAL( x[0], 0.0 );
  // The relative function tolerance is measured against a cost dominated by p0's residual, so the
  //  Rosenbrock part stops a little short of (1,1).
  BOOST_CHECK_CLOSE_FRACTION( x[1], 1.0, 2.0e-3 );
  BOOST_CHECK_CLOSE_FRACTION( x[2], 1.0, 2.0e-3 );

  // The caller's (absent) manifold is back: every parameter is in the tangent space again.
  BOOST_CHECK_EQUAL( problem.ParameterBlockTangentSize(x.data()), 3 );
}//BOOST_AUTO_TEST_CASE( polish_escapes_bound_dominated_stall )


BOOST_AUTO_TEST_CASE( polish_keeps_held_constants_and_manifold )
{
  // With p2 held constant by the caller, the polish must neither move it nor drop its manifold.
  std::vector<double> x = { 0.0, -1.0, 1.0 };
  ceres::Problem problem;
  setup_problem( problem, x );
  problem.SetManifold( x.data(), new ceres::SubsetManifold(3, {2}) );

  const ceres::Solver::Options options = relact_like_options();
  ceres::Solver::Summary summary;
  ceres::Solve( options, &problem, &summary );

  const double stalled = sum_sq( x );
  const PolishResult result = polish_bound_blocked( problem, x, {2}, sm_box, options, sum_sq, nullptr );
  BOOST_CHECK( result.first_blocked == std::vector<int>({0}) );
  BOOST_CHECK_EQUAL( x[2], 1.0 );
  BOOST_CHECK_EQUAL( x[0], 0.0 );
  BOOST_CHECK_LE( sum_sq(x), stalled );
  BOOST_CHECK_EQUAL( problem.ParameterBlockTangentSize(x.data()), 2 );
}//BOOST_AUTO_TEST_CASE( polish_keeps_held_constants_and_manifold )


BOOST_AUTO_TEST_CASE( polish_is_a_no_op_at_the_optimum )
{
  const ceres::Solver::Options options = relact_like_options();

  // At the constrained optimum p0 is blocked, but holding it gains nothing: bit-identical result.
  {
    std::vector<double> x = { 0.0, 1.0, 1.0 };
    ceres::Problem problem;
    setup_problem( problem, x );
    const std::vector<double> before = x;

    const PolishResult result = polish_bound_blocked( problem, x, {}, sm_box, options, sum_sq, nullptr );
    BOOST_CHECK( result.first_blocked == std::vector<int>({0}) );
    BOOST_CHECK_EQUAL( result.rounds_adopted, 0u );
    BOOST_CHECK( x == before );
  }

  // Nothing on a bound: no solve is made at all.
  {
    std::vector<double> x = { 0.5, -1.0, 1.0 };
    ceres::Problem problem;
    setup_problem( problem, x );
    const std::vector<double> before = x;

    const PolishResult result = polish_bound_blocked( problem, x, {}, sm_box, options, sum_sq, nullptr );
    BOOST_CHECK_EQUAL( result.rounds, 0u );
    BOOST_CHECK( result.first_blocked.empty() );
    BOOST_CHECK( x == before );
  }
}//BOOST_AUTO_TEST_CASE( polish_is_a_no_op_at_the_optimum )


namespace
{
/** An efficiency-like fit, exp(a + b*t + c*t^2) to Poisson data made with c = -1, where c is limited to
 [-0.3, 0.3]: c is pressed against its lower bound while a and b adjust to it; d is independent and wants
 10.  Ceres ends the plain solve with c on its bound and chi2 about 4% above the constrained optimum. */
struct PressedEfficiencyResidual
{
  static constexpr int sm_num_points = 20;

  template<typename T>
  bool operator()( const T * const p, T *r ) const
  {
    for( int i = 0; i < sm_num_points; ++i )
    {
      const double t = -1.0 + 2.0*i/(sm_num_points - 1.0);
      const double y = 1000.0*std::exp( 1.0 + 0.5*t - t*t );
      const T model = T(1000.0) * exp( p[1] + p[2]*t + p[0]*t*t );
      r[i] = (model - T(y)) / T(std::sqrt(y));
    }
    r[sm_num_points] = p[3] - T(10.0);
    return true;
  }
};//struct PressedEfficiencyResidual


const std::vector<std::pair<double,double>> sm_efficiency_box
                     = { {-0.3, 0.3}, {sm_lowest, sm_max}, {sm_lowest, sm_max}, {sm_lowest, sm_max} };


void setup_efficiency_problem( ceres::Problem &problem, std::vector<double> &x )
{
  const int nres = PressedEfficiencyResidual::sm_num_points + 1;
  problem.AddResidualBlock( new ceres::AutoDiffCostFunction<PressedEfficiencyResidual,ceres::DYNAMIC,4>(
                                                         new PressedEfficiencyResidual(), nres ),
                            nullptr, x.data() );
  problem.SetParameterLowerBound( x.data(), 0, -0.3 );
  problem.SetParameterUpperBound( x.data(), 0, 0.3 );
}


/** The chi2 of `setup_efficiency_problem`'s problem at `p`. */
double efficiency_chi2( const std::vector<double> &p )
{
  std::vector<double> r( PressedEfficiencyResidual::sm_num_points + 1 );
  PressedEfficiencyResidual()( p.data(), r.data() );
  double chi2 = 0.0;
  for( const double v : r )
    chi2 += v*v;
  return chi2;
}


/** The constrained optimum of `setup_efficiency_problem`'s problem: c held on its bound, the rest solved. */
double efficiency_constrained_optimum( const ceres::Solver::Options &options )
{
  std::vector<double> held = { -0.3, 0.0, 0.0, 0.0 };
  ceres::Problem problem;
  setup_efficiency_problem( problem, held );
  problem.SetManifold( held.data(), new ceres::SubsetManifold(4, {0}) );
  ceres::Solver::Summary summary;
  ceres::Solve( options, &problem, &summary );
  return efficiency_chi2( held );
}


/** p0 in [0,10] is tied to p1 (which wants 5); p2 is independent and wants 10. */
struct PressedResidual
{
  template<typename T>
  bool operator()( const T * const p, T *r ) const
  {
    r[0] = T(100.0) * (p[0] - p[1]);
    r[1] = p[1] - T(5.0);
    r[2] = p[2] - T(10.0);
    return true;
  }
};//struct PressedResidual


double pressed_sum_sq( const std::vector<double> &p )
{
  const double r0 = 100.0*(p[0] - p[1]), r1 = p[1] - 5.0, r2 = p[2] - 10.0;
  return r0*r0 + r1*r1 + r2*r2;
}


/** What `solve_stopping_creep` does when it may stop at a pressed bound: a solve ended by a
 `BoundPressDetector` over the parameters the manifold leaves free (bounds from `sm_efficiency_box`). */
bool solve_stopping_at_pressed_bound( const ceres::Solver::Options &options, ceres::Problem &problem,
                                      ceres::Solver::Summary &summary, const bool may_stop_early )
{
  std::vector<double *> blocks;
  problem.GetParameterBlocks( &blocks );
  BOOST_REQUIRE_EQUAL( blocks.size(), 1u );
  double * const x = blocks[0];
  const int npar = problem.ParameterBlockSize( x );

  std::vector<std::size_t> free;
  const ceres::Manifold * const manifold = problem.GetManifold( x );
  const int ntan = manifold ? manifold->TangentSize() : npar;
  std::vector<double> plus_jac( static_cast<size_t>(npar*ntan), 0.0 );
  if( manifold )
    manifold->PlusJacobian( x, plus_jac.data() );
  for( int row = 0; row < npar; ++row )
  {
    bool is_free = !manifold;
    for( int col = 0; !is_free && (col < ntan); ++col )
      is_free = (plus_jac[static_cast<size_t>(row*ntan + col)] != 0.0);
    if( is_free )
      free.push_back( static_cast<size_t>(row) );
  }

  struct Callback : public ceres::IterationCallback
  {
    BoundPressDetector detector;
    const double *x;
    bool stopped = false;
    Callback( BoundPressDetector d, const double *p ) : detector( std::move(d) ), x( p ) {}
    ceres::CallbackReturnType operator()( const ceres::IterationSummary &it ) override
    {
      if( detector.iteration( x, it ) )
      {
        stopped = true;
        return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
      }
      return ceres::SOLVER_CONTINUE;
    }
  };//struct Callback

  Callback callback( BoundPressDetector( free, sm_efficiency_box, 10, 1.0e-6, 1.0 ), x );
  ceres::Solver::Options stop_options = options;
  if( may_stop_early )
  {
    stop_options.update_state_every_iteration = true;
    stop_options.callbacks.push_back( &callback );
  }
  ceres::Solve( stop_options, &problem, &summary );
  return callback.stopped;
}//solve_stopping_at_pressed_bound(...)
}//namespace


BOOST_AUTO_TEST_CASE( pressed_bound_is_handed_to_the_polish )
{
  const ceres::Solver::Options options = relact_like_options();

  const double optimum = efficiency_constrained_optimum( options );

  // The precondition: plain Ceres stops well above it.  If a Ceres upgrade ever makes this fail, the
  //  bound handling this works around may have changed - re-measure.
  std::vector<double> plain = { 0.0, 0.0, 0.0, 0.0 };
  size_t plain_iterations = 0;
  {
    ceres::Problem problem;
    setup_efficiency_problem( problem, plain );
    ceres::Solver::Summary summary;
    ceres::Solve( options, &problem, &summary );
    plain_iterations = summary.iterations.size();
    BOOST_TEST_MESSAGE( "Plain solve: " << plain_iterations << " iterations, chi2 " << efficiency_chi2(plain)
                        << " (constrained optimum " << optimum << ")" );
    BOOST_REQUIRE_GT( efficiency_chi2(plain), optimum + 10.0 );
  }

  // Stopping once c is pressed against its bound, and polishing, reaches the optimum.
  std::vector<double> x = { 0.0, 0.0, 0.0, 0.0 };
  ceres::Problem problem;
  setup_efficiency_problem( problem, x );
  ceres::Solver::Summary summary;
  const bool stopped = solve_stopping_at_pressed_bound( options, problem, summary, true );
  BOOST_CHECK( stopped );
  BOOST_CHECK_LT( summary.iterations.size(), plain_iterations );

  const PolishResult result = polish_bound_blocked( problem, x, {}, sm_efficiency_box, options, efficiency_chi2,
                                                    nullptr, PolishSettings(), &solve_stopping_at_pressed_bound,
                                                    stopped );
  BOOST_CHECK( result.first_blocked == std::vector<int>({0}) );
  BOOST_CHECK_GE( result.rounds_adopted, 1u );
  BOOST_CHECK_LT( efficiency_chi2(x) - optimum, 1.0e-6*optimum );
  BOOST_CHECK_EQUAL( x[0], -0.3 );
  BOOST_CHECK_CLOSE_FRACTION( x[3], 10.0, 1.0e-6 );
  BOOST_CHECK_EQUAL( problem.ParameterBlockTangentSize(x.data()), 4 );
}//BOOST_AUTO_TEST_CASE( pressed_bound_is_handed_to_the_polish )


BOOST_AUTO_TEST_CASE( early_stop_with_nothing_to_hold_is_finished )
{
  // A solve that ended early, but whose parameter on the bound turns out not to be blocked (here p0 sits
  //  on its bound with zero gradient, though its optimum, following p1, is 5), is finished by a plain solve.
  const ceres::Solver::Options options = relact_like_options();
  std::vector<double> x = { 0.0, 0.0, 0.0 };   // p1 = 0 and p0 = 0 on its bound, p2 far from 10
  ceres::Problem problem;
  problem.AddResidualBlock( new ceres::AutoDiffCostFunction<PressedResidual,3,3>( new PressedResidual() ),
                            nullptr, x.data() );
  problem.SetParameterLowerBound( x.data(), 0, 0.0 );
  problem.SetParameterUpperBound( x.data(), 0, 10.0 );
  const std::vector<std::pair<double,double>> box = { {0.0, 10.0}, {sm_lowest, sm_max}, {sm_lowest, sm_max} };

  // With p0 == p1, p0's gradient is zero, so nothing is blocked.
  const PolishResult result = polish_bound_blocked( problem, x, {}, box, options, pressed_sum_sq, nullptr,
                                                    PolishSettings(), {}, true );
  BOOST_CHECK( result.first_blocked.empty() );
  BOOST_CHECK_GE( result.rounds_adopted, 1u );
  BOOST_CHECK_LT( pressed_sum_sq(x), 1.0e-6 );
  BOOST_CHECK_CLOSE_FRACTION( x[2], 10.0, 1.0e-6 );
}//BOOST_AUTO_TEST_CASE( early_stop_with_nothing_to_hold_is_finished )


BOOST_AUTO_TEST_CASE( early_stop_is_finished_when_the_next_round_fails )
{
  // An early stop whose next round is not kept (here: it too "stops early", but comes back worse, as
  //  round-off between Ceres' cost and the objective can make it) must still be finished, not returned.
  const ceres::Solver::Options options = relact_like_options();
  const double optimum = efficiency_constrained_optimum( options );

  std::vector<double> x = { 0.0, 0.0, 0.0, 0.0 };
  ceres::Problem problem;
  setup_efficiency_problem( problem, x );
  ceres::Solver::Summary summary;
  BOOST_REQUIRE( solve_stopping_at_pressed_bound( options, problem, summary, true ) );
  const double early_chi2 = efficiency_chi2( x );
  BOOST_REQUIRE_GT( early_chi2, optimum + 1.0 );

  int calls = 0;
  const auto solve = [&calls]( const ceres::Solver::Options &opts, ceres::Problem &prob,
                               ceres::Solver::Summary &summ, const bool may_stop_early ) -> bool {
    calls += 1;
    if( calls == 1 )
    {
      BOOST_CHECK( may_stop_early );
      std::vector<double *> blocks;
      prob.GetParameterBlocks( &blocks );
      // d away from its 10 by one more: chi2 up by at least 1, far more than holding c on its bound gains
      //  (d need not be near 10 yet - the creep slows every parameter).
      blocks[0][3] += (blocks[0][3] < 10.0) ? -1.0 : 1.0;
      summ.termination_type = ceres::USER_SUCCESS;
      return true;
    }
    BOOST_CHECK( !may_stop_early );   // The forced finish allows no early stop.
    return solve_stopping_at_pressed_bound( opts, prob, summ, may_stop_early );
  };

  const PolishResult result = polish_bound_blocked( problem, x, {}, sm_efficiency_box, options, efficiency_chi2,
                                                    nullptr, PolishSettings(), solve, true );
  BOOST_CHECK_EQUAL( calls, 2 );
  BOOST_CHECK_EQUAL( result.rounds_adopted, 1u );   // the forced finish only
  BOOST_CHECK( !result.stopped_short );
  BOOST_CHECK_LT( efficiency_chi2(x) - optimum, 1.0e-6*optimum );
  BOOST_CHECK_EQUAL( x[0], -0.3 );
}//BOOST_AUTO_TEST_CASE( early_stop_is_finished_when_the_next_round_fails )


BOOST_AUTO_TEST_CASE( early_stopped_rounds_are_not_counted_against_max_rounds )
{
  const ceres::Solver::Options options = relact_like_options();
  const double optimum = efficiency_constrained_optimum( options );

  // A solve that ends "early" without moving whenever it may, and otherwise runs to convergence.
  int calls = 0, early_allowed = 0;
  const auto solve = [&]( const ceres::Solver::Options &opts, ceres::Problem &prob,
                          ceres::Solver::Summary &summ, const bool may_stop_early ) -> bool {
    calls += 1;
    if( may_stop_early )
    {
      early_allowed += 1;
      summ.termination_type = ceres::USER_SUCCESS;
      return true;
    }
    return solve_stopping_at_pressed_bound( opts, prob, summ, false );
  };

  const auto early_point = [&options]( ceres::Problem &problem, std::vector<double> &x ) {
    setup_efficiency_problem( problem, x );
    ceres::Solver::Summary summary;
    BOOST_REQUIRE( solve_stopping_at_pressed_bound( options, problem, summary, true ) );
  };

  // With room for a single full round, an early-stopped round still leaves that round to finish the fit.
  {
    std::vector<double> x = { 0.0, 0.0, 0.0, 0.0 };
    ceres::Problem problem;
    early_point( problem, x );

    PolishSettings settings;
    settings.max_rounds = 1;
    const PolishResult result = polish_bound_blocked( problem, x, {}, sm_efficiency_box, options,
                                                      efficiency_chi2, nullptr, settings, solve, true );
    BOOST_CHECK_EQUAL( calls, 2 );
    BOOST_CHECK_EQUAL( early_allowed, 1 );
    BOOST_CHECK_EQUAL( result.rounds, 2u );
    BOOST_CHECK( !result.stopped_short );
    BOOST_CHECK_LT( efficiency_chi2(x) - optimum, 1.0e-6*optimum );
  }

  // With no early-stopped rounds allowed, none is offered.
  {
    calls = early_allowed = 0;
    std::vector<double> x = { 0.0, 0.0, 0.0, 0.0 };
    ceres::Problem problem;
    early_point( problem, x );

    PolishSettings settings;
    settings.max_early_stopped_rounds = 0;
    const PolishResult result = polish_bound_blocked( problem, x, {}, sm_efficiency_box, options,
                                                      efficiency_chi2, nullptr, settings, solve, true );
    BOOST_CHECK_GE( calls, 1 );
    BOOST_CHECK_EQUAL( early_allowed, 0 );
    BOOST_CHECK_LT( efficiency_chi2(x) - optimum, 1.0e-6*optimum );
  }
}//BOOST_AUTO_TEST_CASE( early_stopped_rounds_are_not_counted_against_max_rounds )



BOOST_AUTO_TEST_CASE( added_by_freedom_matches_conditional_variance )
{
  using namespace RelActCalcAutoImp::BoundEffects;

  // Inside the covariance: holding parameter 1 of x0 + x1 leaves the conditional variance of x0,
  //  4 - 2^2/3 = 8/3, against the full 11.
  {
    const vector<vector<double>> cov = { {4.0, 2.0, 0.0}, {2.0, 3.0, 0.0}, {0.0, 0.0, 1.0} };
    const vector<double> g = { 1.0, 1.0, 0.0 };
    BOOST_CHECK_CLOSE_FRACTION( quadratic_form(cov, g), 11.0, 1.0e-14 );
    BOOST_CHECK_CLOSE_FRACTION( uncertainty_added_by_freedom(cov, g, 1), std::sqrt(11.0/(8.0/3.0)) - 1.0, 1.0e-12 );
    // A parameter the quantity does not depend on, and that is uncorrelated with it, adds nothing.
    BOOST_CHECK_SMALL( uncertainty_added_by_freedom(cov, g, 2), 1.0e-14 );
    // A parameter carrying all of the variance.
    BOOST_CHECK( std::isinf( uncertainty_added_by_freedom(cov, {0.0, 0.0, 1.0}, 2) ) );
  }
}//BOOST_AUTO_TEST_CASE( added_by_freedom_matches_conditional_variance )


BOOST_AUTO_TEST_CASE( fold_in_held_columns_matches_full_inverse )
{
  using namespace RelActCalcAutoImp::BoundEffects;

  // Residual Jacobian columns j0, j1 kept, and two held columns a0, a1; folding the held ones into the
  //  kept covariance must reproduce the inverse of the full information matrix.
  Eigen::MatrixXd J( 5, 2 ), A( 5, 2 );
  J << 1, 0,
       0, 1,
       1, 1,
       1, 0,
       0, 2;
  A << 1, 0,
       1, 1,
       0, 1,
       0, 0,
       2, 1;

  const Eigen::MatrixXd cov = (J.transpose() * J).inverse();
  const FoldedCovariance folded = fold_in_held_columns( cov, J, A );

  Eigen::MatrixXd full( 5, 4 );
  full << J, A;
  const Eigen::MatrixXd expected = (full.transpose() * full).inverse();

  BOOST_REQUIRE_EQUAL( folded.cov11.rows(), 2 );
  BOOST_REQUIRE_EQUAL( folded.cov22.rows(), 2 );
  for( int r = 0; r < 2; ++r )
  {
    for( int c = 0; c < 2; ++c )
    {
      BOOST_CHECK_CLOSE_FRACTION( folded.cov11(r,c), expected(r,c), 1.0e-10 );
      BOOST_CHECK_CLOSE_FRACTION( folded.cov12(r,c), expected(r,2+c), 1.0e-10 );
      BOOST_CHECK_CLOSE_FRACTION( folded.cov22(r,c), expected(2+r,2+c), 1.0e-10 );
    }
  }

  // A held column the kept ones fully explain cannot be given a variance: it gets none, and leaves the
  //  kept covariance unchanged.
  Eigen::MatrixXd explained( 5, 1 );
  explained.col(0) = J.col(0) + 2.0*J.col(1);
  const FoldedCovariance degenerate = fold_in_held_columns( cov, J, explained );
  BOOST_CHECK_SMALL( degenerate.cov22(0,0), 1.0e-12 );
  BOOST_CHECK_SMALL( (degenerate.cov11 - cov).norm(), 1.0e-10 );
}//BOOST_AUTO_TEST_CASE( fold_in_held_columns_matches_full_inverse )


namespace
{
string g_test_file_dir;

void set_data_dir()
{
  static bool s_have_set = false;
  if( s_have_set )
    return;
  s_have_set = true;

  string datadir;
  for( int i = 1; i < framework::master_test_suite().argc; ++i )
  {
    const string arg = framework::master_test_suite().argv[i];
    if( SpecUtils::istarts_with( arg, "--datadir=" ) )
      datadir = arg.substr( 10 );
    if( SpecUtils::istarts_with( arg, "--testfiledir=" ) )
      g_test_file_dir = arg.substr( 14 );
  }

  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file( SpecUtils::append_path(datadir, "sandia.decay.xml") ),
                         "Pass --datadir=<InterSpec data dir>" );
  BOOST_REQUIRE_NO_THROW( InterSpec::setStaticDataDirectory( datadir ) );
  BOOST_REQUIRE( DecayDataBaseServer::database() );
}//void set_data_dir()
}//namespace


/** Fits `options` with and without `extra` (20 y, age not fit) on the same fixed ROIs, and checks the
 nested-model property: zero activity of `extra` reproduces the fit without it, so the fit with it must
 not end at a higher chi2 - which it did, by up to 6.5x while reporting success, when the solver stalled
 against that nuclide's zero bound.  If `extra` ends at zero, the user must be told it was not
 detected (with an upper limit), not that uncertainties are "under-stated". */
void check_nested_nuclide( RelActCalcAuto::Options options, const string &extra,
                           shared_ptr<const SpecUtils::Measurement> foreground,
                           shared_ptr<const DetectorPeakResponse> det,
                           const PeakFitUtils::CoarseResolutionType det_type )
{
  const SandiaDecay::Nuclide * const extra_nuc = DecayDataBaseServer::database()->nuclide( extra );
  BOOST_REQUIRE( extra_nuc );
  BOOST_REQUIRE( !options.rel_eff_curves.empty() && !options.rel_eff_curves[0].nuclides.empty() );

  vector<RelActCalcAuto::NucInputInfo> &nucs = options.rel_eff_curves[0].nuclides;
  const auto is_extra = [extra_nuc]( const RelActCalcAuto::NucInputInfo &nuc ) -> bool {
    return RelActCalcAuto::nuclide( nuc.source ) == extra_nuc;
  };
  if( std::none_of( begin(nucs), end(nucs), is_extra ) )
  {
    RelActCalcAuto::NucInputInfo nuc = nucs.front();
    nuc.source = extra_nuc;
    nuc.age = 20.0*PhysicalUnits::year;
    nuc.fit_age = false;
    nuc.fit_age_min.reset();
    nuc.fit_age_max.reset();
    nuc.starting_rel_act.reset();
    nuc.min_rel_act.reset();
    nuc.max_rel_act.reset();
    nucs.push_back( nuc );
  }

  {
    const RelActCalcAuto::RelActAutoSolution first
          = RelActCalcAuto::solve( options, foreground, nullptr, det, {}, det_type, nullptr );
    BOOST_REQUIRE( RelActCalcAuto::RelActAutoSolution::is_usable_status(first.m_status) );
    options.rois = first.m_final_roi_ranges;
    for( RelActCalcAuto::RoiRange &roi : options.rois )
      roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
  }

  RelActCalcAuto::Options without = options;
  {
    vector<RelActCalcAuto::NucInputInfo> &without_nucs = without.rel_eff_curves[0].nuclides;
    without_nucs.erase( std::remove_if( begin(without_nucs), end(without_nucs), is_extra ), end(without_nucs) );
  }

  const RelActCalcAuto::RelActAutoSolution without_sol
          = RelActCalcAuto::solve( without, foreground, nullptr, det, {}, det_type, nullptr );
  const RelActCalcAuto::RelActAutoSolution with_sol
          = RelActCalcAuto::solve( options, foreground, nullptr, det, {}, det_type, nullptr );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(without_sol.m_status),
                         "Fit without " << extra << " failed: " << without_sol.m_error_message );
  BOOST_REQUIRE_MESSAGE( RelActCalcAuto::RelActAutoSolution::is_usable_status(with_sol.m_status),
                         "Fit with " << extra << " failed: " << with_sol.m_error_message );

  cout << "check_nested_nuclide(" << extra << "): chi2 without " << without_sol.m_chi2
       << " (data " << without_sol.m_chi2_data << "), with " << with_sol.m_chi2
       << " (data " << with_sol.m_chi2_data << ")" << endl;
  for( const string &warning : with_sol.m_warnings )
    cout << "  warning: " << warning << endl;

  BOOST_CHECK_LE( with_sol.m_chi2, 1.005*without_sol.m_chi2 );
  BOOST_CHECK_LE( with_sol.m_chi2_data, 1.005*without_sol.m_chi2_data );
  BOOST_CHECK( !with_sol.m_optimizer_returned_seed );
  BOOST_CHECK( !without_sol.m_optimizer_returned_seed );

  for( const string &warning : with_sol.m_warnings )
    BOOST_CHECK_MESSAGE( warning.find("under-stated") == string::npos, "Old bound message: " << warning );

  bool at_zero = false;
  for( const RelActCalcAuto::NuclideRelAct &act : with_sol.m_rel_activities.at(0) )
    at_zero = at_zero || ((act.name() == extra) && !(act.rel_activity > 0.0));
  // Whether an absent nuclide lands on zero is a coin flip on the counting noise; when it does, the user
  //  is told it was not detected (a tiny positive activity at the bound may be described the same way).
  if( at_zero )
  {
    bool told = false;
    for( const string &warning : with_sol.m_warnings )
      told = told || ((warning.find(extra + " activity fit to zero (not detected") != string::npos)
                      && (warning.find("upper limit") != string::npos));
    BOOST_CHECK_MESSAGE( told, extra << " at zero without a 'not detected' message" );
  }
}//void check_nested_nuclide(...)


BOOST_AUTO_TEST_CASE( added_nuclide_is_a_nested_model )
{
  set_data_dir();

  const string spec_path = SpecUtils::append_path( SpecUtils::append_path(g_test_file_dir, "RelActAutoReport"),
                                                   "U235_Unshielded_6000.n42" );
  BOOST_REQUIRE_MESSAGE( SpecUtils::is_file(spec_path), "Pass --testfiledir=<target/testing/test_data>" );

  auto meas = make_shared<SpecMeas>();
  BOOST_REQUIRE( meas->load_file( spec_path, SpecUtils::ParserType::Auto ) );

  shared_ptr<const SpecUtils::Measurement> foreground;
  for( const shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
  {
    if( m && (m->source_type() != SpecUtils::SourceType::Background) && !foreground )
      foreground = m;
  }
  BOOST_REQUIRE( foreground );
  const shared_ptr<const DetectorPeakResponse> det = meas->detector();
  BOOST_REQUIRE( det );
  unique_ptr<RelActCalcAuto::RelActAutoGuiState> state = meas->getRelActAutoGuiState();
  BOOST_REQUIRE( state );
  const PeakFitUtils::CoarseResolutionType det_type
                        = PeakFitUtils::effective_det_type( meas->peakFitDetPrefs(), foreground, meas );

  // U232 is in this item (its chain lines are clearly seen); Ba133 is not.
  check_nested_nuclide( state->options, "U232", foreground, det, det_type );
  check_nested_nuclide( state->options, "Ba133", foreground, det, det_type );
}//BOOST_AUTO_TEST_CASE( added_nuclide_is_a_nested_model )
