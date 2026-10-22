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
}//BOOST_AUTO_TEST_CASE( blocked_parameter_classification )


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
