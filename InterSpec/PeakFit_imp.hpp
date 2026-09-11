#ifndef PeakFit_imp_h
#define PeakFit_imp_h

#include <cmath>
#include <atomic>
#include <cstdio>
#include <vector>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <type_traits>

#include "Eigen/Dense"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"

namespace PeakFit
{

/** Half-width, in units of the channel uncertainty, of the smooth one-sided floor applied to the
 continuum-fit target `y = max(data - peaks, 0)/uncert` in #fit_continuum.

 The old hard `(data > peaks) ? (data - peaks)/uncert : 0` branch fired in EVERY empty channel
 carrying any model peak mass (`peaks >= 0` always, and at `data == 0` the uncertainty floors to
 1), so its zero-gradient side was the COMMON case in sparse spectra - deforming the autodiff
 Jacobian (AD-vs-finite-difference disagreements at converged fits traced here).  The quadratic
 hinge keeps the same intent (dont let overshooting peaks drag the continuum negative) while making
 the target C1: values differ from the old clamp by at most `hinge/4` (= this constant / 4, in
 sigma units) inside the +-hinge blend zone, and are exactly equal outside it.
 */
static const double ns_cont_target_hinge_nsigma = 0.25;

/** C1 smooth `max(x, 0)` with compact support (exact outside |x| >= r; quadratic blend inside;
 `>= max(x,0)` everywhere with excess <= r/4).  A local copy of RelActCalc::qmax_hinge (kept
 header-local to avoid layering PeakFit on RelActCalc) - keep the two in sync.
 */
template<typename T>
T smooth_max_zero( const T &x, const double r )
{
  assert( r > 0.0 );

  if( x >= r )
    return x;
  if( x <= -r )
    return T(0.0);
  return (x + r)*(x + r) / (4.0*r);
}//smooth_max_zero(...)


// Defined in PeakDists_imp.hpp; pulled in here so the (many) unqualified uses below still resolve.
using PeakDists::unit_pdf_to_cdf;


/** This function fits the polynomial continuum for a region with a number of fixed amplitude peaks.
 *
 * When `ScalarType` is a `ceres::Jet<>` and the design matrix does not depend on the fit
 * parameters (all non-peak-CDF continua, with a derivative-free `ref_energy` - e.g. RelActCalcAuto
 * with its fixed ROI channel ranges), the matrix is factored in plain `double` and the
 * factorization applied to the Jet-valued right-hand side.  The LLS solution is linear in the RHS,
 * so this is mathematically identical to solving in Jet arithmetic, but ~an order of magnitude
 * cheaper (a Jet<double,32> multiply costs ~33 double multiplies), and it sidesteps the NaN
 * derivative values the Jet-valued `svd.solve(y)` could produce for peaks with a skew.
 * The peak-CDF step continua build the design matrix from the peaks themselves, so they always
 * use the Jet-valued solve.
 */
template<typename PeakType, typename ScalarType>
void fit_continuum( const float * const x,
                    const float * const data,
                    const float * const data_uncert,
                    const size_t nbin,
                    const PeakContinuum::OffsetType cont_type,
                    const ScalarType ref_energy,
                    const std::vector<PeakType> &fixedAmpPeaks,
                    const bool multithread,
                    ScalarType *continuum_coeffs,
                    ScalarType *peak_counts )
{
  using namespace std;

  static const double MIN_CHANNEL_UNCERT = 1.0;

  const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_linear_fit_pars( cont_type ) );
  const bool step_continuum = PeakContinuum::is_step_continuum( cont_type );
  const bool cdf_step = PeakContinuum::is_peak_cdf_step_continuum( cont_type );
  // Number of peak-CDF step coefficients: 1 for FlatStepCDF/LinearStepCDF, 2 for BiLinearStepCDF.
  //  Here every peak has a fixed amplitude, so these are linear and get their own LLS columns.
  const Eigen::Index num_step_terms
             = static_cast<Eigen::Index>( PeakContinuum::num_cdf_step_pars( cont_type ) );

  // Loosely following:
  //   https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  // using the SVD (slowest, but most accurate) method

  const Eigen::Index num_poly_terms = static_cast<Eigen::Index>( num_polynomial_terms  );

  // The step coefficients are solved by LLS here, since the peak amplitudes are known.
  const Eigen::Index num_lls_terms = num_poly_terms + num_step_terms;

  Eigen::VectorX<ScalarType> y( static_cast<Eigen::Index>(nbin) );
  std::vector<double> uncerts( nbin, 1.0 );

  double roi_data_sum = 0.0;

#if( PEAK_CONTINUUM_DATA_STEP_SUBTRACT )
  double min_data_val = static_cast<double>( data[0] );
  for( size_t row = 0; row < nbin; ++row )
  {
    roi_data_sum += (std::max)( data[row], 0.0f );
    min_data_val = (std::min)( min_data_val, static_cast<double>( data[row] ) );
  }
  if( step_continuum && !cdf_step )
  {
    min_data_val = (std::max)( min_data_val, 0.0 );
    roi_data_sum -= min_data_val * nbin;
  }else
  {
    // The CDF step types build their step from the peaks, not from the cumulative data.
    min_data_val = 0.0;
  }
#else
  const double min_data_val = 0.0;
  for( size_t row = 0; row < nbin; ++row )
    roi_data_sum += (std::max)( data[row], 0.0f );
#endif

  // Zero out the destination count array
  for( size_t row = 0; row < nbin; ++row )
    peak_counts[row] = ScalarType(0.0);

  // Add the Gaussian + Skew component of the peaks to destination counts
  const size_t nfixedpeak = fixedAmpPeaks.size();

  if( multithread && (nfixedpeak > 8) ) //8 is arbitrary.
  {
    // TODO: multi-thread computation needs to be evaluated more hollistically both here and in #RelActAutoSolution::eval
    const unsigned nthread = (std::min)( 16, (std::min)( static_cast<int>(nfixedpeak), (std::max)( 1, static_cast<int>( std::thread::hardware_concurrency() ) ) ) );

    vector<vector<ScalarType>> results( nthread );

    SpecUtilsAsync::ThreadPool pool;

    for( size_t thread_index = 0; thread_index < nthread; ++thread_index )
    {
      pool.post( [nbin, thread_index, nfixedpeak, nthread, &results, &fixedAmpPeaks, &x](){
        results[thread_index].resize( nbin, ScalarType(0.0) );

        for( size_t peak_index = thread_index; peak_index < nfixedpeak; peak_index += nthread )
        {
          fixedAmpPeaks[peak_index].gauss_integral( x, results[thread_index].data(), nbin );
        }//for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
      } );
    }//for( size_t thread_index = 0; thread_index < nthread; ++thread_index )

    pool.join();

    /*
    vector<std::thread> threads( nthread );
    for( size_t thread_index = 0; thread_index < nthread; ++thread_index )
    {
      threads[thread_index] = std::thread( [nbin, thread_index, nfixedpeak, nthread, &results, &fixedAmpPeaks, &x](){
        results[thread_index].resize( nbin, ScalarType(0.0) );

        for( size_t peak_index = thread_index; peak_index < nfixedpeak; peak_index += nthread )
        {
          fixedAmpPeaks[peak_index].gauss_integral( x, results[thread_index].data(), nbin );
        }//for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
      } );
    }

    for( size_t thread_index = 0; thread_index < nthread; ++thread_index )
    {
      threads[thread_index].join();
    }
    */

    // TODO: use Eigen to vectorize these sums
    for( size_t thread_index = 0; thread_index < nthread; ++thread_index )
    {
      for( size_t i = 0; i < nbin; ++i )
        peak_counts[i] += results[thread_index][i];
    }
  }else
  {
    for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
    {
      fixedAmpPeaks[peak_index].gauss_integral( x, peak_counts, nbin );
    }//for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
  }//if( multithread && (nfixedpeak > 8) ) / else

  // For CDF step types, precompute the amplitude-weighted peak CDF per channel:
  //   cdf_amp_sum[row] = SUM_j( amp_j * CDFbar_j(chan_center) ) * dx
  // which is the basis the step coefficients multiply.  Step coefficient k additionally carries a
  // factor of (chan_center - ref_energy)^k, applied when the design matrix is built.
  //
  // Unit-area peak integrals come from photopeak_function_integral (templated on ScalarType, so
  // Ceres Jet derivative info is preserved), turned into a CDF by cumulative summation in
  // unit_pdf_to_cdf.  That cumulative sum starts at the ROI's first channel edge and saturates at
  // its last, i.e. it is exactly the ROI-anchored CDFbar the evaluator uses - see
  // PeakContinuum::cdf_step_anchor_energies(...).
  vector<ScalarType> cdf_amp_sum;
  if( cdf_step )
  {
    cdf_amp_sum.resize( nbin, ScalarType(0.0) );

    vector<ScalarType> pdf_per_channel( nbin );
    vector<ScalarType> cdf_per_channel( nbin );

    for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
    {
      const ScalarType &amp_val = fixedAmpPeaks[peak_index].amplitude();
      const ScalarType &mean_val = fixedAmpPeaks[peak_index].mean();
      const ScalarType &sigma_val = fixedAmpPeaks[peak_index].sigma();
      const PeakDef::SkewType skew = fixedAmpPeaks[peak_index].skewType();

      // Get skew parameters as ScalarType pointer
      const ScalarType *skew_pars = nullptr;
      if constexpr ( std::is_same_v<PeakType, PeakDef> )
      {
        skew_pars = (skew != PeakDef::SkewType::NoSkew)
                    ? (fixedAmpPeaks[peak_index].coefficients() + PeakDef::CoefficientType::SkewPar0)
                    : nullptr;
      }else
      {
        skew_pars = (skew != PeakDef::SkewType::NoSkew)
                    ? fixedAmpPeaks[peak_index].skew_parameters()
                    : nullptr;
      }

      // Compute unit-area peak PDF integral per channel, then CDF at channel centers
      std::fill( pdf_per_channel.begin(), pdf_per_channel.end(), ScalarType(0.0) );
      PeakDists::photopeak_function_integral( mean_val, sigma_val, ScalarType(1.0),
                                              skew, skew_pars, nbin, x, pdf_per_channel.data() );
      unit_pdf_to_cdf( pdf_per_channel.data(), cdf_per_channel.data(), nbin );

      for( size_t row = 0; row < nbin; ++row )
      {
        const double dx = x[row+1] - x[row];
        cdf_amp_sum[row] += amp_val * cdf_per_channel[row] * dx;
      }
    }//for( peak_index )
  }//if( cdf_step )

  // Fills the LLS right-hand-side `y` (and `uncerts`), plus the design matrix `A_mat`.
  //  `A_mat`s scalar type is plain `double` when the design matrix does not depend on the fit
  //  parameters (see the function doc-comment), and `ScalarType` otherwise: the design-matrix
  //  entries depend on the fit parameters only through `ref_e` and - for the peak-CDF step
  //  continua - the peak amplitudes/shapes; the channel edges, data, and uncertainties are
  //  all constants.
  const auto fill_design_and_rhs = [&]( auto &A_mat, const auto ref_e )
  {
    // RT is `double` when the design matrix is parameter-independent, and `ScalarType` otherwise
    using RT = std::decay_t<decltype(ref_e)>;

    double step_cumulative_data = 0.0;

    for( size_t row = 0; row < nbin; ++row )
    {
      const double data_counts = data[row];
      const double data_counts_uncert = data_uncert ? data_uncert[row] : sqrt(data[row]);
      const double x0 = x[row];
      const double x1 = x[row+1];
      const RT x0_rel = x0 - ref_e;
      const RT x1_rel = x1 - ref_e;

      const double uncert = (data_counts_uncert > MIN_CHANNEL_UNCERT) ? data_counts_uncert : 1.0;

      uncerts[row] = uncert;

      if( step_continuum && !cdf_step )
        step_cumulative_data += (data_counts - min_data_val);

      // Smooth one-sided floor on the fit target - C1 in the peak parameters, same intent as the
      //  old hard `(data > peaks) ? (data - peaks)/uncert : 0` clamp (see ns_cont_target_hinge_nsigma).
      const ScalarType target_deficit = data_counts - peak_counts[row];
      y(row) = smooth_max_zero( target_deficit, ns_cont_target_hinge_nsigma*uncert ) / uncert;

#if( PERFORM_DEVELOPER_CHECKS )
      {// Opt-in activation statistics for the (former) clamp - how often the floored side is hit.
        static const bool s_clamp_stats = ( std::getenv("RELACT_CONT_CLAMP_STATS") != nullptr );
        if( s_clamp_stats )
        {
          static std::atomic<uint64_t> s_num_channels{ 0 }, s_num_floored{ 0 };
          double deficit_val;
          if constexpr ( std::is_same_v<ScalarType,double> )
            deficit_val = target_deficit;
          else
            deficit_val = target_deficit.a;
          if( deficit_val <= 0.0 )
            s_num_floored += 1;
          const uint64_t nchan = ++s_num_channels;
          if( (nchan & 0x3FFFF) == 0 )  //print every ~262k channels
            fprintf( stderr, "RELACT_CONT_CLAMP_STATS: %llu of %llu channels (%.2f%%) on the floored side\n",
                     static_cast<unsigned long long>(s_num_floored.load()),
                     static_cast<unsigned long long>(nchan),
                     (100.0*s_num_floored.load())/nchan );
        }//if( s_clamp_stats )
      }
#endif //PERFORM_DEVELOPER_CHECKS

      for( Eigen::Index col = 0; col < num_poly_terms; ++col )
      {
        const double exp = col + 1.0;

        if( !cdf_step && step_continuum
            && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
            && (col == (num_polynomial_terms - 1)) )
        {
          // This logic mirrors that of PeakContinuum::offset_integral(...), and code
          // If you change it in one place - change it in here, below, and in offset_integral.
          const double frac_data = (roi_data_sum > 0.0)
              ? (step_cumulative_data - 0.5*(data_counts - min_data_val)) / roi_data_sum : 0.5;
          const double contribution = frac_data * (x1 - x0);

          A_mat(row,col) = RT(contribution / uncert);

          check_jet_for_NaN( A_mat(row,col) );
        }else if( !cdf_step && step_continuum && (num_polynomial_terms == 4) )
        {
          const double frac_data = (roi_data_sum > 0.0)
              ? (step_cumulative_data - 0.5*(data_counts - min_data_val)) / roi_data_sum : 0.5;

          RT contrib( 0.0 );
          switch( col )
          {
            case 0: contrib = (1.0 - frac_data) * (x1_rel - x0_rel);                     break;
            case 1: contrib = 0.5 * (1.0 - frac_data) * (x1_rel*x1_rel - x0_rel*x0_rel); break;
            case 2: contrib = frac_data * (x1_rel - x0_rel);                             break;
            case 3: contrib = 0.5 * frac_data * (x1_rel*x1_rel - x0_rel*x0_rel);         break;
            default: assert( 0 ); break;
          }//switch( col )

          A_mat(row,col) = contrib / uncert;

          check_jet_for_NaN( contrib );
          check_jet_for_NaN( uncert );
          check_jet_for_NaN( A_mat(row,col) );
        }else
        {
          const RT contribution = (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));

          A_mat(row,col) = contribution / uncert;

          check_jet_for_NaN( contribution );
          check_jet_for_NaN( uncert );
          check_jet_for_NaN( A_mat(row,col) );
        }
      }//for( int order = 0; order < maxorder; ++order )

      // Peak-CDF step columns; coefficient k multiplies cdf_amp_sum * (chan_center - ref)^k.
      //  These are built from the peaks themselves, so a CDF step continuum never takes the
      //  constant-design (RT == double) path.
      for( Eigen::Index k = 0; k < num_step_terms; ++k )
      {
        if constexpr ( std::is_same_v<RT,ScalarType> )
        {
          const ScalarType center_rel = 0.5*(x0_rel + x1_rel);
          ScalarType contrib = cdf_amp_sum[row];
          for( Eigen::Index e = 0; e < k; ++e )
            contrib *= center_rel;

          A_mat(row, num_poly_terms + k) = contrib / uncert;

          check_jet_for_NaN( contrib );
          check_jet_for_NaN( A_mat(row, num_poly_terms + k) );
        }else
        {
          assert( 0 ); //A peak-CDF continuum implies a parameter-dependent design matrix
        }
      }//for( Eigen::Index k = 0; k < num_step_terms; ++k )
    }//for( size_t row = 0; row < nbin; ++row )
  };//fill_design_and_rhs lambda


  Eigen::VectorX<ScalarType> coeffs;    //The solved continuum coefficients (num_lls_terms of them)
  Eigen::VectorX<ScalarType> cont_vals; //The continuum counts, divided by `uncerts` (i.e., A*coeffs), per channel

  // The reference solve: build the design matrix in `ScalarType` and solve the LLS problem in
  //  that arithmetic (i.e., through `ceres::Jet` when auto-differentiating).
  const auto solve_generic = [&]()
  {
    Eigen::MatrixX<ScalarType> A( static_cast<Eigen::Index>(nbin), num_lls_terms );

    fill_design_and_rhs( A, ref_energy );

#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
    const Eigen::JacobiSVD<Eigen::MatrixX<ScalarType>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
    const Eigen::BDCSVD<Eigen::MatrixX<ScalarType>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif

    coeffs = svd.solve(y); // coeffs will contain [c_0, c_1, c_2, c_3]

    check_jet_array_for_NaN( coeffs.data(), coeffs.size() );

    //If time is a real issue, we could try other methods, like:
    //Eigen::VectorX<ScalarType> coeffs = A.colPivHouseholderQr().solve(y);

    //If we wanted to get the covariance matrix and/or parameter uncertainties, we could (unchecked):
    //Calculate residuals
    //Eigen::VectorX<ScalarType> residuals = y - A * coeffs;
    //Calculate variance of residuals
    //double sigma_squared = (residuals.squaredNorm()) / (residuals.size() - A.cols());
    //Calculate covariance matrix
    //Eigen::MatrixXd covariance_matrix = sigma_squared * (A_weighted.transpose() * A_weighted).inverse();
    //Extract uncertainties (standard errors)
    //Eigen::VectorXd uncertainties = covariance_matrix.diagonal().array().sqrt();
    //Calculate correlation matrix
    //Eigen::MatrixXd correlation_matrix = covariance_matrix;
    //for (int i = 0; i < correlation_matrix.rows(); ++i) {
    //  for (int j = 0; j < correlation_matrix.cols(); ++j) {
    //    correlation_matrix(i, j) /= (uncertainties(i) * uncertainties(j));
    //  }
    //}

    cont_vals = A * coeffs;
  };//solve_generic lambda

  try
  {
    bool solved_with_const_design = false;

#if( defined(CERES_PUBLIC_JET_H_) )
    if constexpr ( !std::is_same_v<ScalarType,double> )
    {
      static_assert( std::is_same_v<ScalarType,ceres::Jet<double,ScalarType::DIMENSION>>,
                    "fit_continuum: a non-double ScalarType must be ceres::Jet<double,N>" );

      // The design matrix depends on the fit parameters only through `ref_energy` (for the
      //  non-peak-CDF continua), so when `ref_energy` carries no derivative information we can
      //  factor the matrix in plain doubles, and solve for the value and every derivative
      //  direction of the Jet RHS through that one factorization - see the function doc-comment.
      bool design_is_const = !cdf_step;
      for( int i = 0; design_is_const && (i < ScalarType::DIMENSION); ++i )
        design_is_const = (ref_energy.v[i] == 0.0);

      if( design_is_const )
      {
        constexpr int num_deriv = ScalarType::DIMENSION;

        Eigen::MatrixXd A_d( static_cast<Eigen::Index>(nbin), num_lls_terms );

        fill_design_and_rhs( A_d, ref_energy.a );

#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
        const Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::ComputeThinU | Eigen::ComputeThinV> svd( A_d );
#else
        const Eigen::BDCSVD<Eigen::MatrixXd> svd( A_d, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif

        // Solve the value column and each derivative column with the same (thresholded)
        //  pseudo-inverse, so the derivative parts are the exact derivatives of the value part.
        Eigen::MatrixXd rhs( static_cast<Eigen::Index>(nbin), 1 + num_deriv );
        for( size_t row = 0; row < nbin; ++row )
        {
          rhs(row,0) = y(row).a;
          for( int j = 0; j < num_deriv; ++j )
            rhs(row,1+j) = y(row).v[j];
        }

        const Eigen::MatrixXd solved = svd.solve( rhs );
        const Eigen::MatrixXd predicted = A_d * solved;   // nbin x (1 + num_deriv)

        coeffs.resize( num_lls_terms );
        for( Eigen::Index i = 0; i < num_lls_terms; ++i )
        {
          coeffs(i).a = solved(i,0);
          for( int j = 0; j < num_deriv; ++j )
            coeffs(i).v[j] = solved(i,1+j);
        }

        cont_vals.resize( static_cast<Eigen::Index>(nbin) );
        for( size_t row = 0; row < nbin; ++row )
        {
          cont_vals(row).a = predicted(row,0);
          for( int j = 0; j < num_deriv; ++j )
            cont_vals(row).v[j] = predicted(row,1+j);
        }

        check_jet_array_for_NaN( coeffs.data(), coeffs.size() );

        solved_with_const_design = true;

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        // Cross-check the double-factorization against the Jet-valued SVD solve; opt-in via
        //  environment variable, since it roughly doubles the cost of every continuum fit.
        static const bool s_check_const_design_lls = ( std::getenv("RELACT_CHECK_CONT_LLS") != nullptr );
        if( s_check_const_design_lls )
        {
          const Eigen::VectorX<ScalarType> fast_coeffs = coeffs;
          const Eigen::VectorX<ScalarType> fast_cont_vals = cont_vals;

          solve_generic(); //overwrites `coeffs` and `cont_vals` with the Jet-SVD reference values

          for( Eigen::Index i = 0; i < num_lls_terms; ++i )
          {
            // The Jet-valued SVD can produce NaN derivatives (see doc-comment); the double
            //  factorization is exact there, so only compare where the reference is finite.
            if( std::isfinite(coeffs(i).a) && std::isfinite(fast_coeffs(i).a) )
            {
              const double tol = 1.0E-6 * (std::max)( 1.0, (std::max)( fabs(coeffs(i).a), fabs(fast_coeffs(i).a) ) );
              assert( fabs(coeffs(i).a - fast_coeffs(i).a) < tol );
            }

            for( int j = 0; j < num_deriv; ++j )
            {
              if( std::isfinite(coeffs(i).v[j]) && std::isfinite(fast_coeffs(i).v[j]) )
              {
                const double tol = 1.0E-4 * (std::max)( 1.0, (std::max)( fabs(coeffs(i).v[j]), fabs(fast_coeffs(i).v[j]) ) );
                assert( fabs(coeffs(i).v[j] - fast_coeffs(i).v[j]) < tol );
              }
            }
          }//for( loop over coefficients )

          coeffs = fast_coeffs;       //keep the fast-path values, so behavior doesnt
          cont_vals = fast_cont_vals; //  depend on the environment variable
        }//if( s_check_const_design_lls )
#endif //PERFORM_DEVELOPER_CHECKS && !NDEBUG
      }//if( design_is_const )
    }//if constexpr( ScalarType is a ceres::Jet )
#endif //defined(CERES_PUBLIC_JET_H_)

    if( !solved_with_const_design )
      solve_generic();

    assert( coeffs.size() == num_lls_terms );
    assert( coeffs.rows() == num_lls_terms );

    // Copy polynomial coefficients
    for( Eigen::Index i = 0; i < num_poly_terms; ++i )
      continuum_coeffs[i] = coeffs(i);

    // The step coefficients follow the polynomial terms, matching PeakContinuum's parameter order.
    //  e.g. FlatStepCDF: [constant, step0]; LinearStepCDF: [constant, linear, step0];
    //       BiLinearStepCDF: [constant, linear, step0, step1]
    for( Eigen::Index k = 0; k < num_step_terms; ++k )
      continuum_coeffs[num_poly_terms + k] = coeffs(num_poly_terms + k);

    for( size_t bin = 0; bin < nbin; ++bin )
    {
      // NOTE: smoothing this hard clamp with the same hinge as the fit-target floor was tried
      //  (2026-07) and REVERTED: the idb_enrichment_check gate showed constrained shallow-surface
      //  fits reshuffling into worse basins (one file chi2 +135), while the unconstrained corpus
      //  slightly improved - net not neutral-or-better.  The zero-gradient side here is also far
      //  rarer than the fit-target floors was (the fitted continuum seldom dips negative).
      ScalarType y_continuum = cont_vals(bin) * uncerts[bin];
      double scalar_continuum = 0.0;
      if constexpr ( std::is_same_v<ScalarType,double> )
        scalar_continuum = y_continuum;
      else
        scalar_continuum = y_continuum.a;

      if( scalar_continuum < 0.0 )
        y_continuum = ScalarType(0.0);
      peak_counts[bin] += y_continuum;
    }
  }catch( std::exception &e )
  {
    cerr << "RelActCalcAuto::fit_continuum(...): caught: " << e.what() << endl;

    throw runtime_error( "RelActCalcAuto::fit_continuum(...): trouble finding coeffs." );
  }//try / catch
}//void fit_continuum(...)


/** Fits the continuum and amplitude of peaks with specified means and sigmas, over the data range specified.
 Uses a matrix based linear regression fitter to perform the least linear squares

 @param energies The lower-channel energies of ROI.  ROI defined by energies[0] to energies[nbin]. Must be of at least length nbin+1.
 @param data The channel counts of the ROI.  Must be of at least length nbin.
 @param variances Optional variances (uncertainties squared) of the data.  If specified - you must make sure each value is positive and non-zero.
        If not specified, will use the square root of the data value, or 1.0, whichever is larger.
 @param nbin The number of channels in the ROI.
 @param num_polynomial_terms The number of polynomial continuum terms to fit for.
        0 is no continuum (untested), 1 is constant, 2 is linear sloped continuum, etc
 @param cont_type The continuum type; determines polynomial vs step vs CDF step behavior.
 @param step_coeffs For the CDF step types, the `PeakContinuum::num_cdf_step_pars(cont_type)` step
        coefficients, in the order they appear in the continuum's parameters (constant term first,
        then the energy slope for BiLinearStepCDF).  These are bilinear with the peak amplitudes, so
        they cannot be solved here and are taken as known inputs from the non-linear solver.
        May be null, which means "all step coefficients are zero"; ignored for non-CDF types.
 @param means The peak means, in keV
 @param sigmas The peak sigmas, in keV
 @param fixedAmpPeaks The fixed amplitude peaks in the ROI, that we are not fitting for
 @param skew_type The skew type to use for peaks that are being fit
 @param skew_parameters The parameters that specify the skew; must have `PeakDef::num_skew_parameters(SkewType)` number
        of entries; if no skew is used, may be `nullptr`.
 @param[out] amplitudes The fit peak amplitudes
 @param[out] continuum_coeffs The fit continuum coefficients
 @param[out] amplitudes_uncerts The (statistical) uncertainties for the amplitudes
 @param[out] continuum_coeffs_uncerts The (statistical) uncertainties for the continuum coefficients
 @param[out] peak_counts Optional array (at least of length `nbin`) to place the summed counts of the peaks plus continuum;
             the sums of each channel are added to each array element (i.e., the elements arent set equal to sum, so you should
             zero thier values, if you want that)

 @returns The chi2 of the ROI

 Skew uncertainties are also not taken into account in determining amplitude or continuum uncertainties.

 Throws exception upon ill-posed input.
 */
template<typename PeakType, typename ScalarType>
ScalarType fit_amp_and_offset_imp( const float *x,
                                  const float *data,
                                  const float *variances,
                                  const size_t nbin,
                          const PeakContinuum::OffsetType cont_type,
                          const std::type_identity_t<ScalarType> * const step_coeffs,
                          const ScalarType ref_energy,
                          const std::vector<ScalarType> &means,
                          const std::vector<ScalarType> &sigmas,
                          const std::vector<PeakType> &fixedAmpPeaks,
                          const PeakDef::SkewType skew_type,
                          const ScalarType *skew_parameters,
                          std::vector<ScalarType> &amplitudes,
                          std::vector<ScalarType> &continuum_coeffs,
                          std::vector<ScalarType> &amplitudes_uncerts,
                          std::vector<ScalarType> &continuum_coeffs_uncerts,
                          ScalarType * const peak_counts )
{
  using namespace std;

  const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_linear_fit_pars( cont_type ) );
  const bool step_continuum = PeakContinuum::is_step_continuum( cont_type );
  const bool cdf_step = PeakContinuum::is_peak_cdf_step_continuum( cont_type );
  // Number of peak-CDF step coefficients; these are bilinear with the peak amplitudes, so they are
  //  passed in as known values by the non-linear solver rather than being solved here.
  const size_t num_step_terms = PeakContinuum::num_cdf_step_pars( cont_type );

  if( sigmas.size() != means.size() )
    throw runtime_error( "fit_amp_and_offset_imp: invalid input" );

  assert( (skew_type == PeakDef::SkewType::NoSkew) || skew_parameters );
  if( !skew_parameters && (skew_type != PeakDef::SkewType::NoSkew) )
    throw std::logic_error( "fit_amp_and_offset_imp: skew pars not provided" );

  // A null `step_coeffs` means all step coefficients are zero; several callers deliberately fit a
  //  CDF step continuum's polynomial with no step (e.g. a null-hypothesis continuum).
  const bool have_step_coeffs = (num_step_terms > 0) && step_coeffs;

  const size_t npeaks = sigmas.size();

  const Eigen::Index num_poly_terms = static_cast<Eigen::Index>( num_polynomial_terms  );
  const Eigen::Index nfit_terms = static_cast<int>( num_poly_terms + npeaks );

  Eigen::MatrixX<ScalarType> A( static_cast<Eigen::Index>(nbin), nfit_terms );
  Eigen::VectorX<ScalarType> y( static_cast<Eigen::Index>(nbin) );
  Eigen::VectorX<ScalarType> uncerts( static_cast<Eigen::Index>(nbin) );


  ScalarType roi_data_sum( 0.0 ), step_cumulative_data( 0.0 );

  const double roi_lower = x[0];
  const double roi_upper = x[nbin];

#if( PEAK_CONTINUUM_DATA_STEP_SUBTRACT )
  double min_data_val = static_cast<double>( data[0] );
  for( size_t row = 0; row < nbin; ++row )
  {
    roi_data_sum += (std::max)( static_cast<double>( data[row] ), 0.0 );
    min_data_val = (std::min)( min_data_val, static_cast<double>( data[row] ) );
  }
  if( step_continuum && !cdf_step )
  {
    min_data_val = (std::max)( min_data_val, 0.0 );
    roi_data_sum -= ScalarType( min_data_val * static_cast<double>( nbin ) );
  }else
  {
    min_data_val = 0.0;
  }
#else
  const double min_data_val = 0.0;
  for( size_t row = 0; row < nbin; ++row )
    roi_data_sum += (std::max)( static_cast<double>(data[row]), 0.0 );
#endif

  const ScalarType avrg_data_val = roi_data_sum / static_cast<double>(nbin);

  // TODO: implement multi-threaded computation
  // 20250127: Using SpecUtilsAsync::ThreadPool, it looks like calling `pool.join()` is causing significant and unreasonable delays
  //           on macOS (using GCD, at least).  This is likely a problem with
  //           `SpecUtilsAsync::ThreadPool` - I would guess when creating ThreadPools inside of
  //           other ThreadPools, but for the moment will just do this single threaded, which is
  //           like 20 times faster for an example problem

  // Step coefficient at a channel's centre: s0 + s1*E' + ... , E' relative to `ref_energy`.
  //  Averaging the relative edges (rather than the absolute ones) matters: `x` is float and the
  //  absolute energies are ~1e3, so summing there loses several digits of E'.
  const auto step_coeff_for_channel = [&]( const size_t channel ) -> ScalarType {
    if( !have_step_coeffs )
      return ScalarType( 0.0 );

    const ScalarType center_rel = ScalarType(0.5) * ((ScalarType(x[channel]) - ref_energy)
                                                      + (ScalarType(x[channel+1]) - ref_energy));
    ScalarType step_coeff = step_coeffs[0];
    ScalarType energy_pow = center_rel;
    for( size_t k = 1; k < num_step_terms; ++k )
    {
      step_coeff += step_coeffs[k] * energy_pow;
      energy_pow *= center_rel;
    }

    return step_coeff;
  };//step_coeff_for_channel

  const size_t nfixedpeak = fixedAmpPeaks.size();
  vector<ScalarType> fixed_peak_contrib( nfixedpeak ? nbin : size_t(0), ScalarType(0.0) );

  // For the CDF step types the fixed-amplitude peaks contribute `amp_j * CDFbar_j` to the step just
  //  as the peaks being fit do - `PeakContinuum::offset_integral(...)` and `fit_continuum(...)` both
  //  sum over every peak in the ROI, so leaving them out here would make the fitted continuum
  //  disagree with the drawn one on any ROI holding a peak amplitude fixed.  Their amplitudes are
  //  known, so this is a constant per channel rather than a design-matrix column.
  vector<ScalarType> fixed_step_contrib( (nfixedpeak && cdf_step) ? nbin : size_t(0), ScalarType(0.0) );

  if( nfixedpeak )
  {
    ScalarType * const fixed_contrib = &(fixed_peak_contrib[0]);
    for( size_t peak_index = 0; peak_index < fixedAmpPeaks.size(); ++peak_index )
      fixedAmpPeaks[peak_index].gauss_integral( x, fixed_contrib, nbin );
  }//if( nfixedpeak )

  if( !fixed_step_contrib.empty() )
  {
    vector<ScalarType> pdf_per_channel( nbin ), cdf_per_channel( nbin );

    for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
    {
      const ScalarType &amp_val = fixedAmpPeaks[peak_index].amplitude();
      const PeakDef::SkewType skew = fixedAmpPeaks[peak_index].skewType();

      const ScalarType *fixed_skew_pars = nullptr;
      if( skew != PeakDef::SkewType::NoSkew )
      {
        if constexpr ( std::is_same_v<PeakType, PeakDef> )
          fixed_skew_pars = fixedAmpPeaks[peak_index].coefficients() + PeakDef::CoefficientType::SkewPar0;
        else
          fixed_skew_pars = fixedAmpPeaks[peak_index].skew_parameters();
      }

      std::fill( begin(pdf_per_channel), end(pdf_per_channel), ScalarType(0.0) );
      PeakDists::photopeak_function_integral( fixedAmpPeaks[peak_index].mean(),
                                              fixedAmpPeaks[peak_index].sigma(), ScalarType(1.0),
                                              skew, fixed_skew_pars, nbin, x, pdf_per_channel.data() );
      unit_pdf_to_cdf( pdf_per_channel.data(), cdf_per_channel.data(), nbin );

      for( size_t channel = 0; channel < nbin; ++channel )
      {
        const ScalarType dx = ScalarType( x[channel+1] - x[channel] );
        fixed_step_contrib[channel] += step_coeff_for_channel(channel)
                                       * amp_val * cdf_per_channel[channel] * dx;
      }
    }//for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
  }//if( !fixed_step_contrib.empty() )

  for( size_t row = 0; row < nbin; ++row )
  {
    ScalarType dataval = ScalarType( static_cast<double>(data[row]) );

    const double x0 = x[row];
    const double x1 = x[row+1];

    const ScalarType x0_rel = x0 - ref_energy;
    const ScalarType x1_rel = x1 - ref_energy;

    //const double uncert = (dataval > 0.0 ? sqrt(dataval) : 1.0);
    // If we are background subtracting a spectrum, we can end up with bins with really
    //  small values, like 0.0007, which, even one of would mess the whole fit up if
    //  we take its uncertainty to be its square-root, so in this case we will, fairly arbitrarily
    //  we want to use an uncert of 1.
    //  However, there are also highly scaled spectra, whose all values are really small, so
    //  in this case we want to do something more reasonable.
    // TODO: evaluate these choices of thresholds and tradeoffs, more better
    ScalarType uncert = variances ? ScalarType( static_cast<double>(sqrt(variances[row])) )
                                  : ((dataval > PEAK_FIT_MIN_CHANNEL_UNCERT) ? sqrt(dataval) : ScalarType(1.0));
    assert( !variances || (variances[row] > 0.0f) );

    uncerts(row) = uncert;

    if( step_continuum && !cdf_step )
      step_cumulative_data += (dataval - ScalarType( min_data_val ));

    if( nfixedpeak )
    {
      assert( fixed_peak_contrib.size() == nbin );
      dataval -= fixed_peak_contrib[row];
    }

    if( !fixed_step_contrib.empty() )
      dataval -= fixed_step_contrib[row];

    y(row) = ((dataval > 0.0 ? dataval : ScalarType(0.0)) / uncert);

    for( Eigen::Index col = 0; col < num_poly_terms; ++col )
    {
      const ScalarType exp = ScalarType(col + 1.0);

      if( !cdf_step && step_continuum
         && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
         && (col == (num_polynomial_terms - 1)) )
      {
        // This logic mirrors that of PeakContinuum::offset_integral(...), and code
        // If you change it in one place - change it in here, below, and in offset_integral.
        const ScalarType frac_data = (roi_data_sum > 0.0)
            ? (step_cumulative_data - ScalarType( 0.5*(data[row] - min_data_val) )) / roi_data_sum
            : ScalarType( 0.5 );
        const ScalarType contribution = frac_data * (x1 - x0);

        A(row,col) = ScalarType( contribution / uncert );
      }else if( !cdf_step && step_continuum && (num_polynomial_terms == 4) )
      {
        const ScalarType frac_data = (roi_data_sum > 0.0)
            ? (step_cumulative_data - ScalarType( 0.5*(data[row] - min_data_val) )) / roi_data_sum
            : ScalarType( 0.5 );

        ScalarType contrib( 0.0 );
        switch( col )
        {
          case 0: contrib = (1.0 - frac_data) * (x1_rel - x0_rel);                     break;
          case 1: contrib = 0.5 * (1.0 - frac_data) * (x1_rel*x1_rel - x0_rel*x0_rel); break;
          case 2: contrib = frac_data * (x1_rel - x0_rel);                             break;
          case 3: contrib = 0.5 * frac_data * (x1_rel*x1_rel - x0_rel*x0_rel);         break;
          default: assert( 0 ); break;
        }//switch( col )

        A(row,col) = contrib / uncert;
      }else
      {
        const ScalarType contribution = (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));

        A(row,col) = contribution / uncert;
      }
    }//for( int order = 0; order < maxorder; ++order )
  }//for( size_t row = 0; row < nbin; ++row )

  //  TODO: multithread if we have more than X peaks, etc.
  //
  // For convienience, we'll keep peak unit-area counts around, but
  //  we could just use `A`, i.e., `unit_peak_counts[i][bin] == A(bin,num_poly_terms + i)* uncerts(bin)`.
  vector<vector<ScalarType>> unit_peak_counts( npeaks, vector<ScalarType>(nbin,ScalarType(0.0)) );

  // For the CDF step types, peak i's basis function carries the step term
  //  `(SUM_k s_k*E'^k) * CDFbar_i(center) * dx` in addition to its own area - the step is bilinear
  //  with the amplitude, so it rides inside peak i's column and the LLS solves both at once.
  //  We keep the step part separately so the fitted continuum can be clamped at zero the same way
  //  `PeakContinuum::offset_integral(...)` clamps it (the two must agree channel-by-channel).
  vector<vector<ScalarType>> peak_step_counts;
  vector<ScalarType> cdf_at_centers( cdf_step ? nbin : size_t(0) );
  if( cdf_step )
    peak_step_counts.assign( npeaks, vector<ScalarType>(nbin,ScalarType(0.0)) );

  for( size_t i = 0; i < npeaks; ++i )
  {
    ScalarType *peak_areas = &(unit_peak_counts[i][0]);
    PeakDists::photopeak_function_integral( means[i], sigmas[i], ScalarType(1.0),
                                             skew_type, skew_parameters,
                                             nbin, x, peak_areas );

    if( cdf_step )
    {
      // Cumulative-summing the unit-area channel integrals gives the ROI-anchored CDF; it starts
      //  at the ROI's first channel edge and saturates at its last, matching the evaluator's
      //  `clamp(CDF, F0, F1) - F0` - see `PeakContinuum::cdf_step_anchor_energies(...)`.
      unit_pdf_to_cdf( peak_areas, cdf_at_centers.data(), nbin );

      for( size_t channel = 0; channel < nbin; ++channel )
      {
        const ScalarType dx = ScalarType( x[channel+1] - x[channel] );

        peak_step_counts[i][channel] = step_coeff_for_channel(channel) * cdf_at_centers[channel] * dx;
        peak_areas[channel] += peak_step_counts[i][channel];
        A(channel,num_poly_terms + i) = peak_areas[channel] / uncerts(channel);
      }
    }else
    {
      for( size_t channel = 0; channel < nbin; ++channel )
        A(channel,num_poly_terms + i) = peak_areas[channel] / uncerts(channel);
    }
  }//for( size_t i = 0; i < npeaks; ++i )




  try
  {
    // The SVD seem to takeup something like 50% of the time creating peaks for `PeakFitLM::fit_peaks_LM(...)`,
    //  There are other faster methods below, but SVD does the best job.
    //
    //  Comparison over a large number of HPGe spectra searching for peaks, and using `PeakFitLM::fit_peaks_LM(...)`
    //    to fit candidates:
    // ----------------------------------------------------------------------------------------------------------------
    //  Method                            CPU Time        Score (lower is better - derived from success rate of fitting truth-known peaks, not accounting for area accuracy)
    //  BDCSVD                            40487 s         9.09444
    //  ColPivHouseholderQR               27316 s         9.1103
    //  PartialPivLU                      25234 s         9.10558
    //  LLT                               24900 s         9.10558
    //  CompleteOrthogonalDecomposition   27552 s         9.1103
    //  Minuit2-based methods               975 s         11.3892
    // -----------------------------------------------------------------------------------------------------------------
    //  (the required function_tolerance was set to 1E-9 I think, or maybe 1E-11)
    //
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
    const Eigen::JacobiSVD<Eigen::MatrixX<ScalarType>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
    const Eigen::BDCSVD<Eigen::MatrixX<ScalarType>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
    const Eigen::VectorX<ScalarType> coeffs = svd.solve(y); // coeffs will contain [c_0, c_1, ..., PeakAmp0, ....]
    // Or:
    //const Eigen::ColPivHouseholderQR<Eigen::MatrixX<ScalarType>> qr(A);
    //const Eigen::VectorX<ScalarType> coeffs = qr.solve(y);
    // Or:
    //const Eigen::PartialPivLU<Eigen::MatrixX<ScalarType>> lu(A);
    //const Eigen::VectorX<ScalarType> coeffs = lu.solve(y);
    // Or (if A is positive definite):
    //const Eigen::MatrixX<ScalarType> AtA = A.transpose() * A;
    //const Eigen::VectorX<ScalarType> Atb = A.transpose() * y;
    //const Eigen::LLT<Eigen::MatrixX<ScalarType>> llt(AtA);
    //const Eigen::VectorX<ScalarType> coeffs = llt.solve(Atb);
    // Or (useful when matrix is rank-defiecient):
    //const Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixX<ScalarType>> cod(A);
    //const Eigen::VectorX<ScalarType> coeffs = cod.solve(y);

    check_jet_array_for_NaN( coeffs.data(), coeffs.size() );

    assert( coeffs.size() == (num_poly_terms + npeaks) );
    assert( coeffs.rows() == (num_poly_terms + npeaks) );


    // Compute the covariance matrix of the fit parameters
    const Eigen::MatrixX<ScalarType> covariance_matrix = (A.transpose() * A).inverse();

    // Extract the uncertainties (standard deviations) from the diagonal of the covariance matrix
    const Eigen::VectorX<ScalarType> uncertainties = covariance_matrix.diagonal().array().sqrt();

    continuum_coeffs.resize( num_poly_terms );
    continuum_coeffs_uncerts.resize( num_poly_terms );
    for( size_t coef = 0; coef < num_poly_terms; ++coef )
    {
      continuum_coeffs[coef] = coeffs(coef);
      continuum_coeffs_uncerts[coef] = uncertainties[coef];
    }//for( int coef = 0; coef < poly_terms; ++coef )

    amplitudes.resize( npeaks );
    amplitudes_uncerts.resize( npeaks );

    for( size_t i = 0; i < npeaks; ++i )
    {
      const size_t coef = num_poly_terms + i;
      amplitudes[i] = coeffs(coef);
      amplitudes_uncerts[i] = uncertainties[coef];
    }//for( size_t i = 0; i < npeaks; ++i )

    ScalarType chi2( 0.0 );
    for( size_t bin = 0; bin < nbin; ++bin )
    {
      // The continuum is the polynomial plus, for the CDF step types, the step term that rides
      //  inside the peak columns.  Both belong to the continuum, so both must be inside the
      //  non-negativity clamp - `PeakContinuum::offset_integral(...)` clamps `poly + step`, and the
      //  model fitted here has to be the one InterSpec draws and integrates.
      //  NOTE: smoothing this hard clamp was tried (2026-07) and REVERTED - see the note in
      //  `fit_continuum(...)`.  Keep it hard; just apply it to the right quantity.
      ScalarType continuum_bin( 0.0 );
      for( size_t col = 0; col < num_poly_terms; ++col )
        continuum_bin += coeffs(col) * A(bin,col) * uncerts(bin);

      if( cdf_step )
      {
        for( size_t i = 0; i < npeaks; ++i )
          continuum_bin += coeffs(num_poly_terms + i) * peak_step_counts[i][bin];

        if( !fixed_step_contrib.empty() )
          continuum_bin += fixed_step_contrib[bin];
      }

      if( continuum_bin < 0.0 )
        continuum_bin = ScalarType(0.0);

      ScalarType y_pred = continuum_bin;

      for( size_t i = 0; i < npeaks; ++i )
      {
        const size_t col = num_poly_terms + i;

        // Gaussian part only; the step part of this column is already in `continuum_bin`.
        y_pred += coeffs(col) * (cdf_step ? (unit_peak_counts[i][bin] - peak_step_counts[i][bin])
                                          : unit_peak_counts[i][bin]);

        // We could get rid of keeping `unit_peak_counts[][]` around, as `A` has this same info.
        // For CDF step types, values can be large (step_coeff * CDF * dx), so use relative tolerance.
        assert( abs(unit_peak_counts[i][bin] - A(bin,num_poly_terms + i)* uncerts(bin))
                < (std::max)( ScalarType(1.0E-4), ScalarType(1.0E-6) * abs(unit_peak_counts[i][bin]) ) );
      }

      if( nfixedpeak )
      {
        assert( fixed_peak_contrib.size() == nbin );
        y_pred += fixed_peak_contrib[bin];
      }

      if( peak_counts )
        peak_counts[bin] += y_pred;

      chi2 += pow( (y_pred - static_cast<double>(data[bin])) / uncerts(bin), ScalarType(2.0) );
    }//for( int bin = 0; bin < nbin; ++bin )

    return chi2;
  }catch( std::exception &e )
  {

    cerr << "RelActCalcAuto::fut_amp_and_offset_imp(...): caught: " << e.what() << endl;

#ifndef NDEBUG
    if constexpr ( std::is_same_v<ScalarType, double> )
    {
      cerr << "fit_amp_and_offset(...): caught: " << e.what() << endl;
      cerr << "For means = {";
      for( double m : means )
        cerr << m << ", ";
      cerr << "}, sigmas={";
      for( double m : sigmas )
        cerr << m << ", ";
      cerr << "}" << endl;


      printf( "\nA=\n" );

      for( size_t row = 0; row < nbin; ++row )
      {
        for( size_t col = 0; col < nfit_terms; ++col )
          printf( "%12.2f, ", A(row,col) );
        printf( "\n" );
      }
      printf( "\n\n\n" );

      cerr << "For means = {";
      for( double m : means )
        cerr << m << ", ";
      cerr << "}, sigmas={";
      for( double m : sigmas )
        cerr << m << ", ";
      cerr << "}" << endl;
    }//if constexpr ( std::is_same_v<ScalarType, double> )
#endif //#ifndef NDEBUG

    throw runtime_error( "RelActCalcAuto::fit_continuum(...): trouble finding coeffs and amplitudes." );
  }//try / catch

  assert( 0 );
  throw std::logic_error( "shouldnt get here" );
}//double fit_amp_and_offset(...)


} //namespace PeakFit


#endif //PeakFit_imp_h
