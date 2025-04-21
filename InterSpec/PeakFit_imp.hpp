#ifndef PeakFit_imp_h
#define PeakFit_imp_h

#include <vector>
#include <exception>

#include "Eigen/Dense"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"

namespace PeakFit
{

/** This function fits the polynomial continuum for a region with a number of fixed amplitude peaks.
 *
 * Note that when ScalarType is a `ceres::Jet<>`, and peaks SkewType is not PeakDef::SkewType::NoSkew,
 * then some of the derivative values may become NaN during the `svd.solve(y)` step
 */
template<typename PeakType, typename ScalarType>
void fit_continuum( const float *x, const float *data, const size_t nbin,
                   const int num_polynomial_terms,
                   const bool step_continuum,
                   const ScalarType ref_energy,
                   const std::vector<PeakType> &fixedAmpPeaks,
                   const bool multithread,
                   ScalarType *continuum_coeffs,
                   ScalarType *peak_counts )
{
  using namespace std;

  if( step_continuum && ((num_polynomial_terms < 2) || (num_polynomial_terms > 4)) )
    throw std::runtime_error( "fit_continuum: Only 2 to 4 terms are supported for step continuums" );

  if( num_polynomial_terms < 0 )
    throw std::runtime_error( "fit_continuum: continuum must have at least 0 (e.g., no continuum) terms" );

  if( num_polynomial_terms > 4 )
    throw std::runtime_error( "fit_continuum: you asked for a higher order polynomial continuum than reasonable" );

  // Loosely following:
  //   https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  // using the SVD (slowest, but most accurate) method

  const Eigen::Index num_poly_terms = static_cast<Eigen::Index>( num_polynomial_terms  );

  Eigen::MatrixX<ScalarType> A( static_cast<Eigen::Index>(nbin), num_poly_terms );
  Eigen::VectorX<ScalarType> y( static_cast<Eigen::Index>(nbin) );
  Eigen::VectorX<ScalarType> uncerts( static_cast<Eigen::Index>(nbin) );

  double roi_data_sum = 0.0, step_cumulative_data = 0.0;
  for( size_t row = 0; row < nbin; ++row )
    roi_data_sum += std::max( data[row], 0.0f );

  // Zero out the destination count array
  for( size_t row = 0; row < nbin; ++row )
    peak_counts[row] = ScalarType(0.0);

  // Add the Gaussian + Skew component of the peaks to destination counts
  const size_t nfixedpeak = fixedAmpPeaks.size();

  if( multithread && (nfixedpeak > 8) ) //8 is arbitrary.
  {
    // TODO: multi-thread computation needs to be evaluated more hollistically both here and in #RelActAutoSolution::eval
    const unsigned nthread = std::min( 16, std::min( static_cast<int>(nfixedpeak), std::max( 1, static_cast<int>( std::thread::hardware_concurrency() ) ) ) );

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

  for( size_t row = 0; row < nbin; ++row )
  {
    const double data_counts = data[row];
    const double x0 = x[row];
    const double x1 = x[row+1];
    const ScalarType x0_rel = x0 - ref_energy;
    const ScalarType x1_rel = x1 - ref_energy;

    const double uncert = (data_counts > PEAK_FIT_MIN_CHANNEL_UNCERT) ? sqrt(data_counts) : 1.0;

    uncerts(row) = ScalarType(uncert);

    if( step_continuum )
      step_cumulative_data += data_counts;

    y(row) = (data_counts > peak_counts[row]) ? (data_counts - peak_counts[row])/uncert : ScalarType(0.0);

    for( Eigen::Index col = 0; col < num_poly_terms; ++col )
    {
      const double exp = col + 1.0;

      if( step_continuum
         && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
         && (col == (num_polynomial_terms - 1)) )
      {
        // This logic mirrors that of PeakContinuum::offset_integral(...), and code
        // If you change it in one place - change it in here, below, and in offset_integral.
        const double frac_data = (step_cumulative_data - 0.5*data_counts) / roi_data_sum;
        const double contribution = frac_data * (x1 - x0);

        A(row,col) = ScalarType(contribution / uncert);

        check_jet_for_NaN( A(row,col) );
      }else if( step_continuum && (num_polynomial_terms == 4) )
      {
        const double frac_data = (step_cumulative_data - 0.5*data_counts) / roi_data_sum;

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

        check_jet_for_NaN( contrib );
        check_jet_for_NaN( uncert );
        check_jet_for_NaN( A(row,col) );
      }else
      {
        const ScalarType contribution = (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));

        A(row,col) = contribution / uncert;

        check_jet_for_NaN( contribution );
        check_jet_for_NaN( uncert );
        check_jet_for_NaN( A(row,col) );
      }
    }//for( int order = 0; order < maxorder; ++order )
  }//for( size_t row = 0; row < nbin; ++row )


  try
  {
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
    const Eigen::JacobiSVD<Eigen::MatrixX<ScalarType>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
    const Eigen::BDCSVD<Eigen::MatrixX<ScalarType>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif

    const Eigen::VectorX<ScalarType> coeffs = svd.solve(y); // coeffs will contain [c_0, c_1, c_2, c_3]

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

    assert( coeffs.size() == num_poly_terms );
    assert( coeffs.rows() == num_poly_terms );

    for( Eigen::Index i = 0; i < num_poly_terms; ++i )
      continuum_coeffs[i] = coeffs(i);

    const Eigen::VectorX<ScalarType> cont_vals = (A * coeffs).array() * uncerts.array();

    for( size_t bin = 0; bin < nbin; ++bin )
    {
      ScalarType y_continuum = cont_vals(bin);
      if( y_continuum < 0.0 )
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
 @param nbin The number of channels in the ROI.
 @param num_polynomial_terms The number of polynomial continuum terms to fit for.
        0 is no continuum (untested), 1 is constant, 2 is linear sloped continuum, etc
 @param step_continuum Specifies whether or not a step in the continuum is to be used.
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
ScalarType fit_amp_and_offset_imp( const float *x, const float *data, const size_t nbin,
                          const int num_polynomial_terms,
                          const bool step_continuum,
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

  if( sigmas.size() != means.size() )
    throw runtime_error( "fit_amp_and_offset_imp: invalid input" );

  if( step_continuum && ((num_polynomial_terms < 2) || (num_polynomial_terms > 4)) )
    throw runtime_error( "fit_amp_and_offset_imp: Only 2 to 4 terms are supported for step continuums" );

  if( num_polynomial_terms < 0 )
    throw runtime_error( "fit_amp_and_offset_imp: continuum must have at least 0 (e.g., no continuum) terms" );

  if( num_polynomial_terms > 4 )
    throw runtime_error( "fit_amp_and_offset_imp: you asked for a higher order polynomial continuum than reasonable" );

  assert( (skew_type == PeakDef::SkewType::NoSkew) || skew_parameters );
  if( !skew_parameters && (skew_type != PeakDef::SkewType::NoSkew) )
    throw std::logic_error( "fit_amp_and_offset_imp: skew pars not provided" );


  using namespace boost::numeric;
  const size_t npeaks = sigmas.size();

  const Eigen::Index num_poly_terms = static_cast<Eigen::Index>( num_polynomial_terms  );
  const Eigen::Index nfit_terms = static_cast<int>( num_poly_terms + npeaks );

  Eigen::MatrixX<ScalarType> A( static_cast<Eigen::Index>(nbin), nfit_terms );
  Eigen::VectorX<ScalarType> y( static_cast<Eigen::Index>(nbin) );
  Eigen::VectorX<ScalarType> uncerts( static_cast<Eigen::Index>(nbin) );


  ScalarType roi_data_sum( 0.0 ), step_cumulative_data( 0.0 );

  const double roi_lower = x[0];
  const double roi_upper = x[nbin];

  for( size_t row = 0; row < nbin; ++row )
    roi_data_sum += std::max( static_cast<double>(data[row]), 0.0 );

  const ScalarType avrg_data_val = roi_data_sum / static_cast<double>(nbin);

  // TODO: implement multi-threaded computation
  // 20250127: Using SpecUtilsAsync::ThreadPool, it looks like calling `pool.join()` is causing significant and unreasonable delays
  //           on macOS (using GCD, at least).  This is likely a problem with
  //           `SpecUtilsAsync::ThreadPool` - I would guess when creating ThreadPools inside of
  //           other ThreadPools, but for the moment will just do this single threaded, which is
  //           like 20 times faster for an example problem

  const size_t nfixedpeak = fixedAmpPeaks.size();
  vector<ScalarType> fixed_peak_contrib( nfixedpeak ? nbin : size_t(0), ScalarType(0.0) );

  if( nfixedpeak )
  {
    ScalarType * const fixed_contrib = &(fixed_peak_contrib[0]);
    for( size_t peak_index = 0; peak_index < fixedAmpPeaks.size(); ++peak_index )
      fixedAmpPeaks[peak_index].gauss_integral( x, fixed_contrib, nbin );
  }//if( nfixedpeak )


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
    ScalarType uncert = (dataval > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt(dataval) : ScalarType(1.0));

    uncerts(row) = uncert;

    if( step_continuum )
      step_cumulative_data += dataval;

    if( nfixedpeak )
    {
      assert( fixed_peak_contrib.size() == nbin );
      dataval -= fixed_peak_contrib[row];
    }

    y(row) = ((dataval > 0.0 ? dataval : ScalarType(0.0)) / uncert);

    for( Eigen::Index col = 0; col < num_poly_terms; ++col )
    {
      const ScalarType exp = ScalarType(col + 1.0);

      if( step_continuum
         && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
         && (col == (num_polynomial_terms - 1)) )
      {
        // This logic mirrors that of PeakContinuum::offset_integral(...), and code
        // If you change it in one place - change it in here, below, and in offset_integral.
        const ScalarType frac_data = (step_cumulative_data - 0.5*data[row]) / roi_data_sum;
        const ScalarType contribution = frac_data * (x1 - x0);

        A(row,col) = ScalarType( contribution / uncert );
      }else if( step_continuum && (num_polynomial_terms == 4) )
      {
        const ScalarType frac_data = (step_cumulative_data - 0.5*data[row]) / roi_data_sum;

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

  for( size_t i = 0; i < npeaks; ++i )
  {
    ScalarType *peak_areas = &(unit_peak_counts[i][0]);
    PeakDists::photopeak_function_integral( means[i], sigmas[i], ScalarType(1.0),
                                             skew_type, skew_parameters,
                                             nbin, x, peak_areas );

    for( size_t channel = 0; channel < nbin; ++channel )
    {
      A(channel,num_poly_terms + i) = peak_areas[channel] / uncerts(channel);
    }//for( size_t channel = 0; channel < nbin; ++channel )
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
      ScalarType y_pred( 0.0 );
      for( size_t col = 0; col < num_poly_terms; ++col )
        y_pred += coeffs(col) * A(bin,col) * uncerts(bin);

      if( y_pred < 0.0 )
        y_pred = ScalarType(0.0);

      for( size_t i = 0; i < npeaks; ++i )
      {
        const size_t col = num_poly_terms + i;
        y_pred += coeffs(col) * unit_peak_counts[i][bin];

        // We could get rid of keeping `unit_peak_counts[][]` around, as `A` has this same info.
        assert( abs(unit_peak_counts[i][bin] - A(bin,num_poly_terms + i)* uncerts(bin)) < 1.0E-4 );
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
