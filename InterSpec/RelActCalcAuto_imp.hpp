#ifndef RelActCalcAuto_imp_h
#define RelActCalcAuto_imp_h

#include <vector>

#include <thread>


#include "Eigen/Dense"


#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDists_imp.hpp"  //for `check_jet_for_NaN(...)`


namespace RelActCalcAuto
{


/** This function fits the polynomial continuum for a region with a number of fixed amplitude peaks.
 * 
 * Note that when ScalarType is a `ceres::Jet<>`, and peaks SkewType is not PeakDef::SkewType::NoSkew,
 * then some of the derivative values may become NaN during the `svd.solve(y)` step
 */
template<typename PeakType, typename ScalarType>
void fit_continuum( const float * const x, 
                    const float * const data, 
                    const float * const data_uncert,
                    const size_t nbin,
                    const int num_polynomial_terms,
                    const bool step_continuum,
                    const ScalarType ref_energy,
                    const std::vector<PeakType> &fixedAmpPeaks,
                    const bool multithread,
                    ScalarType *continuum_coeffs,
                    ScalarType *peak_counts )
{
  using namespace std;
  
  static const double MIN_CHANNEL_UNCERT = 1.0;
  
  if( step_continuum && ((num_polynomial_terms < 2) || (num_polynomial_terms > 4)) )
    throw std::runtime_error( "fit_amp_and_offset: Only 2 to 4 terms are supported for step continuums" );
  
  if( num_polynomial_terms < 0 )
    throw std::runtime_error( "fit_amp_and_offset: continuum must have at least 0 (e.g., no continuum) terms" );
  
  if( num_polynomial_terms > 4 )
    throw std::runtime_error( "fit_amp_and_offset: you asked for a higher order polynomial continuum than reasonable" );
  
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
    const double data_counts_uncert = data_uncert ? data_uncert[row] : sqrt(data_counts);
    const double x0 = x[row];
    const double x1 = x[row+1];
    const ScalarType x0_rel = x0 - ref_energy;
    const ScalarType x1_rel = x1 - ref_energy;
    
    const double uncert = (data_counts_uncert > MIN_CHANNEL_UNCERT) ? data_counts_uncert : 1.0;

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

}//namespace RelActCalcAuto

#endif //RelActCalcAuto_imp_h
