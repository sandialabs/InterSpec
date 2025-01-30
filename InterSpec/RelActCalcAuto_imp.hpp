#include <vector>



#include "Eigen/Dense"


namespace RelActCalcAuto
{


/** This function fits the polynomial continuum for a region with a number of fixed amplitude peaks.
 * 
 * 
 */
template<typename PeakType, typename ScalarType>
void fit_continuum( const float *x, const float *data, const size_t nbin,
                                  const int num_polynomial_terms,
                                  const bool step_continuum,
                                  const ScalarType ref_energy,
                                  const std::vector<PeakType> &fixedAmpPeaks,
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
  
  for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
  {      
    fixedAmpPeaks[peak_index].gauss_integral( x, peak_counts, nbin );
  }//for( size_t peak_index = 0; peak_index < nfixedpeak; ++peak_index )
  
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double data_counts = data[row];
    const double x0 = x[row];
    const double x1 = x[row+1];
    const ScalarType x0_rel = x0 - ref_energy;
    const ScalarType x1_rel = x1 - ref_energy;
    
    const double uncert = (data_counts > MIN_CHANNEL_UNCERT) ? sqrt(data_counts) : 1.0;

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
      }else
      {
        const ScalarType contribution = (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));
        
        A(row,col) = contribution / uncert;
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
