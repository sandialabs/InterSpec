#include <vector>

#include <thread>

#include <boost/math/distributions/normal.hpp>

#include "Eigen/Dense"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"

namespace RelActCalcAuto
{
/** A stand-in for the `PeakDef` class to allow auto-differentiation, and also simplify things  */
template<typename T>
struct PeakDefImp
{
  T m_mean = T(0.0);
  T m_sigma = T(0.0);
  T m_amplitude = T(0.0);
  T m_skew_pars[4] = { T(0.0), T(0.0), T(0.0), T(0.0) };

  const SandiaDecay::Nuclide *m_parent_nuclide = nullptr;
  const SandiaDecay::Transition *m_transition = nullptr;
  size_t m_rad_particle_index = 0;

  const SandiaDecay::Element *m_xray_element = nullptr;
  const ReactionGamma::Reaction *m_reaction = nullptr;

  /** True energy of the source gamma or x-ray. */
  double m_src_energy = 0.0;

  PeakDef::SkewType m_skew_type = PeakDef::SkewType::NoSkew;
  PeakDef::SourceGammaType m_gamma_type = PeakDef::SourceGammaType::NormalGamma;

  size_t m_rel_eff_index = std::numeric_limits<size_t>::max();

  // Some functions to be compatible with PeakDef, in templated functions
  const T &mean() const { return m_mean; }
  const T &sigma() const { return m_sigma; }
  const T &amplitude() const { return m_amplitude; }

  void setMean( const T &mean ) { m_mean = mean; }
  void setSigma( const T &sigma ) { m_sigma = sigma; }
  void setAmplitude( const T &amp ) { m_amplitude = amp; }

  inline void setSkewType( const PeakDef::SkewType &type )
  {
    m_skew_type = type;
  }

  inline void set_coefficient( T val, const PeakDef::CoefficientType &coef )
  {
    const int index = static_cast<int>( coef - PeakDef::CoefficientType::SkewPar0 );
    if( index >= (sizeof(m_skew_pars) / sizeof(m_skew_pars[0])) )
    {
      assert( coef == PeakDef::CoefficientType::Chi2DOF );
      return; //Like Chi2Dof
    }

    assert( (index >= 0) && (index < 4) );
    m_skew_pars[index] = val;
  }

  void gauss_integral( const float *energies, T *channels, const size_t nchannel ) const
  {
    check_jet_array_for_NaN( channels, nchannel );

    switch( m_skew_type )
    {
      case PeakDef::SkewType::NoSkew:
        PeakDists::gaussian_integral( m_mean, m_sigma, m_amplitude, energies, channels, nchannel );
        break;

      case PeakDef::SkewType::Bortel:
        PeakDists::bortel_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::CrystalBall:
        PeakDists::crystal_ball_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::DoubleSidedCrystalBall:
        PeakDists::double_sided_crystal_ball_integral( m_mean, m_sigma, m_amplitude,
                                           m_skew_pars[0], m_skew_pars[1],
                                           m_skew_pars[2], m_skew_pars[3],
                                           energies, channels, nchannel );
        break;

      case PeakDef::SkewType::GaussExp:
        PeakDists::gauss_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], energies, channels, nchannel );
        break;

      case PeakDef::SkewType::ExpGaussExp:
        PeakDists::exp_gauss_exp_integral( m_mean, m_sigma, m_amplitude, m_skew_pars[0], m_skew_pars[1], energies, channels, nchannel );
        break;
    }//switch( skew_type )

    check_jet_array_for_NaN( channels, nchannel );
  }//void gauss_integral( const float *energies, double *channels, const size_t nchannel ) const


  /** Gives lower and upper energies that contain `1.0 - missing_frac` of the peak,

   For skewed distributions, particularly Crystal Bal, the energy range given by the desired
   peak coverage can get huge, so you can also specify the max number of FWHM to use, which will
   truncate the range in these cases of very large skew.

   @param missing_frac The fraction of the peaks area you are okay not including; half this
          amount will be missing from both lower and upper distribution tails.
          A typical value might be like 1.0E-4 to get 99.99% of the peak.
   @param max_num_fwhm If greater than zero, the energy range will be truncated to be this
          amount of FWHMs from the mean, if the coverage would have the energy range go really
          far out.
   */
  std::pair<double,double> peak_coverage_limits( const double missing_frac, const double max_num_fwhm )
  {
    using namespace std;

    double skew_pars[4] = { 0.0 };

    double mean, sigma;
    if constexpr ( !std::is_same_v<T, double> )
    {
      mean = m_mean.a;
      sigma = m_sigma.a;
    }else
    {
      mean = m_mean;
      sigma = m_sigma;
    }

    const size_t nskew_par = PeakDef::num_skew_parameters(m_skew_type);
    for( size_t i = 0; i < nskew_par; ++i )
    {
      if constexpr ( !std::is_same_v<T, double> )
        skew_pars[i] = m_skew_pars[i].a;
      else
        skew_pars[i] = m_skew_pars[i];
    }


    pair<double,double> vis_limits;

    switch( m_skew_type )
    {
      case PeakDef::SkewType::NoSkew:
      {
        const boost::math::normal_distribution gaus_dist( 1.0 );
        vis_limits.first = mean +  sigma*boost::math::quantile( gaus_dist, 0.5*missing_frac );
        vis_limits.second = mean + sigma*boost::math::quantile( gaus_dist, 1.0 - 0.5*missing_frac );
        break;
      }

      case PeakDef::SkewType::Bortel:
        vis_limits = PeakDists::bortel_coverage_limits( mean, sigma, skew_pars[0], missing_frac );
        break;

      case PeakDef::SkewType::GaussExp:
        vis_limits = PeakDists::gauss_exp_coverage_limits( mean, sigma, skew_pars[0], missing_frac );
        break;

      case PeakDef::SkewType::CrystalBall:
        try
      {
        vis_limits = PeakDists::crystal_ball_coverage_limits( mean, sigma, skew_pars[0], skew_pars[1], missing_frac );
      }catch( std::exception & )
      {
        // CB dist can have really long tail, causing the coverage limits to fail, because
        //  of unreasonable values - in this case we'll just go way out
        vis_limits.first = mean - 15.0*sigma;

        const boost::math::normal_distribution gaus_dist( 1.0 );
        vis_limits.second = mean + sigma*boost::math::quantile( gaus_dist, 1.0 - missing_frac );
      }
        break;

      case PeakDef::SkewType::ExpGaussExp:
        vis_limits = PeakDists::exp_gauss_exp_coverage_limits( mean, sigma, skew_pars[0],
                                                              skew_pars[1], missing_frac );
        break;

      case PeakDef::SkewType::DoubleSidedCrystalBall:
        try
      {
        vis_limits = PeakDists::double_sided_crystal_ball_coverage_limits( mean, sigma, skew_pars[0],
                                                                          skew_pars[1], skew_pars[2], skew_pars[3], missing_frac );
      }catch( std::exception & )
      {
        // CB dist can have really long tail, causing the coverage limits to fail, because
        //  of unreasonable values - in this case we'll just go way out
        vis_limits.first  = mean - 15.0*sigma;
        vis_limits.second = mean + 15.0*sigma;
      }//try / catch
        break;
    }//switch( m_skew_type )

    if( max_num_fwhm > 0.0 )
    {
      vis_limits.first = std::max( vis_limits.first, mean - 2.35482*max_num_fwhm*sigma );
      vis_limits.second = std::min( vis_limits.second, mean + 2.35482*max_num_fwhm*sigma );
    }//if( max_num_fwhm > 0.0 )

    return vis_limits;
  }//pair<double,double> peak_coverage_limits( const double missing_frac = 1.0E-6 )

};//struct PeakDefImp

template<typename T>
struct PeakContinuumImp
{
  PeakContinuum::OffsetType m_type = PeakContinuum::OffsetType::NoOffset;

  T m_lower_energy = T(0.0);
  T m_upper_energy = T(0.0);
  T m_reference_energy = T(0.0);
  std::array<T,5> m_values = { T(0.0), T(0.0), T(0.0), T(0.0), T(0.0) };

  // Some functions for compatibility with `PeakContinuum` in PeakDef.h.
  void setRange( const T lowerenergy, const T upperenergy )
  {
    m_lower_energy = lowerenergy;
    m_upper_energy = upperenergy;
  }

  T lowerEnergy() const { return m_lower_energy; }
  T upperEnergy() const { return m_upper_energy; }

  void setType( PeakContinuum::OffsetType type )
  {
    m_type = type;
  }

  PeakContinuum::OffsetType type() const
  {
    return m_type;
  }

  T referenceEnergy() const { return m_reference_energy; }
  const std::array<T,5> &parameters() const { return m_values; }

  void setParameters( T referenceEnergy, const T *parameters, const T *uncertainties [[maybe_unused]] )
  {
    m_reference_energy = referenceEnergy;
    const size_t npar = PeakContinuum::num_parameters(m_type);
    for( size_t i = 0; i < npar; ++i )
      m_values[i] = parameters[i];
  }

  std::shared_ptr<const SpecUtils::Measurement> externalContinuum() const { return nullptr; }
};//struct PeakContinuumImp


template<typename T>
struct PeaksForEnergyRangeImp
{
  std::vector<PeakDefImp<T>> peaks;

  PeakContinuumImp<T> continuum;

  size_t first_channel;
  size_t last_channel;
  bool no_gammas_in_range;
  bool forced_full_range;

  /** Peak plus continuum counts for [first_channel, last_channel] */
  std::vector<T> peak_counts;
};//struct PeaksForEnergyRangeImp



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
  
  static const double MIN_CHANNEL_UNCERT = 1.0;
  
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
