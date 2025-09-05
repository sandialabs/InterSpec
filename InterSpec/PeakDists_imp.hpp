#ifndef PeakDists_imp_h
#define PeakDists_imp_h

#include <vector>
#include <numeric>
#include <algorithm>
#include <exception>

#include <boost/math/constants/constants.hpp>
#include <unsupported/Eigen/SpecialFunctions>

#include "SpecUtils/SpecFile.h" //Needed for `offset_integral(...)`

// Had a little trouble with the auto-derivative when using Jet - so will define some functions
//  here to help find the issues - but make them be no-ops for non-debug builds
#if( !defined(NDEBUG) && PERFORM_DEVELOPER_CHECKS && defined(CERES_PUBLIC_JET_H_) )
inline void check_jet_for_NaN( const double &jet )
{
  //no-op
}

inline void check_jet_array_for_NaN( const double * const jet, const size_t nelements )
{
  //no-op
}

template<typename T, int N>
inline void check_jet_for_NaN( const ceres::Jet<T,N> &jet )
{
  using namespace std;

  const Eigen::Matrix<double, N, 1> &matrix = jet.v;
  for( size_t par = 0; par < matrix.size(); ++par )
  {
    const double *vals = matrix.data();
    const double &val = vals[par];
    if( isnan(val) || isinf(val) )
    {
      cerr << "For par " << par << " val=" << val << endl;
      assert( !isnan(val) );
      assert( !isinf(val) );
      //val = 0.0;
    }
  }
}//void check_jet_for_NaN( ceres::Jet<T,N> &jet )

template<typename T, int N>
inline void check_jet_array_for_NaN( const ceres::Jet<T,N> * const jets, const size_t nelements )
{
  for( size_t i = 0; i < nelements; ++i )
    check_jet_for_NaN( jets[i] );
}
#else

#define check_jet_for_NaN(a){}
#define check_jet_array_for_NaN(a,b){}

#endif



/** Templating the peak distributions on calculation type is for `RelActAuto`, so they can be computed with
 `ceres::Jet<>` instead of `double`, to allow automatic differentiation.
 */
namespace PeakDists
{
    
template<typename T>
void gaussian_integral( const T peak_mean,
                              const T peak_sigma,
                              const T peak_amplitude,
                              const float * const energies,
                              T *channels,
                              const size_t nchannel )
{
  if( peak_sigma == 0.0 )
    return;

  if constexpr ( std::is_same_v<T, double> )
  {
    // We dont want to return for zero-amplitude peaks if we are usign a Ceres::Jet,
    //  as we may need to take into account the derivatives at zero
    if( peak_amplitude == 0.0 )
      return;
  }//if constexpr ( std::is_same_v<T, double> )

  const double zero_amp_point_nsigma = 8.0;
  const T start_energy = peak_mean - zero_amp_point_nsigma*peak_sigma;
  const T stop_energy = peak_mean + zero_amp_point_nsigma*peak_sigma;
  
  size_t channel = 0;
  while( (channel < nchannel) && (static_cast<double>(energies[channel+1]) < start_energy) )
  {
    channel += 1;
  }
    
  if( channel == nchannel )
    return;
    
  const double sqrt2 = boost::math::constants::root_two<double>();
  const T sqrt2sigma = sqrt2 * peak_sigma;
  const T amp_mult = 0.5 * peak_amplitude;


  // TODO: it looks like Eigen might support vectorized `erf`; should update to take advantage of.
  //       In which case we could maybe vectorize this function.  However, on Apple NEON there is
  //       apparently a bug (although maybe I'm mis-reading the Eigen source code), so it isnt
  //       supported there.  But for for x64, should try and see if this could speed things up.

  // We will keep track of the channels lower value of erf, so we dont have to re-compute
  //  it for each channel (this is the who advantage of )
  T erfarg = (static_cast<double>(energies[channel]) - peak_mean) / sqrt2sigma;


  // We will keep track of the channels lower value of erf, so we dont have to re-compute
  //  it for each channel (this is the who advantage of )
  T erflow;
  if constexpr ( !std::is_same_v<T, double> )
  {
    erflow = erf( erfarg );
  }else
  {
    // `Eigen::numext::erf<float>(erfarg)` is a little faster than `boost_erf_imp(erfarg)`,
    //    but a little less accurate (it even clamps values to +-1 outside of fabs(erfarg) of 4),
    //    so we'll only use it where we wont notice this slight reduction of accuracy
    erflow = (amp_mult > 1.0E5)
                          ? boost_erf_imp(erfarg)
                          : static_cast<double>( Eigen::numext::erf(static_cast<float>(erfarg)) );
  }

#ifndef NDEBUG
  double boost_erflow, eigen_erflow;
  if constexpr ( std::is_same_v<T, double> )
  {
    // We'll check that things area as accurate as we expect
    // Should maybe enable this check using PERFORM_DEVELOPER_CHECKS
    boost_erflow = boost_erf_imp(erfarg);
    eigen_erflow = Eigen::numext::erf(static_cast<float>(erfarg));
  }
#endif

  while( (channel < nchannel) && (energies[channel] < stop_energy) )
  {
    erfarg = (static_cast<double>(energies[channel+1]) - peak_mean)/(sqrt2*peak_sigma);

#ifndef NDEBUG
    if constexpr ( std::is_same_v<T, double> )
    {
      const double boost_erfhigh = boost_erf_imp( erfarg );
      const double eigen_erfhigh = Eigen::numext::erf(static_cast<float>(erfarg));
      const double boost_val = amp_mult * (boost_erfhigh - boost_erflow);
      const double eigen_val = amp_mult * (eigen_erfhigh - eigen_erflow);

      assert( fabs(boost_val - eigen_val) < 1.0E-5*abs(amp_mult) || fabs(boost_val - eigen_val) < 0.001 );

      boost_erflow = boost_erfhigh;
      eigen_erflow = eigen_erfhigh;
    }//if constexpr ( std::is_same_v<T, double> )
#endif

    T erfhigh;
    if constexpr ( !std::is_same_v<T, double> )
    {
      erfhigh = erf( erfarg );
    }else
    {
      erfhigh = (amp_mult > 1.0E5) ? boost_erf_imp( erfarg )
                                   : static_cast<double>( Eigen::numext::erf(static_cast<float>(erfarg)) );
    }

    channels[channel] += amp_mult * (erfhigh - erflow);
    channel += 1;
    erflow = erfhigh;
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
}//gaus_integral(...)


/** Returns the normalization so the GaussExp distribution has unit area. */
template<typename T>
T gauss_exp_norm( const T sigma, const T skew )
{
  static const double sqrt_pi = boost::math::constants::root_pi<double>(); //1.7724538509055160272981
  static const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.707106781186547524400
    
  if constexpr ( !std::is_same_v<T, double> )
  {
    return 1.0 / ((sigma/skew)*exp(-0.5*skew*skew)
                  + (sqrt_pi*one_div_root_two*(erf(skew*one_div_root_two)+1.0)*sigma));
  }else
  {
    return 1.0 / ((sigma/skew)*std::exp(-0.5*skew*skew)
                  + (sqrt_pi*one_div_root_two*(boost_erf_imp(skew*one_div_root_two)+1)*sigma));
  }
}
  
  
template<typename T>
void gauss_exp_integral( const T peak_mean,
                          const T peak_sigma,
                          const T peak_amplitude,
                          const T skew,
                          const float * const energies,
                          T *channels,
                          const size_t nchannel )
{
  using namespace std;
  
#define USE_SIMPLE_GAUSS_EXP_IMP 0
    
#if( USE_SIMPLE_GAUSS_EXP_IMP )
    //compiled in debug mode, this implementation takes about 5 times as long as the more optimized version.
    
#ifdef _MSC_VER
#pragma message( "PeakDef::gauss_exp_integral is not properly coded" )
#else
#warning "PeakDef::gauss_exp_integral is not properly coded"
#endif
    
#if( PERFORM_DEVELOPER_CHECKS )
  T dist_sum = 0.0;
#endif
    
    for( size_t i = 0; i < nchannel; ++i )
    {
      const float x0 = energies[i];
      const float x1 = energies[i+1];
      const double val = peak_amplitude*PeakDists::gauss_exp_integral( peak_mean, peak_sigma, skew, x0, x1 );
      channels[i] += val;
      
#if( PERFORM_DEVELOPER_CHECKS )
      dist_sum += val;
      
      if( isinf(channels[i]) || isnan(channels[i]) )
      {
        cerr << "Found GausExp invalid counts, " << channels[i] << " from [" << x0 << ", " << x1 << "]:\n"
        << "\t" << setw(14) << "range:" << "[" << energies[0] << ", " << energies[nchannel] << "]\n"
        << "\t" << setw(14) << "min_energy =" << energies[0] << "\n"
        << "\t" << setw(14) << "max_energy =" << energies[nchannel] << "\n"
        << "\t" << setw(14) << "nchannel =" << nchannel << "\n";
        if constexpr ( std::is_same_v<T, double> )
          cerr << "\t" << setw(14) << "mean =" << peak_mean << "\n"
          << "\t" << setw(14) << "sigma =" << peak_sigma << "\n"
          << "\t" << setw(14) << "amp =" << peak_amplitude << "\n"
          << "\t" << setw(14) << "skew =" << skew << "\n";
        
        cerr << endl << endl;
      }//if( skew_type_t != PeakDef::NoSkew )
      
      //log_developer_error( __func__, "Invalid CSS color called back " );
#endif //PERFORM_DEVELOPER_CHECKS
    }
    
#else  //USE_SIMPLE_GAUSS_EXP_IMP
    
  if( peak_sigma == 0.0 )
    return;

  if constexpr ( std::is_same_v<T, double> )
  {
    // We dont want to return for zero-amplitude peaks if we are usign a Ceres::Jet,
    //  as we may need to take into account the derivatives at zero
    if( peak_amplitude == 0.0 )
      return;
  }//if constexpr ( std::is_same_v<T, double> )


  // TODO: estimate where we should actually start and stop computing values for, using `gauss_exp_coverage_limits(...)`, but need to check if it actually saves time
  const double zero_amp_point_nsigma = 8.0;
  const T start_energy( static_cast<double>(energies[0]) ); //static_cast<float>( peak_mean - zero_amp_point_nsigma*peak_sigma );
  const T stop_energy = peak_mean + zero_amp_point_nsigma*peak_sigma;
  
  size_t channel = 0;
  while( (channel < nchannel) && (static_cast<double>(energies[channel+1]) < start_energy) )
  {
    channel += 1;
  }
    
  if( channel == nchannel )
    return;
    
    
  const auto tail_indefinite_non_norm = [&peak_mean,&peak_sigma,&skew]( const T x ) -> T {
      const T t = (x - peak_mean) / peak_sigma;
      assert( (t - 1.0E-8) <= -skew );
      return (peak_sigma/skew)*exp((skew/peak_sigma)*(0.5*skew*peak_sigma - peak_mean + x));
  };
    
  const auto gaus_indefinite_non_norm = [&peak_mean,&peak_sigma,&skew]( const T x ) -> T {
    const T t = (x - peak_mean) / peak_sigma;
    assert( (t + 1.0E-6) >= -skew );
      
    const double root_half_pi = boost::math::constants::root_half_pi<double>();
    static const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.707106781186547524400
      
    if constexpr ( !std::is_same_v<T, double> )
      return peak_sigma*root_half_pi * erf(t*one_div_root_two);
    else
      return peak_sigma*root_half_pi * boost_erf_imp(t*one_div_root_two);
  };
    
  const T norm = peak_amplitude * gauss_exp_norm( peak_sigma, skew );
  const T tail_end = peak_mean - peak_sigma*skew;
    
  T lower_energy;
  
  if( energies[channel] < tail_end )
  {
    lower_energy = T( static_cast<double>(energies[channel]) );
    T indefinite_low = tail_indefinite_non_norm( lower_energy );
      
    while( (channel < nchannel) && (energies[channel] < tail_end) )
    {
      const T upper_energy( static_cast<double>(energies[channel+1]) );
        
      T indefinite_high;
      if( upper_energy > tail_end )
      {
        indefinite_high = tail_indefinite_non_norm( tail_end );
        const T val = norm * (indefinite_high - indefinite_low);
        channels[channel] += val;
          
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double simple_answer = peak_amplitude * PeakDists::gauss_exp_integral( peak_mean, peak_sigma, skew, lower_energy, tail_end );
          const double diff = fabs(val - simple_answer);
          const double frac_diff = diff / std::max(val, simple_answer);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
        }
#endif
          
        lower_energy = tail_end;
        break;
      }else
      {
        indefinite_high = tail_indefinite_non_norm( upper_energy );
        const T val = norm * (indefinite_high - indefinite_low);
        channels[channel] += val;
        
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double simple_answer = peak_amplitude * PeakDists::gauss_exp_integral( peak_mean, peak_sigma, skew, lower_energy, upper_energy );
          const double diff = fabs(val - simple_answer);
          const double frac_diff = diff / std::max(val, simple_answer);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
        }
#endif
          
        lower_energy = upper_energy;
        indefinite_low = indefinite_high;
        channel += 1;
      }
    }//while( (channel < nchannel) && (energies[channel] < tail_end) )
  }//if( energies[channel] < tail_end )
    
  if( channel >= nchannel )
    return;
    
  assert( energies[channel+1] >= tail_end );
  lower_energy = max( T( static_cast<double>(energies[channel]) ), tail_end );
  T indefinite_low = gaus_indefinite_non_norm( lower_energy );
    
  while( (channel < nchannel) && (static_cast<double>(energies[channel]) < stop_energy) )
  {
    const T upper_energy( static_cast<double>(energies[channel+1]) );
    const T indefinite_high = gaus_indefinite_non_norm( upper_energy );
    const T val = norm * (indefinite_high - indefinite_low);
    channels[channel] += val;
    
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
    if constexpr ( std::is_same_v<T, double> )
    {
      const double simple_answer = peak_amplitude * PeakDists::gauss_exp_integral( peak_mean, peak_sigma, skew, lower_energy, upper_energy );
      const double diff = fabs(val - simple_answer);
      const double frac_diff = diff / std::max(val, simple_answer);
      assert( (frac_diff < 1.0E-5) || (diff < 1.0E-7) || (fabs(peak_mean - 0.5*(lower_energy + upper_energy)) > (5.0*peak_sigma)) );
    }
#endif
    
    lower_energy = upper_energy;
    indefinite_low = indefinite_high;
    channel += 1;
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
#endif  //USE_SIMPLE_GAUSS_EXP_IMP
}//void gauss_exp_integral(...)
  
  
template<typename T>
T bortel_indefinite_integral( const double x, const T mean, const T sigma, const T skew )
{
  const double one_div_root_two = boost::math::constants::one_div_root_two<double>();
  
  const T t = (x - mean) / sigma;
  const T erf_arg = one_div_root_two*t;
  const T exp_arg = sigma*(2.0*skew*t + sigma) / (2.0*skew*skew);
  const T erfc_arg = one_div_root_two*(t + (sigma/skew));
  
  if( (skew <= 0.0) || (exp_arg > 87.0) || (erfc_arg > 10.0) )
  {
    if constexpr ( !std::is_same_v<T, double> )
      return 0.5 * erf(erf_arg);
    else
      return 0.5 * boost_erf_imp(erf_arg);
  }
  
  if constexpr ( !std::is_same_v<T, double> )
    return 0.5*(erf(erf_arg) + (exp(exp_arg) * erfc(erfc_arg)));
  else
    return 0.5*(boost_erf_imp(erf_arg) + (std::exp(exp_arg) * boost_erfc_imp(erfc_arg)));
}//double bortel_indefinite
  
  
template<typename T>
void bortel_integral( const T mean, const T sigma, const T amp, const T skew,
                       const float * const energies, T *channels, const size_t nchannel )
{
  assert( sigma > 0.0 );
  if( (sigma <= 0.0) || (amp <= 0.0) || !nchannel )
    return;
  
  const double zero_amp_point_nsigma_lower = 12.0; // TODO: Use the skew to determine lower energy
  const double zero_amp_point_nsigma_upper = 8.0;
  const T start_energy = mean - zero_amp_point_nsigma_lower*sigma;
  const T stop_energy = mean + zero_amp_point_nsigma_upper*sigma;
  
  size_t channel = 0;
  while( (channel < nchannel) && (static_cast<double>(energies[channel+1]) < start_energy) )
  {
    channel += 1;
  }
  
  if( channel == nchannel )
    return;
  
  // We will keep track of the channels lower value indefinite integral, so we dont have to
  //  re-compute it for each channel
  T val_low = bortel_indefinite_integral( static_cast<double>(energies[channel]), mean, sigma, skew );
  
  while( (channel < nchannel) && (energies[channel] < stop_energy) )
  {
    const T val_high = bortel_indefinite_integral( static_cast<double>(energies[channel+1]), mean, sigma, skew );
    const T val = amp*(val_high - val_low);
    
    channels[channel] += val;
    val_low = val_high;
    channel += 1;
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
}//bortel_integral( to array values)
  
  
template<typename T>
void crystal_ball_integral( const T peak_mean,
                             const T peak_sigma,
                             const T peak_amplitude,
                             const T alpha,
                             const T power_law,
                             const float * const energies,
                             T *channels,
                             const size_t nchannel )
{
  using namespace std;
  
#if( PERFORM_DEVELOPER_CHECKS )
  T dist_sum( 0.0 );
#endif
  
  if( peak_sigma == 0.0 )
    return;

  if constexpr ( std::is_same_v<T, double> )
  {
    // We dont want to return for zero-amplitude peaks if we are usign a Ceres::Jet,
    //  as we may need to take into account the derivatives at zero
    if( peak_amplitude == 0.0 )
      return;
  }//if constexpr ( std::is_same_v<T, double> )

  
  // TODO: estimate where we should actually start and stop computing values for, using `crystal_ball_coverage_limits(...)`, but need to check if it actually saves time
  const double zero_amp_point_nsigma = 8.0;
  const T start_energy( energies[0] ); //static_cast<float>( peak_mean - zero_amp_point_nsigma*peak_sigma );
  const T stop_energy = peak_mean + zero_amp_point_nsigma*peak_sigma;
  
  size_t channel = 0;
  while( (channel < nchannel) && (static_cast<double>(energies[channel+1]) < start_energy) )
  {
    channel += 1;
  }
  
  if( channel == nchannel )
    return;
  
  
  const T exp_aa = exp(-0.5*alpha*alpha);
  const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.70710678118654752440
  const double root_half_pi = boost::math::constants::root_half_pi<double>();
  const double sqrt_2pi = boost::math::constants::root_two_pi<double>();
  
  const T A = pow(power_law/alpha, power_law) * exp_aa;
  const T B = (power_law / alpha) - alpha;
  const T C = (power_law / alpha) * (1.0/(power_law - 1.0)) * exp_aa;
  T D;
  if constexpr ( !std::is_same_v<T, double> )
    D = root_half_pi * (1.0 + erf( one_div_root_two * alpha ));
  else
    D = root_half_pi * (1.0 + boost_erf_imp( one_div_root_two * alpha ));
  const T N = 1.0 / (peak_sigma * (C + D));
  const T tail_amp = peak_amplitude * N * A * peak_sigma / (power_law - 1.0);
  const T gauss_indef_amp = 0.5 * peak_amplitude * sqrt_2pi / (C + D);

  check_jet_for_NaN( N );
  check_jet_for_NaN( tail_amp );
  check_jet_for_NaN( gauss_indef_amp );
  check_jet_for_NaN( A );
  check_jet_for_NaN( B );
  check_jet_for_NaN( C );
  check_jet_for_NaN( D );
  check_jet_for_NaN( exp_aa );
  check_jet_for_NaN( start_energy );
  check_jet_for_NaN( stop_energy );

  // Brief implementation of crystal_ball_tail_indefinite_t
  auto tail_indefinite = [peak_mean,peak_sigma,alpha,power_law,B,tail_amp]( const T x ) -> T {
    const T t = (x - peak_mean) / peak_sigma;
    assert( ((t - 1.0E-6) <= -alpha) && (alpha > 0.0) && (power_law > 1.0) );

    T answer = tail_amp * pow( B - t, 1.0 - power_law );

    check_jet_for_NaN( t );
    check_jet_for_NaN( answer );

    return answer;
  };
  
  auto gauss_indefinite = [peak_mean,peak_sigma,gauss_indef_amp,one_div_root_two]( const T x ) -> T {
    const T t = (x - peak_mean) / peak_sigma;
    
    if constexpr ( !std::is_same_v<T, double> )
      return gauss_indef_amp * erf( one_div_root_two * t );
    else
      return gauss_indef_amp * boost_erf_imp( one_div_root_two * t );
  };
  
  
  const T tail_end = peak_mean - peak_sigma*alpha;
  
  T lower_energy;
  if( static_cast<double>(energies[channel]) < tail_end )
  {
    lower_energy = T(static_cast<double>(energies[channel]));
    T indefinite_low = tail_indefinite( lower_energy );
    
    while( (channel < nchannel) && (energies[channel] < tail_end) )
    {
      assert( energies[channel] < energies[channel+1] );
      
      const T upper_energy( static_cast<double>(energies[channel+1]) );
      
      if( upper_energy > tail_end )
      {
        const T indefinite_high = tail_indefinite( tail_end );
        const T val = (indefinite_high - indefinite_low);

        check_jet_for_NaN( indefinite_high );
        check_jet_for_NaN( val );
        check_jet_for_NaN( channels[channel] );

        channels[channel] += val;

        check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS )
        dist_sum += val;
#endif
        lower_energy = upper_energy;
        break;
      }else
      {
        const T indefinite_high = tail_indefinite( T(upper_energy) );
        const T val = (indefinite_high - indefinite_low);

        check_jet_for_NaN( indefinite_high );
        check_jet_for_NaN( val );
        check_jet_for_NaN( channels[channel] );

        channels[channel] += val;

        check_jet_for_NaN( channels[channel] );

        indefinite_low = indefinite_high;
        channel += 1;
        lower_energy = upper_energy;
      }
    }//while( (channel < nchannel) && (energies[channel] < tail_end) )
  }//if( energies[channel] < tail_end )
  
  if( channel >= nchannel )
    return;
  
  assert( static_cast<double>(energies[channel+1]) >= tail_end );
  lower_energy = max( T( static_cast<double>(energies[channel]) ), tail_end );
  T indefinite_low = gauss_indefinite( lower_energy );

  check_jet_for_NaN( lower_energy );
  check_jet_for_NaN( indefinite_low );

  while( (channel < nchannel) && (static_cast<double>(energies[channel]) < stop_energy) )
  {
    assert( energies[channel] < energies[channel+1] );
    const T upper_energy( static_cast<double>(energies[channel+1]) );
    
    const T indefinite_high = gauss_indefinite( upper_energy );
    const T val = (indefinite_high - indefinite_low);

    check_jet_for_NaN( indefinite_high );
    check_jet_for_NaN( val );
    check_jet_for_NaN( upper_energy );
    check_jet_for_NaN( channels[channel] );

    channels[channel] += val;

    check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
    if constexpr ( std::is_same_v<T, double> )
    {
      const double simple_answer = peak_amplitude * PeakDists::crystal_ball_integral( peak_mean,
                                    peak_sigma, alpha, power_law, lower_energy, upper_energy );
      const double diff = fabs(val - simple_answer);
      const double frac_diff = diff / std::max(val, simple_answer);
      
      // We'll be pragmatic here; we wont ever care about 1E-9 counts, so we will use this as our
      //   limit, otherwise we will actually trigger asserts.
      // TODO: figure out what is causing the in-accuracy of this (probably just numeric accuracy?)
      assert( (frac_diff < 1.0E-4) || (diff < 1.0E-9) );
    }
#endif //PERFORM_DEVELOPER_CHECKS
    
    lower_energy = upper_energy;
    indefinite_low = indefinite_high;
    channel += 1;
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
  
#if( PERFORM_DEVELOPER_CHECKS )
  /*
   const double frac = dist_sum / peak_amplitude;
   cerr << "Crystal Ball sum over [" << energies[0] << "," << energies[nchannel] << "] is "
   << frac << " (should be near 1)" << endl;
   
   if( fabs(1 - frac) > 0.01 )
   {
   cerr << "alpha(x) -> " << (peak_mean - alpha*peak_sigma) << endl;
   for( size_t i = 0; i < nchannel; ++i )
   cout << setw(5) << i << setw(12) << std::fixed << setprecision(2) << energies[i]
   << setw(12) << std::scientific << (channels[i]/peak_amplitude) << endl;
   
   cerr << "Single step integral is "
   << PeakDef::crystal_ball_integral( peak_mean, peak_sigma, 1.0, alpha, power_law, energies[0], energies[nchannel] )
   << " and over all area: "
   << PeakDef::crystal_ball_integral( peak_mean, peak_sigma, 1.0, alpha, power_law, peak_mean - 50*peak_sigma, peak_mean + 20*peak_sigma )
   << endl;
   }//if( fabs(1 - frac) > 0.01 )
   */
#endif
}//crystal_ball_integral(...)
  
  
template<typename T>
T exp_gauss_exp_norm( const T sigma, const T skew_left, const T skew_right )
{
  using namespace std;
  static const double sqrt_pi = boost::math::constants::root_pi<double>(); //1.7724538509055160272981
  static const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.707106781186547524400
  
  if constexpr ( !std::is_same_v<T, double> )
    return 1.0 / ((sigma/skew_left)*exp(-0.5*skew_left*skew_left)
                  + (one_div_root_two*sqrt_pi*sigma*(erf(one_div_root_two*skew_right)+erf(one_div_root_two*skew_left)))
                  + (sigma/skew_right)*exp(-0.5*skew_right*skew_right) );
  else
    return 1.0 / ((sigma/skew_left)*exp(-0.5*skew_left*skew_left)
                  + (one_div_root_two*sqrt_pi*sigma*(boost_erf_imp(one_div_root_two*skew_right)+boost_erf_imp(one_div_root_two*skew_left)))
                  + (sigma/skew_right)*std::exp(-0.5*skew_right*skew_right) );
}//exp_gauss_exp_norm(...)
  
  
template<typename T>
void exp_gauss_exp_integral( const T peak_mean,
                              const T peak_sigma,
                              const T peak_amplitude,
                              const T skew_left,
                              const T skew_right,
                              const float * const energies,
                              T *channels,
                              const size_t nchannel )
{
  using namespace std;
  
  if( peak_sigma == 0.0 )
    return;

  if constexpr ( std::is_same_v<T, double> )
  {
    // We dont want to return for zero-amplitude peaks if we are usign a Ceres::Jet,
    //  as we may need to take into account the derivatives at zero
    if( peak_amplitude == 0.0 )
      return;
  }//if constexpr ( std::is_same_v<T, double> )


  // TODO: estimate where we should actually start and stop computing values for, using `exp_gauss_exp_coverage_limits(...)`, but need to check if it actually saves time
  //const double zero_amp_point_nsigma = 8.0;
  const T start_energy( static_cast<double>(energies[0]) ); //peak_mean - zero_amp_point_nsigma*peak_sigma;
  const T stop_energy( static_cast<double>(energies[nchannel]) ); //peak_mean + zero_amp_point_nsigma*peak_sigma;
  
  size_t channel = 0;
  while( (channel < nchannel) && (static_cast<double>(energies[channel+1]) < start_energy) )
  {
    channel += 1;
  }
  
  if( channel == nchannel )
    return;
  
  
  auto left_tail_indefinite_non_norm = [peak_mean,peak_sigma,skew_left]( const T x ) -> T {
    return (peak_sigma/skew_left)*exp((skew_left/peak_sigma)*(0.5*skew_left*peak_sigma - peak_mean + x));
  };
  
  // This next line looks suspect, with both a multiplication by 0.5, and division by 2.0, but
  //  the developers checks would catch it if it was a real issue, since the alternative
  //  calculation is fairly independent.
  const T r_const = (exp(-0.5*skew_right*skew_right/2.0)*peak_sigma)/skew_right;
  auto right_tail_indefinite_non_norm = [peak_mean,peak_sigma,skew_right,r_const]( const T x ) -> T {
    return (r_const-(peak_sigma*exp((skew_right*peak_mean)/peak_sigma-(x*skew_right)/peak_sigma+0.5*skew_right*skew_right))/skew_right);
  };
  
  auto gauss_indefinite_non_norm = [peak_mean,peak_sigma]( const T x ) -> T {
    static const double root_half_pi = boost::math::constants::root_half_pi<double>(); //1.2533141373155002512078826424
    static const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.707106781186547524400
    
    const T t = (x - peak_mean) / peak_sigma;
    
    if constexpr ( !std::is_same_v<T, double> )
      return peak_sigma * root_half_pi * erf( one_div_root_two*t );
    else
      return peak_sigma * root_half_pi * boost_erf_imp( one_div_root_two*t );
  };
  
  const T norm = peak_amplitude * exp_gauss_exp_norm( peak_sigma, skew_left, skew_right );
  
  const T left_tail_end = peak_mean - peak_sigma*skew_left;
  const T right_tail_start = peak_mean + peak_sigma*skew_right;
  
  T lower_energy;
  if( static_cast<double>(energies[channel]) < left_tail_end )
  {
    lower_energy = T(static_cast<double>(energies[channel]));
    
    T indefinite_low = left_tail_indefinite_non_norm( lower_energy );
    
    while( (channel < nchannel) && (static_cast<double>(energies[channel]) < left_tail_end) )
    {
      const T upper_energy( static_cast<double>(energies[channel+1]) ) ;
      
      if( upper_energy > left_tail_end )
      {
        const T indefinite_high = left_tail_indefinite_non_norm( left_tail_end );
        const T val = norm * (indefinite_high - indefinite_low);

        check_jet_for_NaN( indefinite_high );
        check_jet_for_NaN( val );
        check_jet_for_NaN( channels[channel] );

        channels[channel] += val;
        
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
      if constexpr ( std::is_same_v<T, double> )
      {
        const double simple_answer = peak_amplitude * PeakDists::exp_gauss_exp_integral( peak_mean, peak_sigma, skew_left, skew_right, lower_energy, left_tail_end );
        const double diff = fabs(val - simple_answer);
        const double frac_diff = diff / std::max(val, simple_answer);
        assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
      }
#endif
        
        lower_energy = left_tail_end;
        break;
      }else
      {
        const T indefinite_high = left_tail_indefinite_non_norm( upper_energy );
        const T val = norm * (indefinite_high - indefinite_low);;

        check_jet_for_NaN( indefinite_high );
        check_jet_for_NaN( val );
        check_jet_for_NaN( channels[channel] );

        channels[channel] += val;
        
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
      if constexpr ( std::is_same_v<T, double> )
      {
        const double simple_answer = peak_amplitude * PeakDists::exp_gauss_exp_integral( peak_mean, peak_sigma, skew_left, skew_right, lower_energy, upper_energy );
        const double diff = fabs(val - simple_answer);
        const double frac_diff = diff / std::max(val, simple_answer);
        assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
      }
#endif
        lower_energy = upper_energy;
        indefinite_low = indefinite_high;
        channel += 1;
      }
    }//while( (channel < nchannel) && (energies[channel] < tail_end) )
  }//if( energies[channel] < tail_end )
  
  if( channel >= nchannel )
    return;
  
  assert( static_cast<double>(energies[channel+1]) >= left_tail_end );
  lower_energy = max( T(static_cast<double>(energies[channel])), left_tail_end );
  T indefinite_low = gauss_indefinite_non_norm( lower_energy );
  
  while( (channel < nchannel) && (static_cast<double>(energies[channel]) < right_tail_start) )
  {
    const T upper_energy( static_cast<double>(energies[channel+1]) );
    
    if( upper_energy > right_tail_start )
    {
      const T indefinite_high = gauss_indefinite_non_norm( right_tail_start );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );
      check_jet_for_NaN( channels[channel] );

      channels[channel] += val;
      
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
      if constexpr ( std::is_same_v<T, double> )
      {
        const double simple_answer = peak_amplitude * PeakDists::exp_gauss_exp_integral( peak_mean, peak_sigma, skew_left, skew_right, lower_energy, right_tail_start );
        const double diff = fabs(val - simple_answer);
        const double frac_diff = diff / std::max(val, simple_answer);
        assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
      }
#endif
      
      lower_energy = right_tail_start;
      break;
    }else
    {
      const T indefinite_high = gauss_indefinite_non_norm( upper_energy );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );
      check_jet_for_NaN( channels[channel] );

      channels[channel] += val;
      
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
      if constexpr ( std::is_same_v<T, double> )
      {
        const double simple_answer = peak_amplitude * PeakDists::exp_gauss_exp_integral( peak_mean, peak_sigma, skew_left, skew_right, lower_energy, upper_energy );
        const double diff = fabs(val - simple_answer);
        const double frac_diff = diff / std::max(val, simple_answer);
        assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
      }
#endif
      
      lower_energy = upper_energy;
      indefinite_low = indefinite_high;
      channel += 1;
    }
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
  
  
  if( channel >= nchannel )
    return;
  
  assert( static_cast<double>(energies[channel+1]) >= right_tail_start );
  
  lower_energy = max( T(static_cast<double>(energies[channel])), right_tail_start);
  indefinite_low = right_tail_indefinite_non_norm( lower_energy );
  
  while( (channel < nchannel) && (static_cast<double>(energies[channel]) < stop_energy) )
  {
    const T upper_energy( static_cast<double>(energies[channel+1]) );
    
    if( upper_energy > stop_energy )
    {
      const T indefinite_high = right_tail_indefinite_non_norm( right_tail_start );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );
      check_jet_for_NaN( channels[channel] );

      channels[channel] += val;
      
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
      if constexpr ( std::is_same_v<T, double> )
      {
        const double simple_answer = peak_amplitude * PeakDists::exp_gauss_exp_integral( peak_mean, peak_sigma, skew_left, skew_right, lower_energy, right_tail_start );
        const double diff = fabs(val - simple_answer);
        const double frac_diff = diff / std::max(val, simple_answer);
        assert( (frac_diff < 1.0E-6) || (diff < 1.0E-12) );
      }
#endif
      
      lower_energy = right_tail_start;
      break;
    }else
    {
      const T indefinite_high = right_tail_indefinite_non_norm( upper_energy );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );
      check_jet_for_NaN( channels[channel] );

      channels[channel] += val;
      
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
      if constexpr ( std::is_same_v<T, double> )
      {
        const double simple_answer = peak_amplitude * PeakDists::exp_gauss_exp_integral( peak_mean, peak_sigma, skew_left, skew_right, lower_energy, upper_energy );
        const double diff = fabs(val - simple_answer);
        const double frac_diff = diff / std::max(val, simple_answer);
        assert( (frac_diff < 1.0E-6) || (diff < 1.0E-9) );
      }
#endif
      
      lower_energy = upper_energy;
      indefinite_low = indefinite_high;
      channel += 1;
    }//if( upper_energy > stop_energy ) / else
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
}//exp_gauss_exp_integral(...)
  
  
template<typename T>
T DSCB_norm( const T alpha_low, const T n_low, const T alpha_high, const T n_high )
{
  using namespace std;
  assert( alpha_low > 0.0 );
  assert( alpha_high > 0.0 );
  assert( n_low > 1.0 );
  assert( n_high > 1.0 );
  
  const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.70710678118654752440
  const double root_pi = boost::math::constants::root_pi<double>();
  
  
  T a = alpha_low;
  T n = n_low;
  // -(e^(-a*a/2)*n*((a*(-t+n/a-a))/n)^(1-n))/(a*(1-n))
  // -(e^(-a^2/2)*n^n*(-t+n/a-a)^(1-n))/((1-n)*a^n)
  //const double ltail = -(std::exp(-a*a/2)*n*std::pow((a*(a +n/a -a))/n,1-n))/(a*(1-n)); //From Integrating with Maxima
  const T ltail = -exp(-0.5*a*a) * n / (a * (1.0 - n)); //Simplifying, and making more numerically stable, using https://herbie.uwplse.org/demo/
  
  // L = alpha_low
  // R = alpha_high
  // (sqrt(pi)*(erf(R/sqrt(2))-erf(L/sqrt(2))))/sqrt(2)
  T mid;
  if constexpr ( !std::is_same_v<T, double> )
    mid = (root_pi*one_div_root_two*(erf(alpha_high*one_div_root_two)
                                                - erf(-alpha_low*one_div_root_two)));
  else
    mid = (root_pi*one_div_root_two*(boost_erf_imp(alpha_high*one_div_root_two)
                                                - boost_erf_imp(-alpha_low*one_div_root_two)));
  
  a = alpha_high;
  n = n_high;
  // Integrate [e^(-a^2/2)* (((a/n)*((n/a)-a+t))^(-n)] dt, from a to inifit
  // (e^(-a^2/2)*n*((a*(t+n/a-a))/n)^(1-n))/(a*(1-n))
  // (e^(-a^2/2)*n^n*(x+n/a-a)^(1-n))/((1-n)*a^n)
  //const double rtail = -(std::exp(-a*a/2)*n*std::pow((a*(a + n/a -a))/n,1-n))/(a*(1-n)); //From Integrating with Maxima
  const T rtail = -exp(-0.5*a*a) * n / (a * (1.0 - n));  //Simplifying, and making more numerically stable, using https://herbie.uwplse.org/demo/
  
  return 1.0 / (ltail + mid + rtail);
}//DSCB_norm( ... )
  
  
template<typename T>
void double_sided_crystal_ball_integral( const T peak_mean,
                                        const T peak_sigma,
                                        const T peak_amplitude,
                                        const T lower_alpha,
                                        const T lower_power_law,
                                        const T upper_alpha,
                                        const T upper_power_law,
                                        const float * const energies,
                                        T *channels,
                                          const size_t nchannel )
{
  using namespace std;
  
  if( peak_sigma == 0.0 )
    return;

  if constexpr ( std::is_same_v<T, double> )
  {
    // We dont want to return for zero-amplitude peaks if we are usign a Ceres::Jet,
    //  as we may need to take into account the derivatives at zero
    if( peak_amplitude == 0.0 )
      return;
  }//if constexpr ( std::is_same_v<T, double> )

  check_jet_for_NaN( peak_mean );
  check_jet_for_NaN( peak_sigma );
  check_jet_for_NaN( peak_amplitude );
  check_jet_for_NaN( lower_alpha );
  check_jet_for_NaN( lower_power_law );
  check_jet_for_NaN( upper_alpha );
  check_jet_for_NaN( upper_power_law );
  check_jet_array_for_NaN( channels, nchannel );


  // TODO: estimate where we should actually start and stop computing values for, using `double_sided_crystal_ball_coverage_limits(...)`, but need to check if it actually saves time
  //const double zero_amp_point_nsigma = 8.0;
  const T start_energy( static_cast<double>(energies[0]) ); // peak_mean - zero_amp_point_nsigma*peak_sigma );
  const T stop_energy( static_cast<double>(energies[nchannel]) ); // peak_mean + zero_amp_point_nsigma*peak_sigma );

  check_jet_for_NaN( start_energy );
  check_jet_for_NaN( stop_energy );

  size_t channel = 0;
  while( (channel < nchannel) && (energies[channel+1] < start_energy) )
  {
    channel += 1;
  }
  
  if( channel == nchannel )
    return;
  
  const T exp_lower_aa = exp(-0.5*lower_alpha*lower_alpha);
  const T exp_upper_aa = exp(-0.5*upper_alpha*upper_alpha);

  check_jet_for_NaN( exp_lower_aa );
  check_jet_for_NaN( exp_upper_aa );

  auto left_tail_indefinite_non_norm = [peak_mean,peak_sigma,lower_alpha,lower_power_law, exp_lower_aa]( const T x ) -> T {
    const T t = (x - peak_mean) / peak_sigma;
    assert( (t - 1.0E-7) <= -lower_alpha );

    const T &a = lower_alpha;
    const T &n = lower_power_law;
    const T t_1 = 1.0 - ((a + t)*a / n);

    T answer = -exp_lower_aa*(t_1 / pow(t_1, n)) / ((a / n) - a); //slightly more stable

    check_jet_for_NaN( t );
    check_jet_for_NaN( t_1 );
    check_jet_for_NaN( answer );

    return answer;
  };
  
  auto right_tail_indefinite_non_norm = [peak_mean,peak_sigma,upper_alpha,upper_power_law,exp_upper_aa]( const T x ) -> T {
    const T t = (x - peak_mean) / peak_sigma;
    assert( (t + 1.0E-7) >= upper_alpha );
    
    const T &a = upper_alpha;
    const T &n = upper_power_law;
    
    return exp_upper_aa*(1.0 / ((a / n) - a)) * pow((1.0 + ((a * (t - a)) / n)), (1.0 - n));
  };
  
  
  auto gauss_indefinite_non_norm = [peak_mean,peak_sigma]( const T x ) -> T {
    const T t = (x - peak_mean) / peak_sigma;
    const double root_half_pi = boost::math::constants::root_half_pi<double>();
    const double one_div_root_two = boost::math::constants::one_div_root_two<double>(); //0.70710678118654752440
    
    if constexpr ( !std::is_same_v<T, double> )
      return root_half_pi * erf( one_div_root_two * t );
    else
      return root_half_pi * boost_erf_imp( one_div_root_two * t );
  };
  
  
  const T norm = peak_amplitude * DSCB_norm( lower_alpha, lower_power_law, upper_alpha, upper_power_law);
  
  const T left_tail_end = peak_mean - peak_sigma*lower_alpha;
  const T right_tail_start = peak_mean + peak_sigma*upper_alpha;

  check_jet_for_NaN( left_tail_end );
  check_jet_for_NaN( right_tail_start );
  check_jet_for_NaN( norm );

  T lower_energy;
  
  if( energies[channel] < left_tail_end )
  {
    lower_energy = T( static_cast<double>(energies[channel]) );
    T indefinite_low = left_tail_indefinite_non_norm( lower_energy );

    check_jet_for_NaN( lower_energy );
    check_jet_for_NaN( indefinite_low );

    while( (channel < nchannel) && (energies[channel] < left_tail_end) )
    {
      const T upper_energy( static_cast<T>(energies[channel+1]) );

      check_jet_for_NaN( upper_energy );

      if( upper_energy > left_tail_end )
      {
        const T indefinite_high = left_tail_indefinite_non_norm( left_tail_end );
        const T val = norm * (indefinite_high - indefinite_low);

        check_jet_for_NaN( indefinite_high );
        check_jet_for_NaN( val );

        channels[channel] += val;
        
#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double alt_val = peak_amplitude * PeakDists::double_sided_crystal_ball_integral( peak_mean, peak_sigma,
                                                                                              lower_alpha, lower_power_law,
                                                                                              upper_alpha, upper_power_law,
                                                                                              lower_energy, left_tail_end );
          
          const double diff = fabs(val - alt_val);
          const double frac_diff = diff / std::max(val, alt_val);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-5) );
        }//if constexpr ( std::is_same_v<T, double> )
#endif
        lower_energy = left_tail_end;
        
        break;
      }else
      {
        const T indefinite_high = left_tail_indefinite_non_norm( upper_energy );
        const T val = norm * (indefinite_high - indefinite_low);

        check_jet_for_NaN( indefinite_high );
        check_jet_for_NaN( val );

        channels[channel] += val;

        check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double alt_val = peak_amplitude * PeakDists::double_sided_crystal_ball_integral( peak_mean, peak_sigma,
                                                                                              lower_alpha, lower_power_law,
                                                                                              upper_alpha, upper_power_law,
                                                                                              lower_energy, upper_energy );
          
          const double diff = fabs(val - alt_val);
          const double frac_diff = diff / std::max(val, alt_val);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-9) );
        }//if constexpr ( std::is_same_v<T, double> )
#endif
        
        lower_energy = upper_energy;
        indefinite_low = indefinite_high;
        channel += 1;
      }
    }//while( (channel < nchannel) && (energies[channel] < tail_end) )
  }//if( energies[channel] < tail_end )
  
  if( channel >= nchannel )
    return;
  
  assert( static_cast<double>(energies[channel+1]) >= left_tail_end );
  
  lower_energy = max( T(static_cast<double>(energies[channel])), left_tail_end );
  T indefinite_low = gauss_indefinite_non_norm( lower_energy );

  check_jet_for_NaN( lower_energy );
  check_jet_for_NaN( indefinite_low );

  while( (channel < nchannel) && (static_cast<double>(energies[channel]) < right_tail_start) )
  {
    const T upper_energy( static_cast<double>(energies[channel+1]) );

    check_jet_for_NaN( upper_energy );

    if( upper_energy > right_tail_start )
    {
      const T indefinite_high = gauss_indefinite_non_norm( right_tail_start );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );

      channels[channel] += val;

      check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double alt_val = peak_amplitude * PeakDists::double_sided_crystal_ball_integral( peak_mean, peak_sigma,
                                                                                              lower_alpha, lower_power_law,
                                                                                              upper_alpha, upper_power_law,
                                                                                              lower_energy, right_tail_start );
          
          const double diff = fabs(val - alt_val);
          const double frac_diff = diff / std::max(val, alt_val);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-5) );
        }//if constexpr ( std::is_same_v<T, double> )
#endif
      
      lower_energy = right_tail_start;
      break;
    }else
    {
      const T indefinite_high = gauss_indefinite_non_norm( upper_energy );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );

      channels[channel] += val;

      check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double alt_val = peak_amplitude * PeakDists::double_sided_crystal_ball_integral( peak_mean, peak_sigma,
                                                                                              lower_alpha, lower_power_law,
                                                                                              upper_alpha, upper_power_law,
                                                                                              lower_energy, upper_energy );
          
          const double diff = fabs(val - alt_val);
          const double frac_diff = diff / std::max(val, alt_val);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-9) );
        }//if constexpr ( std::is_same_v<T, double> )
#endif
      
      lower_energy = upper_energy;
      indefinite_low = indefinite_high;
      channel += 1;
    }
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )
  
  
  if( channel >= nchannel )
    return;
  
  assert( energies[channel+1] >= right_tail_start );
  lower_energy = max( T(static_cast<double>(energies[channel])), right_tail_start );
  indefinite_low = right_tail_indefinite_non_norm( lower_energy );

  check_jet_for_NaN( indefinite_low );
  check_jet_for_NaN( lower_energy );

  while( (channel < nchannel) && (energies[channel] < stop_energy) )
  {
    const T upper_energy( static_cast<double>(energies[channel+1]) );

    check_jet_for_NaN( upper_energy );

    if( upper_energy > stop_energy )
    {
      const T indefinite_high = right_tail_indefinite_non_norm( right_tail_start );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );

      channels[channel] += val;

      check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double alt_val = peak_amplitude * PeakDists::double_sided_crystal_ball_integral( peak_mean, peak_sigma,
                                                                                              lower_alpha, lower_power_law,
                                                                                              upper_alpha, upper_power_law,
                                                                                              lower_energy, right_tail_start );
          
          const double diff = fabs(val - alt_val);
          const double frac_diff = diff / std::max(val, alt_val);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-9) );
        }//if constexpr ( std::is_same_v<T, double> )
#endif
      
      lower_energy = right_tail_start;
      break;
    }else
    {
      const T indefinite_high = right_tail_indefinite_non_norm( upper_energy );
      const T val = norm * (indefinite_high - indefinite_low);

      check_jet_for_NaN( indefinite_high );
      check_jet_for_NaN( val );

      channels[channel] += val;

      check_jet_for_NaN( channels[channel] );

#if( PERFORM_DEVELOPER_CHECKS && !defined(NDEBUG) )
        if constexpr ( std::is_same_v<T, double> )
        {
          const double alt_val = peak_amplitude * PeakDists::double_sided_crystal_ball_integral( peak_mean, peak_sigma,
                                                                                              lower_alpha, lower_power_law,
                                                                                              upper_alpha, upper_power_law,
                                                                                              lower_energy, upper_energy );
          
          const double diff = fabs(val - alt_val);
          const double frac_diff = diff / std::max(val, alt_val);
          assert( (frac_diff < 1.0E-6) || (diff < 1.0E-9) );
        }//if constexpr ( std::is_same_v<T, double> )
#endif
      
      lower_energy = upper_energy;
      indefinite_low = indefinite_high;
      channel += 1;
    }
  }//while( (channel < nchannel) && (energies[channel] < stop_energy) )

  check_jet_array_for_NaN( channels, nchannel );
}//double_sided_crystal_ball_integral(...)


template<typename T>
void photopeak_function_integral( const T mean,
                                  const T sigma,
                                  const T amp,
                                  const PeakDef::SkewType skew_type,
                                  const T * const skew_parameters,
                                  const size_t nchannel,
                                  const float * const energies,
                                  T *channels )
{
  assert( (skew_type == PeakDef::SkewType::NoSkew) || skew_parameters );

  switch( skew_type )
  {
    case PeakDef::SkewType::NoSkew:
      gaussian_integral( mean, sigma, amp, energies, channels, nchannel );
      break;

    case PeakDef::SkewType::Bortel:
      bortel_integral( mean, sigma, amp, skew_parameters[0], energies, channels, nchannel );
      break;

    case PeakDef::SkewType::CrystalBall:
      crystal_ball_integral( mean, sigma, amp, skew_parameters[0], skew_parameters[1], energies, channels, nchannel );
      break;

    case PeakDef::SkewType::DoubleSidedCrystalBall:
      double_sided_crystal_ball_integral( mean, sigma, amp,
                                         skew_parameters[0], skew_parameters[1],
                                         skew_parameters[2], skew_parameters[3],
                                         energies, channels, nchannel );
      break;

    case PeakDef::SkewType::GaussExp:
      gauss_exp_integral( mean, sigma, amp, skew_parameters[0], energies, channels, nchannel );
      break;

    case PeakDef::SkewType::ExpGaussExp:
      exp_gauss_exp_integral( mean, sigma, amp, skew_parameters[0], skew_parameters[1], energies, channels, nchannel );
      break;
  }//switch( skew_type )
}//void photopeak_function_integral(...)


#if( __cplusplus >= 202002L )
void offset_integral( const ContType &cont,
                     const float *energies,
                     T *channels,
                     const size_t nchannel,
                     const std::shared_ptr<const SpecUtils::Measurement> &data ) requires ContinuumTypeConcept<ContType,T>
#else
template <typename ContType, typename T>
typename std::enable_if<ContinuumTypeConcept<ContType, T>::value, void>::type
offset_integral(const ContType& cont,
                const float* energies,
                T* channels,
                const size_t nchannel,
                const std::shared_ptr<const SpecUtils::Measurement>& data)
#endif
{
  // This function should give the same answer as
  //  `PeakContinuum::offset_integral( double x0, const double x1, data)`, on a channel-by-channel
  //  basis, just be a little faster computationally, especially for stepped continua.

  using namespace std;

  assert( nchannel > 0 );
  if( !nchannel )
    return;

  const T reference_energy = cont.referenceEnergy();

  const auto &pars = cont.parameters();

  const PeakContinuum::OffsetType cont_type = cont.type();

  switch( cont_type )
  {
    case PeakContinuum::OffsetType::NoOffset:
      return;

    case PeakContinuum::OffsetType::Constant:
    case PeakContinuum::OffsetType::Linear:
    case PeakContinuum::OffsetType::Quadratic:
    case PeakContinuum::OffsetType::Cubic:
    {
      for( size_t i = 0; i < nchannel; ++i )
      {
        const T x0 = static_cast<double>(energies[i]) - reference_energy;
        const T x1 = static_cast<double>(energies[i+1]) - reference_energy;

        T answer(0.0);
        switch( cont_type )
        {
          case PeakContinuum::OffsetType::NoOffset:
          case PeakContinuum::OffsetType::External:
          case PeakContinuum::OffsetType::FlatStep:
          case PeakContinuum::OffsetType::LinearStep:
          case PeakContinuum::OffsetType::BiLinearStep:
            assert( 0 );

          case PeakContinuum::OffsetType::Cubic:
            assert( pars.size() >= 4 );
            answer += 0.25*pars[3]*(x1*x1*x1*x1 - x0*x0*x0*x0);
            //fall-through intentional

          case PeakContinuum::OffsetType::Quadratic:
            assert( pars.size() >= 3 );
            answer += 0.333333333333333*pars[2]*(x1*x1*x1 - x0*x0*x0);
            //fall-through intentional

          case PeakContinuum::OffsetType::Linear:
            assert( pars.size() >= 2 );
            answer += 0.5*pars[1]*(x1*x1 - x0*x0);
            //fall-through intentional

          case PeakContinuum::OffsetType::Constant:
            assert( pars.size() >= 1 );
            answer += pars[0]*(x1 - x0);
            break;
        };//switch( type )

        if constexpr ( std::is_same_v<T, double> )
        {
          assert( std::max(answer,0.0) == cont.offset_eqn_integral( &(pars[0]), cont_type, energies[i], energies[i+1], reference_energy ) );
        }

        channels[i] += max( answer, T(0.0) );
      }//for( size_t i = 0; i < nchannel; ++i )

      break;
    }//case Constant: case Linear: case Quadratic: case Cubic:


    case PeakContinuum::OffsetType::FlatStep:
    case PeakContinuum::OffsetType::LinearStep:
    case PeakContinuum::OffsetType::BiLinearStep:
    {
      if( !data || !data->num_gamma_channels() )
        throw std::runtime_error( "PeakContinuum::offset_integral: invalid data spectrum passed in" );

      // To be consistent with how fit_amp_and_offset(...) handles things, we will do our own
      //  summing here, rather than calling Measurement::gamma_integral.
      double lowerEnergy, upperEnergy;

      if constexpr ( std::is_same_v<T, double> )
      {
        lowerEnergy = cont.lowerEnergy();
        upperEnergy = cont.upperEnergy();
      }else
      {
        lowerEnergy = cont.lowerEnergy().a;
        upperEnergy = cont.upperEnergy().a;
      }

      const size_t roi_lower_channel = data->find_gamma_channel( lowerEnergy );
      const size_t roi_upper_channel = data->find_gamma_channel( upperEnergy );

      //const double roi_lower = data->gamma_channel_lower(roi_lower_channel);
      //const double roi_upper = data->gamma_channel_upper(roi_upper_channel);

      const std::vector<float> &counts = *data->gamma_counts();

      assert( roi_lower_channel < counts.size() );
      assert( roi_upper_channel < counts.size() );


      const double roi_data_sum = std::accumulate( begin(counts) + roi_lower_channel, begin(counts) + roi_upper_channel + 1,  0.0 );
      // Might be able to take advantage of vectorized sum using something like:
      //const double roi_data_sum = Eigen::VectorXf::Map( &(counts[roi_lower_channel]), (1 + roi_upper_channel - roi_lower_channel) ).sum();


      const size_t begin_channel = data->find_gamma_channel( energies[0] );
      assert( energies[0] == data->gamma_channel_lower(begin_channel) );
      const size_t end_channel = begin_channel + nchannel; //one past last channel we want

      // Lets check that `energies` points into the lower channel energies of `data`.
      //  We actually only care that the values of the array are the same, but we'll be
      //  a little tighter for development.
      const std::vector<float> &data_energies = *data->channel_energies();

      if( (energies[0] != data->gamma_channel_lower(begin_channel))
         || (energies[nchannel] != data->gamma_channel_lower(begin_channel+nchannel)) )
        throw std::logic_error( "PeakContinuum::offset_integral: for stepped continua" );

      double cumulative_data = 0.0;

      // Incase we are starting
      if( begin_channel > roi_lower_channel )
      {
        for( size_t i = begin_channel; i < begin_channel; ++i )
          cumulative_data += counts[i];
      }//if( begin_channel > lower_channel )


      for( size_t i = begin_channel; i < end_channel; ++i )
      {
        const size_t input_index = i - begin_channel;
        assert( data_energies[i] == energies[input_index] );

        const T x0_rel = static_cast<double>(data_energies[i]) - reference_energy;
        const T x1_rel = static_cast<double>(data_energies[i+1]) - reference_energy;

        if( i >= roi_lower_channel && i <= roi_upper_channel )
          cumulative_data += 0.5*counts[i];

        const double frac_data = cumulative_data / roi_data_sum;

        switch( cont_type )
        {
          case PeakContinuum::OffsetType::FlatStep:
          case PeakContinuum::OffsetType::LinearStep:
          {
            const T offset = pars[0]*(x1_rel - x0_rel);
            const T linear = ((cont_type == PeakContinuum::OffsetType::FlatStep) ? T(0.0) :  0.5*pars[1]*(x1_rel*x1_rel - x0_rel*x0_rel));
            const size_t step_index = ((cont_type == PeakContinuum::OffsetType::FlatStep) ? 1 : 2);
            const T step_contribution = pars[step_index] * frac_data * (x1_rel - x0_rel);

            const T answer = max( T(0.0), offset + linear + step_contribution );

            if constexpr ( std::is_same_v<T, double> )
            {
              assert( answer == cont.offset_integral( data_energies[i], data_energies[i+1], data ) );
            }

            channels[input_index] += answer;
            break;
          }//case FlatStep: case LinearStep:

          case PeakContinuum::OffsetType::BiLinearStep:
          {
            assert( pars.size() == 4 );
            const T left_poly = pars[0]*(x1_rel - x0_rel) + 0.5*pars[1]*(x1_rel*x1_rel - x0_rel*x0_rel);
            const T right_poly = pars[2]*(x1_rel - x0_rel) + 0.5*pars[3]*(x1_rel*x1_rel - x0_rel*x0_rel);
            const T contrib = std::max( T(0.0), ((1.0 - frac_data) * left_poly) + (frac_data * right_poly) );

            if constexpr ( std::is_same_v<T, double> )
            {
              assert( contrib == cont.offset_integral( data_energies[i], data_energies[i+1], data ) );
            }

            channels[input_index] += contrib;
            break;
          }//case BiLinearStep:

          case PeakContinuum::OffsetType::NoOffset:
          case PeakContinuum::OffsetType::Constant:
          case PeakContinuum::OffsetType::Linear:
          case PeakContinuum::OffsetType::Quadratic:
          case PeakContinuum::OffsetType::Cubic:
          case PeakContinuum::OffsetType::External:
            assert( 0 );
            break;
        }//switch( cont_type )

        if( i >= roi_lower_channel && i <= roi_upper_channel )
          cumulative_data += 0.5*counts[i];
      }//for( size_t i = 0; i < channels; ++i )

      break;
    }//case FlatStep/LinearStep/BiLinearStep

    case PeakContinuum::OffsetType::External:
    {
      std::shared_ptr<const SpecUtils::Measurement> ext_cont = cont.externalContinuum();
      assert( ext_cont );
      if( !ext_cont )
        break;

      for( size_t i = 0; i < nchannel; ++i )
        channels[i] += ext_cont ? ext_cont->gamma_integral( energies[i], energies[i+1] ) : 0.0;
      break;
    }//case PeakContinuum::OffsetType::External:
  }//switch( cont_type )
}//void PeakContinuum::offset_integral( ... ) const
}//namespace PeakDists

#endif //PeakDists_imp_h
