#ifndef PeakDists_imp_h
#define PeakDists_imp_h

#include <boost/math/constants/constants.hpp>


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
  if( peak_sigma==0.0 || peak_amplitude==0.0 )
    return;
    
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
    
  // We will keep track of the channels lower value of erf, so we dont have to re-compute
  //  it for each channel (this is the who advantage of )
  T erflow;
  if constexpr ( !std::is_same_v<T, double> )
    erflow = erf( (static_cast<double>(energies[channel]) - peak_mean)/(sqrt2*peak_sigma) );
  else
    erflow = boost_erf_imp( (energies[channel] - peak_mean)/(sqrt2*peak_sigma) );
    
  while( (channel < nchannel) && (energies[channel] < stop_energy) )
  {
    const T erfhigharg = (static_cast<double>(energies[channel+1]) - peak_mean)/(sqrt2*peak_sigma);
    T erfhigh;
    if constexpr ( !std::is_same_v<T, double> )
      erfhigh = erf( erfhigharg );
    else
      erfhigh = boost_erf_imp( erfhigharg );
    
    channels[channel] += 0.5 * peak_amplitude * (erfhigh - erflow);
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
    
  if( (peak_sigma == 0.0) || (peak_amplitude == 0.0) )
    return;
    
    
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
      assert( (frac_diff < 1.0E-5) || (diff < 1.0E-8) );
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
  
  if( (peak_sigma == 0.0) || (peak_amplitude == 0.0) )
    return;
  
  
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
  
  if( (peak_sigma == 0.0) || (peak_amplitude == 0.0) )
    return;
  
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
  
  if( (peak_sigma == 0.0) || (peak_amplitude == 0.0) )
    return;

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
  
}//namespace PeakDists

#endif //PeakDists_imp_h
