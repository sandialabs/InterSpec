#ifndef PeakDists_h
#define PeakDists_h
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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <memory>
#include <vector>
#include <utility>

#include "InterSpec/PeakDef.h" //for PeakDef::SkewType

namespace SpecUtils
{
  class Measurement;
}

/** The functions in this .h/.cpp are for computing skewed and Gaussian photopeak distributions.
 
 */
namespace PeakDists
{
  /** `erfc` implementation based on boost, but sped up to be more efficient.
   
   Surprisingly, the erf() function is the major bottleneck for peak fitting.
   
   20191230: wcjohns extracted the boost::math::erf() function implementation
   from boost 1.65.1 for double precision (53 bit mantissa) into this function,
   boost_erf_imp(). Removing some of the supporting code structure, and
   explicitly writing out the polynomial equation evaluation, and removing some
   very minor corrections, seems to speed things up by about a factor of ~3 over
   calling boost::math::erf() (and agrees with boosts implementation to within
   1E-10%, across the entire range).
   
   In the implementation file, there is a commented out erf_approx() function looks to be about 25% faster than
   this boost version, but I havent carefully checked out the precision implications
   so not switching to it yet.
   */
  double boost_erf_imp( double z );
  
  /** Returns `1-boost_erf_imp(z)`, however, due to rounding, is a separate implementation than erf.
   */
  double boost_erfc_imp( double z );

  /** Function to semi-efficiently integrate the gaussian plus optionally skew distribution over an energy range.
   
   @param mean The peak mean, in keV
   @param sigma The peak sigma, in keV
   @param amplitude The peak amplitude; e.g., use 1.0 for unit-area peak.
   @param skew_type The type of skew of the peak
   @param skew_parameters The values of the skew parameters; must have at least `num_skew_parameters(SkewType)` entries.
          or can be nullptr if no skew.
   @param nchannel The number of channels to sum
   @param lower_energies The lower energies of the channels.  Must have at least `nchannel + 1` entries
   @param[out] peak_count_channels Where the channel sums of the distribution are added to.
          Note that the distribution sum for each channel is _added_ to this array, so you should zero-initialize it.
          Must have at least `nchannel` entries.
   */
template<typename T>
void photopeak_function_integral( const T mean,
                                  const T sigma,
                                  const T amp,
                                  const PeakDef::SkewType skew_type,
                                  const T * const skew_parameters,
                                  const size_t nchannel,
                                  const float * const energies,
                                 T *channels );

extern template void photopeak_function_integral<double>( const double, const double,const double,
                         const PeakDef::SkewType, const double * const, const size_t, const float * const, double * );


  
  
  
  
  
  /** Calculates the area of a Gaussian with specified mean, sigma, and amplitude, between x0 and x1.
   
    Results have approximately 9 decimal digits of accuracy.
   */
  double gaussian_integral( const double peak_mean, const double peak_sigma,
                            const double x0, const double x1 );
  
  /** Slightly CPU optimized method of computing the Gaussian integral over a number of channels.
   
   Cuts the number of calls to the `erf` function (which is what takes the longest in the
   function) in half.
   Also, only calculates values between +-8 sigma of the mean (which is 1 - 1E-15 the total counts).
   
   @param peak_mean
   @param peak_sigma
   @param peak_amplitude
   @param energies Array of lower channel energies; must have at least one more entry than
          `nchannel`
   @param channels Channel count array integrals of Gaussian and Skew will be _added_ to (e.g.,
          will not be zeroed); must have at least `nchannel` entries
   @param nchannel The number of channels to do the integration over.
   */
  template<typename T>
  void gaussian_integral( const T peak_mean,
                                const T peak_sigma,
                                const T peak_amplitude,
                                const float * const energies,
                                T *channels,
                         const size_t nchannel );
  
  extern template void gaussian_integral<double>(const double, const double, const double,
                                  const float * const, double *, const size_t);
  
  //void gaussian_integral( const double peak_mean,
  //                            const double peak_sigma,
  //                            const double peak_amplitude,
  //                            const float * const energies,
  //                            double *channels,
  //                            const size_t nchannel );
  
  
  
  /** Returns the integral of a unit-area Bortel function, between `x1` and `x2`. */
  double bortel_integral( const double mean, const double sigma, const double skew,
                         const double x1, const double x2 );
  
  /** Slightly CPU optimized method of computing the peak area, for Bortel skew over a number of channels.
   
   Cuts the number of calls to the `erf`, `erfc`, and exp functionsn half.
   Also, only calculates values between -12 to +8 sigma of the mean/
   
   @param peak_mean
   @param peak_sigma
   @param peak_amplitude
   @param skew The Bortel skew of the peak (between 0 and 10)
   @param energies Array of lower channel energies; must have at least one more entry than
          `nchannel`
   @param channels Channel count array integrals of Gaussian and Skew will be _added_ to (e.g.,
          will not be zeroed); must have at least `nchannel` entries
   @param nchannel The number of channels to do the integration over.
   */
  template<typename T>
  void bortel_integral( const T peak_mean,
                              const T peak_sigma,
                              const T peak_amplitude,
                              const T skew,
                              const float * const energies,
                              T *channels,
                              const size_t nchannel );
  extern template void bortel_integral<double>( const double, const double, const double, const double,
                                       const float * const, double *, const size_t );
  
  /** Returns the PDF for a unit-area Bortel function.
   */
  double bortel_pdf( const double mean, const double sigma, const double skew_low, const double x );

  /** Returns the indefinite integral (e.g., from negative infinity up to `x`) of a unit-area Bortel function
   */
  template<typename T>
  T bortel_indefinite_integral( const double x, const T mean, const T sigma, const T skew );
  
  extern template double bortel_indefinite_integral<double>( const double, const double,
                                                            const double, const double );

  /** Returns an approximate area between `x1` and `x2` for a unit-area Bortel function.
   
   Just multiplies the x-range by the Bortel PDF value in the middle of the range.
   */
  double bortel_integral_fast( const double mean, const double sigma, const double skew,
                              const double x1, const double x2 );
  
  /** Return the limits so that `1-p` of the Bortel distribution is covered.
   
   @param mean The peak mean
   @param sigma The peak width
   @param skew The Bortel skew of the peak (between 0 and 10)
   @param The fraction of the distribution you want to be outside of the returned range; must
          be between 1.0E-11 and 0.999.
   
   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.
   
   The area above and below the retuned limits is accurate to within 0.02 of requested area (i.e. `0.02*p`).
   
   Throws error on invalid input.
   */
  std::pair<double,double> bortel_coverage_limits( const double mean, const double sigma,
                                         const double skew, const double p );
  
  
  
  double gauss_exp_integral( const double peak_mean,
                                      const double peak_sigma,
                                      const double skew,
                                      const double x0, const double x1 );
  
  
  
  template<typename T>
  void gaussian_integral( const T peak_mean,
                                const T peak_sigma,
                                const T peak_amplitude,
                                const float * const energies,
                                T *channels,
                         const size_t nchannel );
  
  extern template void gaussian_integral<double>(const double, const double, const double,
                                  const float * const, double *, const size_t);
  
  template<typename T>
  void gauss_exp_integral( const T peak_mean,
                          const T peak_sigma,
                          const T peak_amplitude,
                          const T skew,
                          const float * const energies,
                          T *channels,
                          const size_t nchannel );
  
  extern template void gauss_exp_integral<double>( const double, const double, const double,
                                      const double, const float * const, double *, const size_t );
  
  
  /** Returns the normalization so the GaussExp distribution has unit area. */
  template<typename T>
  T gauss_exp_norm( const T sigma, const T skew );
  
  extern template double gauss_exp_norm<double>( const double, const double );
  
  
  double gauss_exp_pdf(const double mean,
                       const double sigma,
                       const double skew,
                       const double x );
  
  double gauss_exp_tail_indefinite(const double mean,
                                  const double sigma,
                                  const double skew,
                                   const double x );
  
  double gauss_exp_indefinite(const double mean,
                             const double sigma,
                             const double skew,
                              const double x );
  
  /** Return the limits so that `1-p` of the GaussExp distribution is covered.
   
   @param mean The peak mean
   @param sigma The peak width
   @param skew The GaussExp skew of the peak (between 0.15 and 3.25)
   @param The fraction of the distribution you want to be outside of the returned range; must
          be between 1.0E-11and 0.999.
   
   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.
   
   The area above and below the retuned limits is computed by inverting the indefinite integral, so should be
   decently accurate up to numeric precision (which hasn't been checked/optimized).
   
   Throws error on invalid input.
   */
  std::pair<double,double> gauss_exp_coverage_limits( const double mean, const double sigma,
                                         const double skew, const double p );
  
  
  
  
  
  
  
  double exp_gauss_exp_integral( const double mean,
                           const double sigma,
                           const double skew_left,
                           const double skew_right,
                           const double x0,
                                const double x1 );
  
  template<typename T>
  void exp_gauss_exp_integral( const T peak_mean,
                               const T peak_sigma,
                               const T peak_amplitude,
                               const T skew_left,
                               const T skew_right,
                               const float * const energies,
                              T *channels,
                               const size_t nchannel );
  
  extern template void exp_gauss_exp_integral<double>( const double, const double, const double,
                                                  const double, const double, const float * const,
                                                  double *, const size_t );
  
  template<typename T>
  T exp_gauss_exp_norm( const T sigma, const T skew_left, const T skew_right );
  
  extern template double exp_gauss_exp_norm<double>( const double, const double, const double );
  
  double exp_gauss_exp_pdf( const double mean,
                           const double sigma,
                           const double skew_left,
                           const double skew_right,
                           const double x );
  
  double exp_gauss_exp_left_tail_indefinite( const double mean,
                                          const double sigma,
                                          const double skew_left,
                                          const double skew_right,
                                            const double x );

  double exp_gauss_exp_right_tail_indefinite( const double mean,
                                           const double sigma,
                                           const double skew_left,
                                           const double skew_right,
                                             const double x );

  double exp_gauss_exp_gauss_indefinite( const double mean,
                                      const double sigma,
                                      const double skew_left,
                                      const double skew_right,
                                        const double x );
  
  double exp_gauss_exp_indefinite( const double mean,
                                const double sigma,
                                const double skew_left,
                                const double skew_right,
                                  const double x );
  
  /** Return the limits so that `1-p` of the ExpGaussExp distribution is covered.
   
   @param mean The peak mean
   @param sigma The peak width
   @param skew_left The left-sided ExpGaussExp skew of the peak (between 0.15 and 3.25)
   @param skew_right The right-sided ExpGaussExp skew of the peak (between 0.15 and 3.25)
   @param The fraction of the distribution you want to be outside of the returned range; must
          be between 1.0E-11and 0.999.
   
   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.
   
   The area above and below the retuned limits is computed by inverting the indefinite integral, so should be
   decently accurate up to numeric precision (which hasn't been checked/optimized).
   
   Throws error on invalid input.
   */
  std::pair<double,double> exp_gauss_exp_coverage_limits( const double mean, const double sigma,
                                const double left_skew, const double right_skew, const double p );
  
  
  
  
  
  
  
  double crystal_ball_integral( const double peak_mean,
                                      const double peak_sigma,
                                      const double alpha,
                                      const double power_law,
                                      const double x0, const double x1 );
  
  template<typename T>
  void crystal_ball_integral( const T peak_mean,
                              const T peak_sigma,
                              const T peak_amplitude,
                              const T alpha,
                              const T power_law,
                              const float * const energies,
                              T *channels,
                              const size_t nchannel );
  
  extern template void crystal_ball_integral<double>( const double, const double, const double,
                        const double, const double, const float * const, double *, const size_t );
  
  
  double crystal_ball_norm( const double sigma,
                           const double alpha,
                           const double n );
  
  double crystal_ball_pdf(const double mean,
                          const double sigma,
                          const double alpha,
                          const double n,
                          const double x );

  /** Returns indefinite integral (negative infinite to `x0` for the power-law component for the unit-area Crystal Ball function.
   */
  double crystal_ball_tail_indefinite_t( const double sigma, const double alpha,
                                        const double n, const double t );
  
  /** Return the limits so that `1-p` of the Crustal Ball distribution is covered.
   
   @param mean The peak mean
   @param sigma The peak width
   @param alpha The Crystal Ball skew of the peak
   @param n The Crystal Ball power-law of the peak
   @param The fraction of the distribution you want to be outside of the returned range; must
          be between 1.0E-11and 0.999.
   
   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.
   
   The area above and below the retuned limits is computed by inverting the indefinite integral, so should be
   decently accurate up to numeric precision (which hasn't been checked/optimized).
   
   Throws error on invalid input.
   */
  std::pair<double,double> crystal_ball_coverage_limits( const double mean, const double sigma,
                                                        const double alpha,
                                                        const double n,
                                                        const double p );
  
  
  
  
  
  
  double double_sided_crystal_ball_integral( const double peak_mean,
                                                   const double peak_sigma,
                                                   const double lower_alpha,
                                                   const double lower_power_law,
                                                   const double upper_alpha,
                                                   const double upper_power_law,
                                                   const double x0, const double x1 );
  
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
                                          const size_t nchannel );
  
  extern template void double_sided_crystal_ball_integral<double>( const double, const double,
                                                 const double, const double, const double,
                                                 const double, const double, const float * const,
                                                 double *, const size_t );
    
  /** Return the limits so that `1-p` of the ExpGaussExp distribution is covered.
   
   @param mean The peak mean
   @param sigma The peak width
   @param alpha_left The left-sided Crystal Ball skew of the peak
   @param n_left The left-sided Crystal Ball power-law of the peak
   @param alpha_right The right-sided Crystal Ball skew of the peak
   @param n_right The right-sided Crystal Ball power-law of the peak
   @param The fraction of the distribution you want to be outside of the returned range; must
          be between 1.0E-11and 0.999.
   
   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.
   
   The area below the retuned limit should be quite accurate, while above may only be good to within 0.02 of requested area (i.e. `0.02*p`).
   
   Throws error on invalid input - or if limit cant be found (which can happen for really larger skews).
   */
  std::pair<double,double> double_sided_crystal_ball_coverage_limits( const double mean, const double sigma,
                                                                     const double left_skew,
                                                                     const double left_n,
                                                                     const double right_skew,
                                                                     const double right_n,
                                                                     const double p );



  // ========== Voigt Plus Bortel Distribution Functions ==========

  /** Control whether to use pseudo-Voigt or true-Voigt for VoigtPlusBortel peaks.
   
   The true-Voigt distribution (using Faddeeva function) does not have a CDF function
   implemented due to accuracy issues with erf(z) for complex z. Therefore:
   - When USE_PSEUDO_VOIGT_DISTRIBUTION=0: Uses true-Voigt PDF with Gauss-Legendre
     quadrature for channel integrals (slower but more accurate for Voigt component)
   - When USE_PSEUDO_VOIGT_DISTRIBUTION=1: Uses pseudo-Voigt (Thompson-Cox-Hastings
     approximation) with analytic CDF for fast channel integrals (faster, slightly
     less accurate for Voigt component)
   
   Note: voigt_exp_coverage_limits() always uses pseudo-Voigt due to CDF requirement.
   */
  #ifndef USE_PSEUDO_VOIGT_DISTRIBUTION
  #define USE_PSEUDO_VOIGT_DISTRIBUTION 1
  #endif

  /** Returns the integral of a Voigt with exponential tail distribution between x0 and x1.

   Note: for `USE_PSEUDO_VOIGT_DISTRIBUTION == 0`, this function uses true-Voigt with adaptive 
     Gauss-Kronrod quadrature (slower but more accurate for Voigt component) to get the integral
     (this is because of the aformentioned accuracy issues with the true-Voigt CDF)

   @param peak_mean The peak mean in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau
   @param x0 Lower integration limit
   @param x1 Upper integration limit

   @returns Integral of the unit-area distribution from x0 to x1
   */
  double voigt_exp_integral( const double peak_mean, const double sigma_gauss,
                             const double gamma_lor, const double tail_ratio,
                             const double tail_slope, const double x0, const double x1 );


  /** Optimized array-filling version of the Voigt with exponential tail integral.

   Uses indefinite integral caching to cut the number of Voigt function evaluations in half.

   @param peak_mean The peak mean in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param peak_amplitude The peak amplitude (use 1.0 for unit-area peak)
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau
   @param energies Array of channel lower energies (must have nchannel+1 entries)
   @param channels Array where distribution values will be added (must have nchannel entries)
   @param nchannel Number of channels to integrate over
   */
  template<typename T>
  void voigt_exp_integral( const T peak_mean, const T sigma_gauss,
                           const T peak_amplitude, const T gamma_lor,
                           const T tail_ratio, const T tail_slope,
                           const float * const energies, T *channels,
                           const size_t nchannel );

  extern template void voigt_exp_integral<double>( const double, const double, const double,
                                                    const double, const double, const double,
                                                    const float * const, double *, const size_t );


  /** Return the limits so that `1-p` of the VoigtExpTail distribution is covered.

   @param mean The peak mean
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau
   @param p The fraction of the distribution you want to be outside of the returned range

   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.

   @note This function always uses the pseudo-Voigt approximation (regardless of the
   USE_PSEUDO_VOIGT_DISTRIBUTION setting) because it requires CDF evaluation, and the
   true-Voigt CDF is not reliably implemented. This may introduce small accuracy differences
   compared to the actual peak shape when USE_PSEUDO_VOIGT_DISTRIBUTION=0.

   Throws error on invalid input.
   */
  std::pair<double,double> voigt_exp_coverage_limits( const double mean, const double sigma_gauss,
                                                       const double gamma_lor, const double tail_ratio,
                                                       const double tail_slope, const double p );


  // ========== Gauss Plus Bortel Distribution Functions ==========

  /** Returns the integral of a Gauss+Bortel distribution between x0 and x1.

   This is a weighted mixture of Gaussian and Bortel distributions:
   PDF(x) = (1-R) * Gaussian(x) + R * Bortel(x)

   When R=0, this equals a pure Gaussian.
   When R=1, this equals a pure Bortel.
   This distribution should match VoigtPlusBortel when gamma_lor=0.

   @param peak_mean The peak mean in keV
   @param sigma Gaussian width (detector resolution), in keV
   @param R Mixing ratio (0=pure Gaussian, 1=pure Bortel)
   @param tau Bortel skew parameter (exponential tail decay constant, in units of sigma)
   @param x0 Lower integration limit
   @param x1 Upper integration limit

   @returns Integral of the unit-area distribution from x0 to x1
   */
  double gauss_plus_bortel_integral( const double peak_mean, const double sigma,
                                     const double R, const double tau,
                                     const double x0, const double x1 );


  /** Optimized array-filling version of the Gauss+Bortel integral.

   @param peak_mean The peak mean in keV
   @param sigma Gaussian width (detector resolution), in keV
   @param peak_amplitude The peak amplitude (use 1.0 for unit-area peak)
   @param R Mixing ratio (0=pure Gaussian, 1=pure Bortel)
   @param tau Bortel skew parameter
   @param energies Array of channel lower energies (must have nchannel+1 entries)
   @param channels Array where distribution values will be added (must have nchannel entries)
   @param nchannel Number of channels to integrate over
   */
  template<typename T>
  void gauss_plus_bortel_integral( const T peak_mean, const T sigma,
                                   const T peak_amplitude, const T R, const T tau,
                                   const float * const energies, T *channels,
                                   const size_t nchannel );

  extern template void gauss_plus_bortel_integral<double>( const double, const double, const double,
                                                           const double, const double,
                                                           const float * const, double *, const size_t );


  /** Return the limits so that `1-p` of the Gauss+Bortel distribution is covered.

   @param mean The peak mean
   @param sigma Gaussian width (detector resolution), in keV
   @param R Mixing ratio (0=pure Gaussian, 1=pure Bortel)
   @param tau Bortel skew parameter
   @param p The fraction of the distribution you want to be outside of the returned range

   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.

   Throws error on invalid input.
   */
  std::pair<double,double> gauss_plus_bortel_coverage_limits( const double mean, const double sigma,
                                                               const double R, const double tau,
                                                               const double p );


  // ========== Double Bortel Distribution Functions ==========

  /** Returns the integral of a DoubleBortel distribution between x0 and x1.

   From Bortels & Collaers 1987 (Eq. 11), this is a weighted sum of two Bortel distributions:
   PDF(x) = (1-eta) * Bortel(mean, sigma, tau1, x) + eta * Bortel(mean, sigma, tau2, x)

   Where tau2 = tau1 + tau2_delta, ensuring tau2 >= tau1.

   When tau2_delta=0, this reduces to a single Bortel distribution (regardless of eta).
   When eta=0, this equals Bortel with tau1.
   When eta=1, this equals Bortel with tau2.

   @param peak_mean The peak mean in keV
   @param sigma Gaussian width (detector resolution), in keV
   @param tau1 First exponential decay constant (in units of sigma)
   @param tau2_delta Non-negative delta so tau2 = tau1 + tau2_delta
   @param eta Weight of second exponential (0 to 1)
   @param x0 Lower integration limit
   @param x1 Upper integration limit

   @returns Integral of the unit-area distribution from x0 to x1
   */
  double double_bortel_integral( const double peak_mean, const double sigma,
                                 const double tau1, const double tau2_delta,
                                 const double eta,
                                 const double x0, const double x1 );


  /** Optimized array-filling version of the DoubleBortel integral.

   @param peak_mean The peak mean in keV
   @param sigma Gaussian width (detector resolution), in keV
   @param peak_amplitude The peak amplitude (use 1.0 for unit-area peak)
   @param tau1 First exponential decay constant
   @param tau2_delta Non-negative delta so tau2 = tau1 + tau2_delta
   @param eta Weight of second exponential (0 to 1)
   @param energies Array of channel lower energies (must have nchannel+1 entries)
   @param channels Array where distribution values will be added (must have nchannel entries)
   @param nchannel Number of channels to integrate over
   */
  template<typename T>
  void double_bortel_integral( const T peak_mean, const T sigma,
                               const T peak_amplitude,
                               const T tau1, const T tau2_delta, const T eta,
                               const float * const energies, T *channels,
                               const size_t nchannel );

  extern template void double_bortel_integral<double>( const double, const double, const double,
                                                       const double, const double, const double,
                                                       const float * const, double *, const size_t );


  /** Return the limits so that `1-p` of the DoubleBortel distribution is covered.

   @param mean The peak mean
   @param sigma Gaussian width (detector resolution), in keV
   @param tau1 First exponential decay constant
   @param tau2_delta Non-negative delta so tau2 = tau1 + tau2_delta
   @param eta Weight of second exponential (0 to 1)
   @param p The fraction of the distribution you want to be outside of the returned range

   @returns limits so that `0.5*p` of the distribution will be below the first element, and
   `0.5*p` will be above the second element, so the fraction of the distribution between the
   returned limits is `1 - p`.

   Throws error on invalid input.
   */
  std::pair<double,double> double_bortel_coverage_limits( const double mean, const double sigma,
                                                           const double tau1, const double tau2_delta,
                                                           const double eta, const double p );


  /** Returns [lower, upper] energy limits such that fraction `p` of the peak's area lies outside
   the range (i.e., `p/2` below lower, `p/2` above upper), dispatching to the appropriate
   distribution's coverage-limits function based on skew type.

   For NoSkew, uses the standard Gaussian quantile.

   @param p           Fraction of total area outside returned range; must be in (0, 1).
   @param skew_type   The peak skew type.
   @param mean        Peak mean in keV.
   @param sigma       Peak Gaussian sigma in keV.
   @param skew_pars   Pointer to skew parameter array (SkewPar0, SkewPar1, ...); may be nullptr
                      for NoSkew.
   @returns  Pair [lower_energy, upper_energy].
   Throws on invalid parameters or if limits cannot be found.
   */
  std::pair<double,double> coverage_limits( const double p,
                                            const PeakDef::SkewType skew_type,
                                            const double mean,
                                            const double sigma,
                                            const double *skew_pars );


    /** 20231109: Currently `DSCB_pdf_non_norm` yields a that can be almost 20% too low, over the entire effective range
     I'm totally not sure why - it could be an error in my integration for normalization, a bug in the indefinite integral functions,
     or likely a mis-understanding on my part, or a bug in `DSCB_pdf_non_norm` that I have just become blind to seeing.
     So for the moment, we wont use the equivalent of `DSCB_pdf_non_norm` in the JS, which is faster and more compact,
     see test\_PeakDists.cpp for some commented out tests that show this issue.  I'm a bit miffed at this inconsistency.
     
     Returns the non-normalized, double sided Crystal Ball PDF value for `x`
     */
    //double DSCB_pdf_non_norm( const double mean, const double sigma,
    //                          const double alpha_low, const double n_low,
    //                          const double alpha_high, const double n_high,
    //                          const double x );
  

  template<typename T>
  T DSCB_norm( const T alpha_low, const T n_low, const T alpha_high, const T n_high );
  
  extern template double DSCB_norm<double>( const double, const double, const double, const double );
    
  double DSCB_left_tail_indefinite_non_norm_t( const double alpha_low, const double n_low,
                                                 const double t );
    
  double DSCB_right_tail_indefinite_non_norm_t( const double alpha_high,
                                                             const double n_high,
                                                             const double t);
    
  double DSCB_gauss_indefinite_non_norm_t( const double t );

#if( __cplusplus >= 202002L )
  template <typename ContType, typename ScalarType>
  concept ContinuumTypeConcept = requires(ContType cont, ScalarType scalar, std::size_t index) {
    // Check `cont.parameters()` returns something like an array, or vector, or something
    { scalar = cont.parameters()[index] };
    { cont.referenceEnergy() } -> std::same_as<ScalarType>;
    { cont.lowerEnergy() } -> std::same_as<ScalarType>;
    { cont.upperEnergy() } -> std::same_as<ScalarType>;
    { cont.type() } -> std::same_as<PeakContinuum::OffsetType>;
    { cont.externalContinuum() } -> std::same_as<std::shared_ptr<const SpecUtils::Measurement>>;
  };


  // This function is just templated version of `PeakContinuum::offset_integral(...)` - need to refactor
  //  both to use the same code
  template<typename ContType, typename ScalarType>
  void offset_integral( const ContType &cont,
                  const float *energies,
                  ScalarType *channels,
                  const size_t nchannel,
                  const std::shared_ptr<const SpecUtils::Measurement> &data ) requires ContinuumTypeConcept<ContType,ScalarType>;
#else

// Helper traits to check if a type satisfies certain conditions
template <typename ContType, typename ScalarType>
struct ContinuumTypeConcept {
private:
    template <typename T>
    static auto check_parameters(T* cont, std::size_t index) -> decltype((*cont).parameters()[index], std::true_type{});

    template <typename T>
    static auto check_referenceEnergy(T* cont) -> decltype((*cont).referenceEnergy(), std::true_type{});

    template <typename T>
    static auto check_lowerEnergy(T* cont) -> decltype((*cont).lowerEnergy(), std::true_type{});

    template <typename T>
    static auto check_upperEnergy(T* cont) -> decltype((*cont).upperEnergy(), std::true_type{});

    template <typename T>
    static auto check_type(T* cont) -> decltype((*cont).type(), std::true_type{});

    template <typename T>
    static auto check_externalContinuum(T* cont) -> decltype((*cont).externalContinuum(), std::true_type{});

public:
    static constexpr bool value =
        std::is_same<decltype(check_parameters(static_cast<ContType*>(nullptr), std::size_t{})), std::true_type>::value &&
        std::is_same<decltype(check_referenceEnergy(static_cast<ContType*>(nullptr))), std::true_type>::value &&
        std::is_same<decltype(check_lowerEnergy(static_cast<ContType*>(nullptr))), std::true_type>::value &&
        std::is_same<decltype(check_upperEnergy(static_cast<ContType*>(nullptr))), std::true_type>::value &&
        std::is_same<decltype(check_type(static_cast<ContType*>(nullptr))), std::true_type>::value &&
        std::is_same<decltype(check_externalContinuum(static_cast<ContType*>(nullptr))), std::true_type>::value;
};


// Enable the function only if the ContinuumTypeConcept is satisfied
template <typename ContType, typename ScalarType>
typename std::enable_if<ContinuumTypeConcept<ContType, ScalarType>::value, void>::type
offset_integral(const ContType& cont,
                const float* energies,
                ScalarType* channels,
                const size_t nchannel,
                const std::shared_ptr<const SpecUtils::Measurement>& data);
#endif //__cplusplus >= 202002L
}//namespace PeakDists

#endif  //PeakDists_h
