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

/** The functions in this .h/.cpp are for computing skewed and Gaussian photopeak distributions.
 
 */
namespace PeakDists
{
  /** `erfc` implementation based on boost, but sped up to be more efficient.
   
   Surprisingly, the erf() function is the major bottleneck for peak fitting.
   
   20191230: wcjohns extracted the boost::math::erf() function implementation
   from boost 1.65.1 for double precision (53 bit mantissa) into this function,
   boost_erf_imp(). Removing some of the supporting code structure, and
   explicitly writing out the polynomial equation evaluation seems to speed
   things up by about a factor of ~3 over calling boost::math::erf().
   
   In the implementation file, there is a commented out erf_approx() function looks to be about 25% faster than
   this boost version, but I havent carefully checked out the precision implications
   so not switching to it yet.
   */
  double boost_erf_imp( double z );
  
  /** Returns `1-boost_erf_imp(z)` */
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
  void photopeak_function_integral( const double mean,
                                          const double sigma,
                                          const double amplitude,
                                          const PeakDef::SkewType skew_type,
                                          const double * const skew_parameters,
                                          const size_t nchannel,
                                          const float * const lower_energies,
                                          double *peak_count_channels );
  
  
  
  
  
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
  void gaussian_integral( const double peak_mean,
                              const double peak_sigma,
                              const double peak_amplitude,
                              const float * const energies,
                              double *channels,
                              const size_t nchannel );
  
  
  
  /** Returns the integral of a unit-area Bortel function, between `x1` and `x2`. */
  double bortel_integral( const double mean, const double sigma, const double skew,
                         const double x1, const double x2 );
  
  /** Slightly CPU optimized method of computing the peak area, for Bortel skew over a number of channels.
   
   Cuts the number of calls to the `erf`, `erfc`, and exp functionsn half.
   Also, only calculates values between -12 to +8 sigma of the mean/
   
   @param peak_mean
   @param peak_sigma
   @param peak_amplitude
   @param skew The Bortel skew of the peak.
   @param energies Array of lower channel energies; must have at least one more entry than
          `nchannel`
   @param channels Channel count array integrals of Gaussian and Skew will be _added_ to (e.g.,
          will not be zeroed); must have at least `nchannel` entries
   @param nchannel The number of channels to do the integration over.
   */
  void bortel_integral( const double peak_mean,
                              const double peak_sigma,
                              const double peak_amplitude,
                              const double skew,
                              const float * const energies,
                              double *channels,
                              const size_t nchannel );
  
  /** Returns the PDF for a unit-area Bortel function.
   */
  double bortel_pdf( const double mean, const double sigma, const double skew_low, const double x );

  /** Returns the indefinite integral (e.g., from negative infinity up to `x`) of a unit-area Bortel function
   */
  double bortel_indefinite_integral( const double x, const double mean,
                                    const double sigma, const double skew );

  /** Returns an approximate area between `x1` and `x2` for a unit-area Bortel function.
   
   Just multiplies the x-range by the Bortel PDF value in the middle of the range.
   */
  double bortel_integral_fast( const double mean, const double sigma, const double skew,
                              const double x1, const double x2 );
  
  
  
  
  
  double gauss_exp_integral( const double peak_mean,
                                      const double peak_sigma,
                                      const double skew,
                                      const double x0, const double x1 );
  
  void gauss_exp_integral( const double peak_mean,
                               const double peak_sigma,
                               const double peak_amplitude,
                               const double skew,
                               const float * const energies,
                               double *channels,
                               const size_t nchannel );
  
  /** Returns the normalization so the GaussExp distribution has unit area. */
  double gauss_exp_norm( const double sigma, const double skew );
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
  
  
  
  
  
  
  
  
  
  double exp_gauss_exp_integral( const double mean,
                           const double sigma,
                           const double skew_left,
                           const double skew_right,
                           const double x0,
                                const double x1 );
  
  void exp_gauss_exp_integral( const double peak_mean,
                               const double peak_sigma,
                               const double peak_amplitude,
                               const double skew_left,
                               const double skew_right,
                               const float * const energies,
                               double *channels,
                               const size_t nchannel );
  
  double exp_gauss_exp_norm( const double sigma, const double skew_left, const double skew_right );
  
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
  
  
  
  
  
  
  
  double crystal_ball_integral( const double peak_mean,
                                      const double peak_sigma,
                                      const double alpha,
                                      const double power_law,
                                      const double x0, const double x1 );
  
  void crystal_ball_integral( const double peak_mean,
                               const double peak_sigma,
                               const double peak_amplitude,
                               const double alpha,
                               const double power_law,
                               const float * const energies,
                               double *channels,
                               const size_t nchannel );
  
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
  
  
  
  
  
  
  
  
  
  
  
  double double_sided_crystal_ball_integral( const double peak_mean,
                                                   const double peak_sigma,
                                                   const double lower_alpha,
                                                   const double lower_power_law,
                                                   const double upper_alpha,
                                                   const double upper_power_law,
                                                   const double x0, const double x1 );
  
  void double_sided_crystal_ball_integral( const double peak_mean,
                                                 const double peak_sigma,
                                                 const double peak_amplitude,
                                                 const double lower_alpha,
                                                 const double lower_power_law,
                                                 const double upper_alpha,
                                                 const double upper_power_law,
                                                 const float * const energies,
                                                 double *channels,
                                                 const size_t nchannel );
    
    

    
    /** Returns the non-normalized, double sided Crystal Ball PDF value for `x` */
    double DSCB_pdf_non_norm(const double mean,
                                      const double sigma,
                                      const double alpha_low,
                                      const double n_low,
                                      const double alpha_high,
                                      const double n_high,
                                         const double x );

    double DSCB_norm( const double alpha_low,
                                          const double n_low,
                                          const double alpha_high,
                                          const double n_high );
    
    double DSCB_left_tail_indefinite_non_norm_t( const double alpha_low, const double n_low,
                                                 const double t );
    
    double DSCB_right_tail_indefinite_non_norm_t( const double alpha_high,
                                                             const double n_high,
                                                             const double t);
    
    double DSCB_gauss_indefinite_non_norm_t( const double t );
  
}//namespace PeakDists

#endif  //PeakDists_h
