#ifndef VOIGT_EXP_TAIL_HPP
#define VOIGT_EXP_TAIL_HPP

#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cstddef>
#include <utility>
#include <vector>
#include "Faddeeva.hpp"
#if __has_include(<ceres/jet.h>)
#include <ceres/jet.h>
#endif

/**
 * @brief Implements the Voigt distribution with an exponential tail (Hypermet), in a templated manner.
  
 Note: unlike other distributions in InterSpec, or even pseudo-Voigt distribution, this
 implementation does NOT use the CDF to fill the array of channels. Instead, it uses the PDF,
 with a Gauss-Legendre quadrature to integrate the PDF over each channel.  This is because 
 implementing the symbolically integrated PDF in a numerically stanble way is not straightforward, 
 and I couldnt get a reasonable result with the time avaiable.

 A side-effect of this is this file does not yet imeplent the `voigt_exp_coverage_limits(...)`
 function - you can use `pseudo_voigt::voigt_exp_coverage_limits(...)` as a decent approximation.
 */

namespace voigt_exp_tail {

// ADL helpers for std / ceres math so Jets keep derivatives
namespace voigt_adl {
template <typename T>
inline T adl_exp(const T& v) {
    using std::exp;
#if __has_include(<ceres/jet.h>)
    using ceres::exp;
#endif
    return exp(v);
}

template <typename T>
inline T adl_sqrt(const T& v) {
    using std::sqrt;
#if __has_include(<ceres/jet.h>)
    using ceres::sqrt;
#endif
    return sqrt(v);
}

template <typename T, typename U>
inline auto adl_pow(const T& a, const U& b) {
    using std::pow;
#if __has_include(<ceres/jet.h>)
    using ceres::pow;
#endif
    return pow(a, b);
}

template <typename T>
inline T adl_atan(const T& v) {
    using std::atan;
#if __has_include(<ceres/jet.h>)
    using ceres::atan;
#endif
    return atan(v);
}

template <typename T>
inline T adl_abs(const T& v) {
    using std::abs;
#if __has_include(<ceres/jet.h>)
    using ceres::abs;
#endif
    return abs(v);
}
} // namespace voigt_adl

// ============================================================================
// Gaussian distribution functions
// ============================================================================

template<typename T>
T gaussian_pdf(const T x, const T mean, const T sigma) {
    const T sqrt2pi = T(2.5066282746310002);
    T diff = (x - mean) / sigma;
    return voigt_adl::adl_exp(-T(0.5) * diff * diff) / (sigma * sqrt2pi);
}

template<typename T>
T gaussian_cdf(const T x, const T mean, const T sigma) {
    const T sqrt2 = T(1.41421356237309504880);
    T z = (x - mean) / (sigma * sqrt2);
    return T(0.5) * (T(1) + FaddeevaT::erf_real(z));
}

// ============================================================================
// Lorentzian distribution functions
// ============================================================================

template<typename T>
T lorentzian_pdf(const T x, const T mean, const T gamma) {
    // Lorentzian (Cauchy) PDF with HWHM gamma
    const T pi = T(3.14159265358979323846);
    T diff = x - mean;
    return (gamma / pi) / (diff * diff + gamma * gamma);
}

template<typename T>
T lorentzian_cdf(const T x, const T mean, const T gamma) {
    // Lorentzian (Cauchy) CDF with HWHM gamma
    const T pi = T(3.14159265358979323846);
    return T(0.5) + voigt_adl::adl_atan((x - mean) / gamma) / pi;
}

// ============================================================================
// Voigt profile functions using true Faddeeva implementation
// ============================================================================

template<typename T>
T voigt_mixing_eta(const T sigma, const T gamma) {
    // Pseudo-Voigt mixing parameter using Thompson-Cox-Hastings formula
    // This approximates the Voigt as η*Lorentzian + (1-η)*Gaussian
    // f_L = 2*gamma is the Lorentzian FWHM
    // f_G = 2*sigma*sqrt(2*ln(2)) is the Gaussian FWHM
    const T ln2 = T(0.6931471805599453);
    T f_G = T(2) * sigma * voigt_adl::adl_sqrt(T(2) * ln2);  // Gaussian FWHM
    T f_L = T(2) * gamma;  // Lorentzian FWHM

    // Approximate total FWHM using empirical formula
    T f5 = voigt_adl::adl_pow(f_G, 5) + T(2.69269) * voigt_adl::adl_pow(f_G, 4) * f_L
         + T(2.42843) * voigt_adl::adl_pow(f_G, 3) * voigt_adl::adl_pow(f_L, 2)
         + T(4.47163) * voigt_adl::adl_pow(f_G, 2) * voigt_adl::adl_pow(f_L, 3)
         + T(0.07842) * f_G * voigt_adl::adl_pow(f_L, 4) + voigt_adl::adl_pow(f_L, 5);
    T f_V = voigt_adl::adl_pow(f5, T(0.2));  // Total Voigt FWHM

    // Mixing parameter
    T eta = T(1.36603) * (f_L / f_V) - T(0.47719) * voigt_adl::adl_pow(f_L / f_V, 2)
          + T(0.11116) * voigt_adl::adl_pow(f_L / f_V, 3);

    // Clamp to [0, 1]
    if (eta < T(0)) eta = T(0);
    if (eta > T(1)) eta = T(1);

    return eta;
}

template<typename T>
T voigt_pdf(const T x, const T mean, const T sigma, const T gamma) {
    // True Voigt profile using Faddeeva function w(z)
    // V(x) = (1/(σ√(2π))) * Re[w(z)] where z = ((x - μ) + iγ)/(σ√2)

    // Special case: pure Gaussian (gamma = 0)
    if (gamma <= T(1e-10) * sigma) {
        return gaussian_pdf(x, mean, sigma);
    }

    // Special case: pure Lorentzian (sigma very small)
    if (sigma <= T(1e-10) * gamma) {
        return lorentzian_pdf(x, mean, gamma);
    }

    // Compute z = ((x - mean) + i*gamma) / (sigma * sqrt(2))
    T sqrt2 = T(1.41421356237309504880);
    T sqrt_pi = T(1.77245385090551602730);
    T sigma_sqrt2 = sigma * sqrt2;
    T z_real = (x - mean) / sigma_sqrt2;
    T z_imag = gamma / sigma_sqrt2;

    // w(z) = Faddeeva function
    Complex<T> z_complex(z_real, z_imag);
    Complex<T> w_val = FaddeevaT::w(z_complex);

    // Voigt PDF = Re[w(z)] / (σ * √π * √2) = Re[w(z)] / (σ√(2π))
    // This matches test_voigt_against_mit_imp.cpp line 29 and scipy formula
    return w_val.real / (sigma * sqrt_pi * sqrt2);
}

// ============================================================================
// Exponentially-modified Gaussian (low-energy tail) functions
// ============================================================================

template<typename T>
T gaussexp_pdf(const T x, const T mean, const T sigma, const T tau) {
    // Exponentially-modified Gaussian PDF (low-energy/left tail)
    // f(x) = (λ/2) * exp(λ(μ - x) + (λσ)^2 / 2) *
    //        erfc( (μ - x)/(σ√2) + λσ/√2 )
    // where λ = 1/τ and τ > 0.
    const T sqrt2 = T(1.41421356237309504880);

    if (tau <= T(0)) {
        // Degenerate to pure Gaussian
        T diff = (x - mean) / sigma;
        return T(1) / (sigma * sqrt2 * T(1.77245385090551602730)) * voigt_adl::adl_exp(-T(0.5) * diff * diff);
    }

    const T lambda = T(1) / tau;
    const T lambda_sigma = lambda * sigma;
    T arg = (mean - x) / (sigma * sqrt2) + lambda_sigma / sqrt2;
    T erfc_val = FaddeevaT::erfc_real(arg);
    T exp_arg = lambda * (mean - x) + T(0.5) * lambda_sigma * lambda_sigma;
    T exp_val = voigt_adl::adl_exp(exp_arg);

    return lambda * T(0.5) * exp_val * erfc_val;
}

template<typename T>
T gaussexp_cdf(const T x, const T mean, const T sigma, const T tau) {
    // Exponentially-modified Gaussian CDF (low-energy/left tail)
    // Derived by mirroring the standard right-tail EMG:
    // Let a = (μ - x)/σ and λ = 1/τ.
    // F(x) = 1 - Φ(a) + exp(-λ(μ - x) + (λσ)^2 / 2) * Φ(a - λσ)
    // This is monotone and stays in [0, 1] for τ > 0.
    const T sqrt2 = T(1.41421356237309504880);

    if (tau <= T(0)) {
        // Degenerate to pure Gaussian CDF
        T diff = (x - mean) / (sigma * sqrt2);
        return T(0.5) * (T(1) + FaddeevaT::erf_real(diff));
    }

    const T lambda = T(1) / tau;
    const T lambda_sigma = lambda * sigma;

    T a = (mean - x) / sigma;
    T phi_a = T(0.5) * (T(1) + FaddeevaT::erf_real(a / sqrt2));          // Φ(a)
    T phi_shift = T(0.5) * (T(1) + FaddeevaT::erf_real((a - lambda_sigma) / sqrt2)); // Φ(a - λσ)
    T exp_arg = -lambda * (mean - x) + T(0.5) * lambda_sigma * lambda_sigma;
    // Guard against overflow/underflow in exp term
    if (exp_arg > T(700)) {
        return T(1); // far right tail, CDF ~ 1
    }
    if (exp_arg < T(-700)) {
        T res = T(1) - phi_a; // exponential term vanishes
        if (res < T(0)) res = T(0);
        return res;
    }
    T exp_term = voigt_adl::adl_exp(exp_arg);

    T result = T(1) - phi_a + exp_term * phi_shift;

    // Guard against tiny numerical drift
    if (result < T(0)) result = T(0);
    if (result > T(1)) result = T(1);

    return result;
}


// ============================================================================
// VoigtExpTail API functions
// ============================================================================

/** Returns the PDF (probability density function) of the Voigt with exponential tail distribution.

   The VoigtExpTail distribution combines a Voigt profile (convolution of Gaussian and
   Lorentzian) with an exponential low-energy tail, used for modeling x-ray peaks measured
   with HPGe detectors where incomplete charge collection creates a tail.

   @param x Energy value at which to evaluate the PDF
   @param mean Peak centroid energy (mean of the distribution), in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail (0 to ~0.3)
   @param tail_slope Exponential tail slope parameter tau, similar to GaussExp skew

   @returns PDF value at the given energy
   */
template<typename T>
T voigt_exp_pdf(const T x, const T mean, const T sigma_gauss, const T gamma_lor,
                const T tail_ratio, const T tail_slope) {
    // Combined PDF: (1-tail_ratio) * Voigt + tail_ratio * GaussExp
    T voigt_pdf_val = voigt_pdf(x, mean, sigma_gauss, gamma_lor);
    T gaussexp_pdf_val = gaussexp_pdf(x, mean, sigma_gauss, tail_slope);

    return (T(1) - tail_ratio) * voigt_pdf_val + tail_ratio * gaussexp_pdf_val;
}

/** Array-filling version of the Voigt with exponential tail integral for each channel.
   Uses Gauss-Legendre quadrature to integrate the PDF over each channel.

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
  void voigt_exp_integral(const T peak_mean, const T sigma_gauss,
                          const T peak_amplitude, const T gamma_lor,
                          const T tail_ratio, const T tail_slope,
                          const float * const energies, T *channels,
                          const size_t nchannel) 
  {
      // Determine characteristic width for adaptive quadrature
      T max_width = (sigma_gauss > gamma_lor) ? sigma_gauss : gamma_lor;
      T threshold = T(0.25) * max_width;

      // 3-point Gauss-Legendre quadrature nodes and weights on [-1, 1]
      const T gl3_nodes[3] = {
          T(-0.77459666924148337704),  // -sqrt(3/5)
          T(0.0),
          T(0.77459666924148337704)    // sqrt(3/5)
      };
      const T gl3_weights[3] = {
          T(0.55555555555555555556),   // 5/9
          T(0.88888888888888888889),   // 8/9
          T(0.55555555555555555556)    // 5/9
      };

      // 5-point Gauss-Legendre quadrature nodes and weights on [-1, 1]
      const T gl5_nodes[5] = {
          T(-0.90617984593866399280),  // -sqrt(5+2*sqrt(10/7))/3
          T(-0.53846931010568309104),  // -sqrt(5-2*sqrt(10/7))/3
          T(0.0),
          T(0.53846931010568309104),   // sqrt(5-2*sqrt(10/7))/3
          T(0.90617984593866399280)    // sqrt(5+2*sqrt(10/7))/3
      };
      const T gl5_weights[5] = {
          T(0.23692688505618908751),   // (322-13*sqrt(70))/900
          T(0.47862867049936646804),   // (322+13*sqrt(70))/900
          T(0.56888888888888888889),   // 128/225
          T(0.47862867049936646804),   // (322+13*sqrt(70))/900
          T(0.23692688505618908751)    // (322-13*sqrt(70))/900
      };

      for (size_t i = 0; i < nchannel; ++i) {
          T a = T(energies[i]);      // Lower edge
          T b = T(energies[i + 1]);  // Upper edge
          T bin_width = b - a;

          // Choose number of quadrature points based on bin width
          int npoints;
          const T* nodes;
          const T* weights;
          
          if (bin_width < threshold) {
              // Narrow bins: use 3-point quadrature
              npoints = 3;
              nodes = gl3_nodes;
              weights = gl3_weights;
          } else {
              // Wide bins: use 5-point quadrature
              npoints = 5;
              nodes = gl5_nodes;
              weights = gl5_weights;
          }

          // Transform from [-1, 1] to [a, b]
          // x = (b-a)/2 * xi + (a+b)/2
          T half_width = T(0.5) * bin_width;
          T center = T(0.5) * (a + b);

          // Integrate PDF over [a, b] using Gauss-Legendre quadrature
          T integral = T(0);
          for (int j = 0; j < npoints; ++j) {
              T x = half_width * nodes[j] + center;
              T pdf_val = voigt_exp_pdf(x, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
              integral += weights[j] * pdf_val;
          }

          // Scale by bin width and amplitude
          // The weights are for integration on [-1, 1], so we multiply by (b-a)/2
          T bin_content = integral * half_width * peak_amplitude;
          if (bin_content < T(0)) bin_content = T(0);
          channels[i] += bin_content;
      }
  }



} // namespace voigt_exp_tail

#endif // VOIGT_EXP_TAIL_HPP