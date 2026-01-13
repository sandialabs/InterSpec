#ifndef VOIGT_EXP_TAIL_HPP
#define VOIGT_EXP_TAIL_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cstddef>
#include <utility>
#include "Faddeeva.hpp"
#if __has_include(<ceres/jet.h>)
#include <ceres/jet.h>
#endif

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
    T z_real = (x - mean) / (sigma * sqrt2);
    T z_imag = gamma / (sigma * sqrt2);

    // w(z) = Faddeeva function
    Complex<T> z_complex(z_real, z_imag);
    Complex<T> w_val = FaddeevaT::w(z_complex);

    // Voigt PDF = (1/(σ√(2π))) * Re[w(z)]
    T sqrt_2pi = T(2.50662827463100050242); // √(2π)
    return w_val.real / (sigma * sqrt_2pi);
}

template<typename T>
T voigt_cdf(const T x, const T mean, const T sigma, const T gamma) {
    // True Voigt CDF using complex erf (like the exact test):
    // F_V(x) = 0.5 + 0.5 * Re[erf(z)] where z = ((x - μ) + iγ)/(σ√2)

    // Special case: pure Gaussian (gamma = 0)
    if (gamma <= T(1e-10) * sigma) {
        return gaussian_cdf(x, mean, sigma);
    }

    // Special case: pure Lorentzian (sigma very small) or Lorentzian-dominated
    // When gamma/sigma is large (> 5), the Voigt is Lorentzian-dominated and
    // the erf formula becomes numerically unstable. Use Lorentzian CDF instead.
    if (sigma <= T(1e-10) * gamma || gamma > T(5.0) * sigma) {
        return lorentzian_cdf(x, mean, gamma);
    }

    // Compute z = ((x - mean) + i*gamma) / (sigma * sqrt(2))
    T sqrt2 = T(1.41421356237309504880);
    T z_real = (x - mean) / (sigma * sqrt2);
    T z_imag = gamma / (sigma * sqrt2);

    // Check if z_imag is too large (would cause numerical issues in erf)
    // When z_imag > 10, erf(z) becomes numerically unstable
    double z_imag_scalar = get_scalar(z_imag);
    if (z_imag_scalar > 10.0) {
        // Use Lorentzian CDF as approximation for Lorentzian-dominated case
        return lorentzian_cdf(x, mean, gamma);
    }

    // erf(z) = complex error function
    Complex<T> z_complex(z_real, z_imag);
    Complex<T> erf_result = FaddeevaT::erf(z_complex);

    // Compute CDF and clamp to [0, 1] to handle numerical issues
    T cdf = T(0.5) + T(0.5) * erf_result.real;
    
    // Clamp to valid CDF range [0, 1] to handle numerical overflow/underflow
    if (cdf < T(0.0)) cdf = T(0.0);
    if (cdf > T(1.0)) cdf = T(1.0);
    
    return cdf;
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

/** Returns the normalization for a unit-area Voigt with exponential tail distribution.

   The VoigtExpTail distribution combines a Voigt profile (convolution of Gaussian and
   Lorentzian) with an exponential low-energy tail, used for modeling x-ray peaks measured
   with HPGe detectors where incomplete charge collection creates a tail.

   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail (0 to ~0.3)
   @param tail_slope Exponential tail slope parameter tau, similar to GaussExp skew

   @returns Normalization constant so the distribution integrates to 1
   */
template<typename T>
T voigt_exp_norm(const T sigma_gauss, const T gamma_lor, const T tail_ratio, const T tail_slope) {
    // Both Voigt and GaussExp are unit-area distributions, so the combination
    // (1 - tail_ratio) * Voigt + tail_ratio * GaussExp is also unit-area
    return T(1.0);
}

/** Returns the indefinite integral (from -infinity to x) of a unit-area Voigt with exponential tail.

   @param x The upper limit of integration
   @param mean The peak mean in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau

   @returns Cumulative distribution value at x
   */
template<typename T>
T voigt_exp_indefinite(const T x, const T mean, const T sigma_gauss,
                       const T gamma_lor, const T tail_ratio, const T tail_slope) {
    T voigt_cdf_val = voigt_cdf(x, mean, sigma_gauss, gamma_lor);
    T gaussexp_cdf_val = gaussexp_cdf(x, mean, sigma_gauss, tail_slope);
    return (T(1) - tail_ratio) * voigt_cdf_val + tail_ratio * gaussexp_cdf_val;
}

// Non-templated version for compatibility
inline double voigt_exp_indefinite(const double x, const double mean, const double sigma_gauss,
                            const double gamma_lor, const double tail_ratio, const double tail_slope) {
    return voigt_exp_indefinite<double>(x, mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
}

/** Returns the integral of a Voigt with exponential tail distribution between x0 and x1.

   @param peak_mean The peak mean in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau
   @param x0 Lower integration limit
   @param x1 Upper integration limit

   @returns Integral of the unit-area distribution from x0 to x1
   */
inline double voigt_exp_integral(const double peak_mean, const double sigma_gauss,
                             const double gamma_lor, const double tail_ratio,
                          const double tail_slope, const double x0, const double x1) {
    double cdf_x1 = voigt_exp_indefinite(x1, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    double cdf_x0 = voigt_exp_indefinite(x0, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    return cdf_x1 - cdf_x0;
}

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
void voigt_exp_integral(const T peak_mean, const T sigma_gauss,
                           const T peak_amplitude, const T gamma_lor,
                           const T tail_ratio, const T tail_slope,
                           const float * const energies, T *channels,
                        const size_t nchannel) {
    // Evaluate CDF at all bin edges
    std::vector<T> cdf_values(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        T x = T(energies[i]);
        cdf_values[i] = voigt_exp_indefinite(x, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    }

    // Compute bin contents as differences
    for (size_t i = 0; i < nchannel; ++i) {
        T bin_content = (cdf_values[i + 1] - cdf_values[i]) * peak_amplitude;
        // Ensure non-negative (clamp to zero for numerical precision issues)
        if (bin_content < T(0)) bin_content = T(0);
        channels[i] += bin_content;
    }
}

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

   Throws error on invalid input.
   */
inline std::pair<double, double> voigt_exp_coverage_limits(const double mean, const double sigma_gauss,
                                                       const double gamma_lor, const double tail_ratio,
                                                     const double tail_slope, const double p) {
    if (p < 0.0 || p >= 1.0) {
        throw std::invalid_argument("p must be in [0, 1)");
    }

    double target_lower = 0.5 * p;
    double target_upper = 1.0 - 0.5 * p;

    // Binary search for lower limit
    // Use a width that accounts for Gaussian, Lorentzian, and the exponential tail.
    // - Gaussian: ~20*sigma covers essentially all mass
    // - Lorentzian: heavy tails decay slowly; derive width from desired tail probability p
    //   For a Cauchy/Lorentzian, width to leave p/2 mass in each tail is:
    //     width = gamma * cot(pi * p / 2) ≈ 2*gamma / (pi*p) for small p
    // - Exponential tail: ~tail_slope * log(1/p) to include tail mass
    double tail_extent = (tail_ratio > 0.0 && tail_slope > 0.0) ? tail_slope * std::log(1.0 / p) : 0.0;

    // Lorentzian-derived width for target p (guard small p with asymptotic)
    double lorentz_width = 0.0;
    if (gamma_lor > 0.0) {
        if (p > 0.0 && p < 1e-3) {
            lorentz_width = (2.0 * gamma_lor) / (M_PI * p); // asymptotic cot(pi*p/2)
        } else if (p > 0.0) {
            lorentz_width = gamma_lor / std::tan(0.5 * M_PI * p);
        } else { // p == 0: impossible; choose a large multiple of gamma
            lorentz_width = 1e6 * gamma_lor;
        }
        // Ensure at least 20*gamma for moderate tails
        lorentz_width = std::max(lorentz_width, 20.0 * gamma_lor);
    }

    double search_width = std::max({20.0 * sigma_gauss, lorentz_width, tail_extent});

    double lower = mean - search_width;
    double upper = mean;
    for (int iter = 0; iter < 100; ++iter) {
        double mid = 0.5 * (lower + upper);
        double cdf_mid = voigt_exp_indefinite(mid, mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
        if (std::abs(cdf_mid - target_lower) < 1e-10) {
            lower = mid;
            break;
        }
        if (cdf_mid < target_lower) {
            lower = mid;
        } else {
            upper = mid;
        }
    }

    // Binary search for upper limit
    double lower2 = mean;
    double upper2 = mean + search_width;  // consistent width for upper tail
    for (int iter = 0; iter < 100; ++iter) {
        double mid = 0.5 * (lower2 + upper2);
        double cdf_mid = voigt_exp_indefinite(mid, mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
        if (std::abs(cdf_mid - target_upper) < 1e-10) {
            upper2 = mid;
            break;
        }
        if (cdf_mid < target_upper) {
            lower2 = mid;
        } else {
            upper2 = mid;
        }
    }

    return std::make_pair(lower, upper2);
}

#endif // VOIGT_EXP_TAIL_HPP