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
 * @def USE_INACCURATE_VOIGT_CDF
 * @brief Enable/disable the Voigt CDF-based functions (default: 0 = disabled)
 * 
 * The exact analytical CDF formula for the Voigt profile (using erf(z) for complex z)
 * has proven problematic - Re[erf(z)] can fall outside [-1, 1] for complex arguments,
 * leading to invalid CDF values (negative or >1) that require clamping, which causes
 * large errors in integral calculations. We were unable to get the exact CDF formula
 * to work reliably, so we default to numerical PDF integration instead, which is slower
 * but more accurate.
 * 
 * When USE_INACCURATE_VOIGT_CDF is 0 (default):
 *   - voigt_cdf() is not available
 *   - voigt_exp_indefinite() is not available (uses voigt_cdf internally)
 *   - voigt_exp_integral() scalar version is not available (uses voigt_exp_indefinite)
 *   - voigt_exp_integral() array version is not available (uses voigt_exp_indefinite)
 *   - voigt_exp_coverage_limits() is not available (uses voigt_exp_indefinite)
 * 
 * When USE_INACCURATE_VOIGT_CDF is 1:
 *   - All CDF-based functions are enabled, but may have accuracy issues
 *   - May revisit in the future to try alternative CDF formulations
 */
#ifndef USE_INACCURATE_VOIGT_CDF
#define USE_INACCURATE_VOIGT_CDF 0
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
    // NOTE: Matching test_voigt_exact.cpp exactly to debug discrepancy with scipy
    T sqrt2 = T(1.41421356237309504880);
    T sqrt_pi = T(1.77245385090551602730);
    T sigma_sqrt2 = sigma * sqrt2;
    T z_real = (x - mean) / sigma_sqrt2;
    T z_imag = gamma / sigma_sqrt2;

    // w(z) = Faddeeva function
    Complex<T> z_complex(z_real, z_imag);
    Complex<T> w_val = FaddeevaT::w(z_complex);

    // Voigt PDF = Re[w(z)] / (σ * √π * √2) = Re[w(z)] / (σ√(2π))
    // This matches test_voigt_exact.cpp line 29 and scipy formula
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

// ============================================================================
// Functions gated by USE_INACCURATE_VOIGT_CDF
// ============================================================================

#if USE_INACCURATE_VOIGT_CDF
/**
 * Numerically stable Voigt CDF.
 * Uses the Faddeeva erfcx based formulation for better stability than erf(z).
 * 
 * WARNING: This function has known accuracy issues. The exact analytical CDF formula
 * (using erf(z) for complex z) can produce invalid values (negative or >1) that require
 * clamping, leading to errors in integral calculations.
 */
template<typename T>
T voigt_cdf(const T x, const T mean, const T sigma, const T gamma) {
    if (gamma <= T(1e-11) * sigma) return gaussian_cdf(x, mean, sigma);
    if (sigma <= T(1e-11) * gamma) return lorentzian_cdf(x, mean, gamma);

    T sqrt2 = T(1.41421356237309504880);
    T z_re = (x - mean) / (sigma * sqrt2);
    T z_im = gamma / (sigma * sqrt2);

    // To avoid the exponential explosion in erf(z) = 1 - exp(-z^2)w(iz),
    // we use Faddeeva's w(z) directly. The CDF of a Voigt profile is:
    // CDF(x) = 0.5 + 0.5 * Re[erf(z)]
    // We utilize FaddeevaT::erf which is implemented in your Faddeeva_impl.hpp.
    // However, if gamma is large, erf(z) is unstable. 
    // We check the scalar value of z_im to decide on strategy.
    
    if (get_scalar(z_im) > 8.0) {
        // High-gamma limit (Lorentzian dominated)
        return lorentzian_cdf(x, mean, gamma);
    }

    Complex<T> z(z_re, z_im);
    Complex<T> res = FaddeevaT::erf(z);

    T val = T(0.5) + T(0.5) * res.real;
    
    // Final safety clamps
    if (val < T(0)) return T(0);
    if (val > T(1)) return T(1);
    return val;
}

/** Returns the indefinite integral (from -infinity to x) of a unit-area Voigt with exponential tail.

   @param x The upper limit of integration
   @param mean The peak mean in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau

   @returns Cumulative distribution value at x
   
   WARNING: This function uses voigt_cdf() which has known accuracy issues.
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
   
   WARNING: This function uses voigt_exp_indefinite() which relies on the inaccurate voigt_cdf().
   */
inline double voigt_exp_integral(const double peak_mean, const double sigma_gauss,
                             const double gamma_lor, const double tail_ratio,
                          const double tail_slope, const double x0, const double x1) {
    double cdf_x1 = voigt_exp_indefinite(x1, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    double cdf_x0 = voigt_exp_indefinite(x0, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    return cdf_x1 - cdf_x0;
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
   
   WARNING: This function uses voigt_exp_indefinite() which relies on the inaccurate voigt_cdf().
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

/** Optimized array-filling version of the Voigt with exponential tail integral.

   @param peak_mean The peak mean in keV
   @param sigma_gauss Gaussian width (from detector resolution), in keV
   @param peak_amplitude The peak amplitude (use 1.0 for unit-area peak)
   @param gamma_lor Lorentzian HWHM (from natural line width), in keV
   @param tail_ratio Fraction of counts in the exponential tail
   @param tail_slope Exponential tail slope parameter tau
   @param energies Array of channel lower energies (must have nchannel+1 entries)
   @param channels Array where distribution values will be added (must have nchannel entries)
   @param nchannel Number of channels to integrate over
   
   WARNING: This function uses voigt_exp_indefinite() which relies on the inaccurate voigt_cdf().
   */
template<typename T>
void voigt_exp_integral(const T peak_mean, const T sigma_gauss,
                        const T peak_amplitude, const T gamma_lor,
                        const T tail_ratio, const T tail_slope,
                        const float * const energies, T *channels,
                        const size_t nchannel) {
    // CDF-based implementation using voigt_exp_indefinite
    std::vector<T> cdf_values(nchannel + 1);
    for (size_t i = 0; i <= nchannel; ++i) {
        T x = T(energies[i]);
        cdf_values[i] = voigt_exp_indefinite(x, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    }

    for (size_t i = 0; i < nchannel; ++i) {
        T bin_content = (cdf_values[i + 1] - cdf_values[i]) * peak_amplitude;
        if (bin_content < T(0)) bin_content = T(0);
        channels[i] += bin_content;
    }
}
#endif // USE_INACCURATE_VOIGT_CDF

#endif // VOIGT_EXP_TAIL_HPP