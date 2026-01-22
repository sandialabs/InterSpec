// Pseudo-Voigt-based Voigt+Exponential Tail implementation
// This mirrors the API of voigt_exp_tail.hpp but uses the pseudo-Voigt
// approximation (Gaussian/Lorentzian mixture) for speed and simplicity.
// Only standard C++ (cmath) is used; no external deps.

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace pseudo_voigt {

// ADL helpers for std / ceres math so Jets keep derivatives
namespace pseudo_voigt_adl {
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
inline T adl_erf(const T& v) {
    using std::erf;
#if __has_include(<ceres/jet.h>)
    using ceres::erf;
#endif
    return erf(v);
}

template <typename T>
inline T adl_erfc(const T& v) {
    using std::erfc;
#if __has_include(<ceres/jet.h>)
    using ceres::erfc;
#endif
    return erfc(v);
}
} // namespace pseudo_voigt_adl

// ===============================================================
// Basic distributions
// ===============================================================

template<typename T>
inline T gaussian_pdf(const T x, const T mean, const T sigma) {
    const T sqrt_2pi = T(2.50662827463100050242); // sqrt(2*pi)
    T diff = (x - mean) / sigma;
    return exp(-T(0.5) * diff * diff) / (sigma * sqrt_2pi);
}

template<typename T>
inline T gaussian_cdf(const T x, const T mean, const T sigma) {
    const T sqrt2 = T(1.41421356237309504880);
    T z = (x - mean) / (sigma * sqrt2);
    return T(0.5) * (T(1) + erf(z));
}

template<typename T>
inline T lorentzian_pdf(const T x, const T mean, const T gamma) {
    const T pi = T(3.14159265358979323846);
    T diff = x - mean;
    return (gamma / pi) / (diff * diff + gamma * gamma);
}

template<typename T>
inline T lorentzian_cdf(const T x, const T mean, const T gamma) {
    const T pi = T(3.14159265358979323846);
    return T(0.5) + atan((x - mean) / gamma) / pi;
}

// ===============================================================
// Pseudo-Voigt helper (Thompson-Cox-Hastings)
// ===============================================================

template<typename T>
struct PseudoVoigtParams {
    T eta;
    T sigma_p;
    T gamma_p;
};

template<typename T>
inline PseudoVoigtParams<T> voigt_pseudo_params(const T sigma, const T gamma) {
    PseudoVoigtParams<T> out{};

    if (gamma <= T(1e-12) * sigma) {
        out.eta = T(0);
        out.sigma_p = sigma;
        out.gamma_p = gamma;
        return out;
    }
    if (sigma <= T(1e-12) * gamma) {
        out.eta = T(1);
        out.sigma_p = sigma;
        out.gamma_p = gamma;
        return out;
    }

    const T ln2 = T(0.6931471805599453);
    const T f_G = T(2) * sigma * pseudo_voigt_adl::adl_sqrt(T(2) * ln2);
    const T f_L = T(2) * gamma;

    const T f5 = pseudo_voigt_adl::adl_pow(f_G, 5) + T(2.69269) * pseudo_voigt_adl::adl_pow(f_G, 4) * f_L
               + T(2.42843) * pseudo_voigt_adl::adl_pow(f_G, 3) * pseudo_voigt_adl::adl_pow(f_L, 2)
               + T(4.47163) * pseudo_voigt_adl::adl_pow(f_G, 2) * pseudo_voigt_adl::adl_pow(f_L, 3)
               + T(0.07842) * f_G * pseudo_voigt_adl::adl_pow(f_L, 4) + pseudo_voigt_adl::adl_pow(f_L, 5);
    const T f_V = pseudo_voigt_adl::adl_pow(f5, T(0.2));

    const T ratio = f_L / f_V;
    T eta = T(1.36603) * ratio - T(0.47719) * ratio * ratio + T(0.11116) * ratio * ratio * ratio;
    if (eta < T(0)) eta = T(0);
    if (eta > T(1)) eta = T(1);

    const T sigma_p = f_V / (T(2) * pseudo_voigt_adl::adl_sqrt(T(2) * ln2));
    const T gamma_p = f_V / T(2);

    out.eta = eta;
    out.sigma_p = sigma_p;
    out.gamma_p = gamma_p;
    return out;
}

template<typename T>
inline T voigt_pdf_pseudo_precomputed(const T x, const T mean, const PseudoVoigtParams<T>& p) {
    return p.eta * lorentzian_pdf(x, mean, p.gamma_p) +
           (T(1) - p.eta) * gaussian_pdf(x, mean, p.sigma_p);
}

// ===============================================================
// Pseudo-Voigt PDF / CDF
// ===============================================================

template<typename T>
inline T voigt_pdf(const T x, const T mean, const T sigma, const T gamma) {
    // Pseudo-Voigt PDF
    if (gamma <= T(1e-10) * sigma) return gaussian_pdf(x, mean, sigma);
    if (sigma <= T(1e-10) * gamma) return lorentzian_pdf(x, mean, gamma);

    auto p = voigt_pseudo_params(sigma, gamma);
    return voigt_pdf_pseudo_precomputed(x, mean, p);
}

template<typename T>
inline T voigt_cdf(const T x, const T mean, const T sigma, const T gamma) {
    // Pseudo-Voigt CDF (exact for the Gaussian/Lorentzian mixture)
    if (gamma <= T(1e-10) * sigma) return gaussian_cdf(x, mean, sigma);
    if (sigma <= T(1e-10) * gamma) return lorentzian_cdf(x, mean, gamma);

    auto p = voigt_pseudo_params(sigma, gamma);
    T cdf = p.eta * lorentzian_cdf(x, mean, p.gamma_p) +
            (T(1) - p.eta) * gaussian_cdf(x, mean, p.sigma_p);
    // Clamp to [0,1] for numerical safety
    if (cdf < T(0)) cdf = T(0);
    if (cdf > T(1)) cdf = T(1);
    return cdf;
}

// ===============================================================
// Exponentially-modified Gaussian (left tail)
// ===============================================================

template<typename T>
inline T gaussexp_pdf(const T x, const T mean, const T sigma, const T tau) {
    const T sqrt2 = T(1.41421356237309504880);
    if (tau <= T(0)) {
        return gaussian_pdf(x, mean, sigma);
    }
    const T lambda = T(1) / tau;
    const T lambda_sigma = lambda * sigma;
    T arg = (mean - x) / (sigma * sqrt2) + lambda_sigma / sqrt2;
    T erfc_val = pseudo_voigt_adl::adl_erfc(arg);
    T exp_arg = lambda * (mean - x) + T(0.5) * lambda_sigma * lambda_sigma;
    T exp_val = pseudo_voigt_adl::adl_exp(exp_arg);
    return lambda * T(0.5) * exp_val * erfc_val;
}

template<typename T>
inline T gaussexp_cdf(const T x, const T mean, const T sigma, const T tau) {
    const T sqrt2 = T(1.41421356237309504880);
    if (tau <= T(0)) {
        T diff = (x - mean) / (sigma * sqrt2);
        return T(0.5) * (T(1) + pseudo_voigt_adl::adl_erf(diff));
    }
    const T lambda = T(1) / tau;
    const T lambda_sigma = lambda * sigma;

    T a = (mean - x) / sigma;
    T phi_a = T(0.5) * (T(1) + pseudo_voigt_adl::adl_erf(a / sqrt2));
    T phi_shift = T(0.5) * (T(1) + pseudo_voigt_adl::adl_erf((a - lambda_sigma) / sqrt2));
    T exp_arg = -lambda * (mean - x) + T(0.5) * lambda_sigma * lambda_sigma;

    if (exp_arg > T(700)) return T(1);
    if (exp_arg < T(-700)) {
        T res = T(1) - phi_a;
        if (res < T(0)) res = T(0);
        return res;
    }
    T exp_term = pseudo_voigt_adl::adl_exp(exp_arg);
    T result = T(1) - phi_a + exp_term * phi_shift;
    if (result < T(0)) result = T(0);
    if (result > T(1)) result = T(1);
    return result;
}

// ===============================================================
// Voigt + exponential tail API (pseudo-Voigt core)
// ===============================================================

template<typename T>
inline T voigt_exp_indefinite(const T x, const T mean, const T sigma_gauss,
                              const T gamma_lor, const T tail_ratio, const T tail_slope) {
    T voigt_cdf_val = voigt_cdf(x, mean, sigma_gauss, gamma_lor);
    T gaussexp_cdf_val = gaussexp_cdf(x, mean, sigma_gauss, tail_slope);
    return (T(1) - tail_ratio) * voigt_cdf_val + tail_ratio * gaussexp_cdf_val;
}

inline double voigt_exp_indefinite(const double x, const double mean, const double sigma_gauss,
                                   const double gamma_lor, const double tail_ratio, const double tail_slope) {
    return voigt_exp_indefinite<double>(x, mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
}

inline double voigt_exp_integral(const double peak_mean, const double sigma_gauss,
                                 const double gamma_lor, const double tail_ratio,
                                 const double tail_slope, const double x0, const double x1) {
    double cdf_x1 = voigt_exp_indefinite(x1, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    double cdf_x0 = voigt_exp_indefinite(x0, peak_mean, sigma_gauss, gamma_lor, tail_ratio, tail_slope);
    return cdf_x1 - cdf_x0;
}

template<typename T>
inline void voigt_exp_integral(const T peak_mean, const T sigma_gauss,
                               const T peak_amplitude, const T gamma_lor,
                               const T tail_ratio, const T tail_slope,
                               const float * const energies, T *channels,
                               const size_t nchannel) {
    // Precompute pseudo-Voigt parameters once for this call
    const bool use_gauss = (gamma_lor <= T(1e-10) * sigma_gauss);
    const bool use_lorentz = (sigma_gauss <= T(1e-10) * gamma_lor);
    const auto p = voigt_pseudo_params(sigma_gauss, gamma_lor);

    // Evaluate CDF at bin edges and take differences, avoiding vector allocation
    // Start with CDF at lower edge of first channel
    T x_lower = T(energies[0]);
    T voigt_cdf_lower;
    if (use_gauss) {
        voigt_cdf_lower = gaussian_cdf(x_lower, peak_mean, sigma_gauss);
    } else if (use_lorentz) {
        voigt_cdf_lower = lorentzian_cdf(x_lower, peak_mean, gamma_lor);
    } else {
        voigt_cdf_lower = p.eta * lorentzian_cdf(x_lower, peak_mean, p.gamma_p) +
                         (T(1) - p.eta) * gaussian_cdf(x_lower, peak_mean, p.sigma_p);
        if (voigt_cdf_lower < T(0)) voigt_cdf_lower = T(0);
        if (voigt_cdf_lower > T(1)) voigt_cdf_lower = T(1);
    }
    T gaussexp_cdf_lower = gaussexp_cdf(x_lower, peak_mean, sigma_gauss, tail_slope);
    T prev_cdf = (T(1) - tail_ratio) * voigt_cdf_lower + tail_ratio * gaussexp_cdf_lower;

    // Process each channel
    for (size_t i = 0; i < nchannel; ++i) {
        // Compute CDF at upper edge of this channel
        T x_upper = T(energies[i + 1]);
        T voigt_cdf_upper;
        if (use_gauss) {
            voigt_cdf_upper = gaussian_cdf(x_upper, peak_mean, sigma_gauss);
        } else if (use_lorentz) {
            voigt_cdf_upper = lorentzian_cdf(x_upper, peak_mean, gamma_lor);
        } else {
            voigt_cdf_upper = p.eta * lorentzian_cdf(x_upper, peak_mean, p.gamma_p) +
                             (T(1) - p.eta) * gaussian_cdf(x_upper, peak_mean, p.sigma_p);
            if (voigt_cdf_upper < T(0)) voigt_cdf_upper = T(0);
            if (voigt_cdf_upper > T(1)) voigt_cdf_upper = T(1);
        }
        T gaussexp_cdf_upper = gaussexp_cdf(x_upper, peak_mean, sigma_gauss, tail_slope);
        T current_cdf = (T(1) - tail_ratio) * voigt_cdf_upper + tail_ratio * gaussexp_cdf_upper;

        // Compute bin content from CDF difference
        T bin_content = (current_cdf - prev_cdf) * peak_amplitude;
        if (bin_content < T(0)) bin_content = T(0);
        channels[i] += bin_content;

        // Update previous CDF for next iteration
        prev_cdf = current_cdf;
    }
}

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
            lorentz_width = (2.0 * gamma_lor) / (3.14159265358979323846264338327950288 * p); // asymptotic cot(pi*p/2)
        } else if (p > 0.0) {
            lorentz_width = gamma_lor / std::tan(0.5 * 3.14159265358979323846264338327950288 * p);
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

} // namespace pseudo_voigt

