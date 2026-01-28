/* Copyright (c) 2012 Massachusetts Institute of Technology
 *
 * Upstream code by Steven G. Johnson, MIT (http://ab-initio.mit.edu/faddeeva/),
 * updated 12 May 2015. This file is a templated adaptation to support ceres::Jet<>
 * performed by Will Johnson (wcjohns@sandia.gov) using an LLM-assisted workflow.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/* Templated version of Faddeeva functions for use with double and ceres::Jet<>
 * Based on the MIT Faddeeva implementation
 *
 * This is a C++-only templated version that works with both double and
 * ceres::Jet<> for automatic differentiation.
 */

#ifndef FADDEEVA_TEMPLATED_HPP
#define FADDEEVA_TEMPLATED_HPP 1

#include <cmath>
#include <limits>
#include <complex>
#include <iostream>

// Allow ADL to pick ceres::Jet overloads while still falling back to std::
using std::exp;
using std::sqrt;
using std::sin;
using std::cos;
using std::atan2;
// Note: abs is not included here to allow ADL to find ceres::Jet overloads

// Forward declaration - Complex<T> should be defined in the including file
// or we can define a simple one here
template<typename T>
struct Complex {
    T real;
    T imag;

    Complex() : real(T(0)), imag(T(0)) {}
    Complex(T r) : real(r), imag(T(0)) {}
    Complex(T r, T i) : real(r), imag(i) {}

    Complex<T> operator+(const Complex<T>& other) const {
        return Complex<T>(real + other.real, imag + other.imag);
    }

    Complex<T> operator-(const Complex<T>& other) const {
        return Complex<T>(real - other.real, imag - other.imag);
    }

    Complex<T> operator*(const Complex<T>& other) const {
        return Complex<T>(real * other.real - imag * other.imag,
                         real * other.imag + imag * other.real);
    }

    Complex<T> operator/(const Complex<T>& other) const {
        T denom = other.real * other.real + other.imag * other.imag;
        return Complex<T>((real * other.real + imag * other.imag) / denom,
                         (imag * other.real - real * other.imag) / denom);
    }

    Complex<T> operator*(const T& scalar) const {
        return Complex<T>(real * scalar, imag * scalar);
    }

    Complex<T> operator/(const T& scalar) const {
        return Complex<T>(real / scalar, imag / scalar);
    }

    Complex<T> operator-() const {
        return Complex<T>(-real, -imag);
    }

    Complex<T> conj() const {
        return Complex<T>(real, -imag);
    }

    T abs2() const {
        return real * real + imag * imag;
    }
};

template<typename T>
Complex<T> complex_exp(const Complex<T>& z) {
    T exp_real = exp(z.real);
    return Complex<T>(exp_real * cos(z.imag), exp_real * sin(z.imag));
}

template<typename T>
Complex<T> complex_sqrt(const Complex<T>& z) {
    T r = sqrt(z.abs2());
    T theta = atan2(z.imag, z.real);
    T sqrt_r = sqrt(r);
    return Complex<T>(sqrt_r * cos(theta / 2), sqrt_r * sin(theta / 2));
}

// Helper to get scalar value from T (works for double and ceres::Jet<>)
template<typename T>
double get_scalar(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return x;
    } else {
        // For ceres::Jet<>, return the scalar part (a)
        return x.a;
    }
}

// Helper to check if value is NaN
template<typename T>
bool is_nan(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return std::isnan(x);
    } else {
        return std::isnan(x.a);
    }
}

// Helper to check if value is Inf
template<typename T>
bool is_inf(const T& x) {
    if constexpr (std::is_same_v<T, double>) {
        return std::isinf(x);
    } else {
        return std::isinf(x.a);
    }
}

namespace FaddeevaT {

// Constants
template<typename T>
constexpr T get_inf() {
    if constexpr (std::is_same_v<T, double>) {
        return std::numeric_limits<double>::infinity();
    } else {
        return T(std::numeric_limits<double>::infinity());
    }
}

template<typename T>
constexpr T get_nan() {
    if constexpr (std::is_same_v<T, double>) {
        return std::numeric_limits<double>::quiet_NaN();
    } else {
        return T(std::numeric_limits<double>::quiet_NaN());
    }
}

// Precomputed table of expa2n2[n-1] = exp(-a2*n*n)
// for double-precision a2 = 0.26865... in w(z), below.
static const double expa2n2[] = {
  7.64405281671221563e-01,
  3.41424527166548425e-01,
  8.91072646929412548e-02,
  1.35887299055460086e-02,
  1.21085455253437481e-03,
  6.30452613933449404e-05,
  1.91805156577114683e-06,
  3.40969447714832381e-08,
  3.54175089099469393e-10,
  2.14965079583260682e-12,
  7.62368911833724354e-15,
  1.57982797110681093e-17,
  1.91294189103582677e-20,
  1.35344656764205340e-23,
  5.59535712428588720e-27,
  1.35164257972401769e-30,
  1.90784582843501167e-34,
  1.57351920291442930e-38,
  7.58312432328032845e-43,
  2.13536275438697082e-47,
  3.51352063787195769e-52,
  3.37800830266396920e-57,
  1.89769439468301000e-62,
  6.22929926072668851e-68,
  1.19481172006938722e-73,
  1.33908181133005953e-79,
  8.76924303483223939e-86,
  3.35555576166254986e-92,
  7.50264110688173024e-99,
  9.80192200745410268e-106,
  7.48265412822268959e-113,
  3.33770122566809425e-120,
  8.69934598159861140e-128,
  1.32486951484088852e-135,
  1.17898144201315253e-143,
  6.13039120236180012e-152,
  1.86258785950822098e-160,
  3.30668408201432783e-169,
  3.43017280887946235e-178,
  2.07915397775808219e-187,
  7.36384545323984966e-197,
  1.52394760394085741e-206,
  1.84281935046532100e-216,
  1.30209553802992923e-226,
  5.37588903521080531e-237,
  1.29689584599763145e-247,
  1.82813078022866562e-258,
  1.50576355348684241e-269,
  7.24692320799294194e-281,
  2.03797051314726829e-292,
  3.34880215927873807e-304,
  0.0 // underflow (also prevents reads past array end, below)
};

// Helper functions
template<typename T>
static inline T sqr(const T& x) { return x*x; }

template<typename T>
static inline T sinc(const T& x, const T& sinx) {
    double abs_x = std::abs(get_scalar(x));
    return abs_x < 1e-4 ? T(1) - (T(0.1666666666666666666667))*x*x : sinx / x;
}

template<typename T>
static inline T sinh_taylor(const T& x) {
    T x2 = x*x;
    return x * (T(1) + x2 * (T(0.1666666666666666666667) + T(0.00833333333333333333333) * x2));
}

// Forward declarations
template<typename T>
T erfcx_real(const T& x);

template<typename T>
T w_im_real(const T& x);

template<typename T>
Complex<T> w(const Complex<T>& z, double relerr = 0.0);

template<typename T>
Complex<T> erf(const Complex<T>& z, double relerr = 0.0);

template<typename T>
Complex<T> erfc(const Complex<T>& z, double relerr = 0.0);

template<typename T>
Complex<T> erfcx(const Complex<T>& z, double relerr = 0.0);

template<typename T>
Complex<T> erfi(const Complex<T>& z, double relerr = 0.0);

template<typename T>
Complex<T> Dawson(const Complex<T>& z, double relerr = 0.0);

// Real-argument versions
template<typename T>
T erf_real(const T& x);

template<typename T>
T erfc_real(const T& x);

template<typename T>
T erfcx_real(const T& x);

template<typename T>
T erfi_real(const T& x);

template<typename T>
T Dawson_real(const T& x);

template<typename T>
T w_im_real(const T& x);

} // namespace FaddeevaT

// Include implementation
#include "Faddeeva_impl.hpp"

#endif // FADDEEVA_TEMPLATED_HPP