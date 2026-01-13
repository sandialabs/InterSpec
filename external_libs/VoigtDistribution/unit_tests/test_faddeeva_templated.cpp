/* Test file for templated Faddeeva implementation
 * Compares templated double version against original MIT Faddeeva implementation
 */

#include <boost/test/unit_test.hpp>

#include "../Faddeeva.hpp"
#include "3rd_party/Faddeeva_MIT/Faddeeva.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <cassert>
#include <random>
#include <limits>

#if __has_include("ceres/jet.h") && __has_include("eigen3/Eigen/Core")
#  define HAS_CERES_JET 1
#  include "eigen3/Eigen/Core"
#  include "ceres/jet.h"
#else
#  define HAS_CERES_JET 0
#endif

using namespace FaddeevaT;
using namespace Faddeeva;

namespace {

double rel_err(std::complex<double> a, std::complex<double> b) {
    double mag = std::abs(a);
    double dr = std::abs(a.real() - b.real());
    double di = std::abs(a.imag() - b.imag());
    double num = std::max(dr, di);
    return (mag > 1e-14) ? num / mag : num;
}

void compare_case(const std::complex<double>& z, double relerr,
                  const std::string& label,
                  int& passed, int& failed, double& max_rel_error) {
    std::complex<double> ref = Faddeeva::w(z, relerr);
    double mag = std::abs(ref);
    if (!std::isfinite(ref.real()) || !std::isfinite(ref.imag()) || mag > 1e150) {
        return; // skip regions where reference over/underflows
    }
    Complex<double> zt(z.real(), z.imag());
    auto wt = FaddeevaT::w(zt, relerr);
    std::complex<double> got(wt.real, wt.imag);

    double err = rel_err(ref, got);
    max_rel_error = std::max(max_rel_error, err);
    constexpr double tol = 1e-11;
    if (err < tol) {
        ++passed;
    } else {
        ++failed;
        std::cout << "FAIL " << label << " z=" << z
                  << " err=" << std::scientific << err
                  << " ref=" << ref << " got=" << got << "\n";
    }
}

void test_real_functions(int& passed, int& failed) {
    std::vector<double> xs = {0.0, 0.5, 1.0, 2.0, 5.0, 10.0, -1.0, -5.0};
    for (double x : xs) {
        double r1 = Faddeeva::erfcx(x);
        double r2 = FaddeevaT::erfcx_real(x);
        double e1 = rel_err({r1,0.0}, {r2,0.0});
        if (e1 >= 1e-11) {
            ++failed;
            std::cout << "FAIL erfcx x=" << x << " err=" << std::scientific << e1 << "\n";
        } else {
            ++passed;
        }
        double w1 = Faddeeva::w_im(x);
        double w2 = FaddeevaT::w_im_real(x);
        double e2 = rel_err({w1,0.0}, {w2,0.0});
        if (e2 >= 1e-11) {
            ++failed;
            std::cout << "FAIL w_im x=" << x << " err=" << std::scientific << e2 << "\n";
        } else {
            ++passed;
        }
    }
}

void test_dense_grid(int& passed, int& failed, double& max_rel_error) {
    std::vector<std::complex<double>> zs = {
        {1e-8, 0.0}, {1e-6, 1e-6}, {1e-4, 1e-4},
        {5e-4, 0.0}, {0.1, 0.05}, {0.5, 0.01},
        {1.0, 0.1}, {2.0, 0.5}, {5.0, 0.05},
        {6.0, 0.05}, {7.0, 0.05}, {8.0, 0.2},
        {10.0, 0.1}, {20.0, 0.1}, {6.0, 7.0},
        {8.0, 7.0}, {28.0, 0.1}, {30.0, 10.0},
        {50.0, 0.0}, {0.0, 50.0}
    };
    for (auto z : zs) compare_case(z, 0.0, "grid", passed, failed, max_rel_error);
}

void test_gamma_relevant(int& passed, int& failed, double& max_rel_error) {
    std::vector<double> ys = {0.02, 0.1, 0.5};
    for (double y : ys) {
        for (int i = -8; i <= 8; ++i) {
            double x = static_cast<double>(i);
            compare_case({x, y}, 0.0, "gamma", passed, failed, max_rel_error);
        }
    }
}

void test_fuzz(int& passed, int& failed, double& max_rel_error) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist_x(-12.0, 12.0);
    std::uniform_real_distribution<double> dist_y(0.0, 10.0);
    for (int i = 0; i < 200; ++i) {
        std::complex<double> z(dist_x(rng), dist_y(rng));
        compare_case(z, 0.0, "fuzz", passed, failed, max_rel_error);
    }
}

void test_extremes(int& passed, int& failed, double& max_rel_error) {
    std::vector<std::complex<double>> zs = {
        {1e3, 1e-3}, {1e4, 1e-6}, {1e3, 50.0},
        {1e5, 0.5}, {1e2, 1e-8}
    };
    for (auto z : zs) compare_case(z, 0.0, "extreme", passed, failed, max_rel_error);
}

void test_derivative_jet(int& passed, int& failed) {
#if !HAS_CERES_JET
    std::cout << "Skipping Jet derivative tests (ceres/jet.h not available).\n";
    return;
#else
std::cout << "Running Jet derivative tests.\n";
    const double h = 1e-6;

    // Test cases relevant to gamma spectroscopy Voigt distributions
    // z = ((x - μ) + i*γ)/(σ*sqrt(2))
    // Typical ranges: Re(z) ∈ [-5, 5], Im(z) ∈ [10^-3, 10^-1] for realistic γ/σ ratios
    // Skip very small Im(z) values that have numerical derivative issues
    std::vector<std::complex<double>> zs = {
        // Gamma spectroscopy realistic cases (moderate Im(z), various Re(z))
        {0.0, 1e-3},   // Peak center, moderate Lorentzian
        {1.0, 1e-3},   // Near peak
        {2.0, 5e-4},   // Shoulder region
        {-1.5, 2e-3},  // Left side
        {3.0, 1e-3},   // Right tail
        {-2.5, 5e-3},  // Far left
        {4.0, 2e-3},   // Far right

        // Edge cases with larger Lorentzian components
        {0.5, 0.01},   // Moderate Lorentzian
        {1.0, 0.05},   // Larger Lorentzian
        {2.0, 0.1},    // Significant Lorentzian broadening
        {-1.0, 0.02},  // Negative real, moderate imaginary

        // Original test cases for completeness
        {0.5, 0.1}, {1.0, 0.5}, {5.0, 0.05}, {2.0, 2.0}
    };

    for (auto z : zs) {
        // Test derivative w.r.t. real part
        {
            ceres::Jet<double,1> xr, yr;
            xr.a = z.real(); xr.v[0] = 1.0;
            yr.a = z.imag(); yr.v[0] = 0.0;
            Complex<ceres::Jet<double,1>> zjet(xr, yr);
            auto wj = FaddeevaT::w(zjet, 0.0);

            // finite difference along real
            std::complex<double> fwd = Faddeeva::w({z.real()+h, z.imag()}, 0.0);
            std::complex<double> bwd = Faddeeva::w({z.real()-h, z.imag()}, 0.0);
            std::complex<double> dfdx = (fwd - bwd) / (2*h);
            double err_r = std::abs(wj.real.v[0] - dfdx.real());
            double err_i = std::abs(wj.imag.v[0] - dfdx.imag());
            // Use adaptive tolerance based on magnitude
            double magnitude = std::abs(dfdx);
            double tol = std::max(1e-6, magnitude * 1e-4); // 1e-6 absolute or 0.01% relative
            if (err_r < tol && err_i < tol) {
                passed += 2;
            } else {
                failed += 2;
                std::cout << "FAIL deriv_real z=" << z
                          << " dfdx_fd=" << dfdx
                          << " jet=(" << wj.real.v[0] << "," << wj.imag.v[0] << ")"
                          << " err_r=" << err_r << " err_i=" << err_i << " tol=" << tol << "\n";
            }
        }

        // Test derivative w.r.t. imaginary part (also important for gamma spectroscopy)
        {
            ceres::Jet<double,1> xr, yr;
            xr.a = z.real(); xr.v[0] = 0.0;
            yr.a = z.imag(); yr.v[0] = 1.0;
            Complex<ceres::Jet<double,1>> zjet(xr, yr);
            auto wj = FaddeevaT::w(zjet, 0.0);

            // finite difference along imaginary
            std::complex<double> fwd = Faddeeva::w({z.real(), z.imag()+h}, 0.0);
            std::complex<double> bwd = Faddeeva::w({z.real(), z.imag()-h}, 0.0);
            std::complex<double> dfdy = (fwd - bwd) / (2*h);
            double err_r = std::abs(wj.real.v[0] - dfdy.real());
            double err_i = std::abs(wj.imag.v[0] - dfdy.imag());
            // Use adaptive tolerance based on magnitude
            double magnitude = std::abs(dfdy);
            double tol = std::max(1e-6, magnitude * 1e-4); // 1e-6 absolute or 0.01% relative
            if (err_r < tol && err_i < tol) {
                passed += 2;
            } else {
                failed += 2;
                std::cout << "FAIL deriv_imag z=" << z
                          << " dfdy_fd=" << dfdy
                          << " jet=(" << wj.real.v[0] << "," << wj.imag.v[0] << ")"
                          << " err_r=" << err_r << " err_i=" << err_i << " tol=" << tol << "\n";
            }
        }
    }
#endif
}

void test_w_mit_list(int& passed, int& failed, double& max_rel_error) {
    // z list from MIT TEST_FADDEEVA (57 entries)
    const std::complex<double> zs[] = {
        {624.2,-0.26123},{-0.4,3.0},{0.6,2.0},{-1.0,1.0},{-1.0,-9.0},{-1.0,9.0},
        {-2.34545e-8,1.1234},{-3.0,5.1},{-53.0,30.1},{0.0,0.12345},{11.0,1.0},{-22.0,-2.0},
        {9.0,-28.0},{21.0,-33.0},{1e5,1e5},{1e14,1e14},{-3001.0,-1000.0},{1e160,-1e159},
        {-6.01,0.01},{-0.7,-0.7},{26.1178,4540.909610972489},
        {0.8e7,0.3e7},{-20.0,-19.8081},{1e-16,-1.1e-16},{2.3e-8,1.3e-8},
        {6.3,-1e-13},{6.3,1e-20},{1e-20,6.3},{1e-20,16.3},{9.0,1e-300},
        {6.01,0.11},{8.01,1.01e-10},{28.01,1e-300},{10.01,1e-200},{10.01,-1e-200},
        {10.01,0.99e-10},{10.01,-0.99e-10},{1e-20,7.01},{-1.0,7.01},{5.99,7.01},
        {1.0,0.0},{55.0,0.0},{-0.1,0.0},{1e-20,0.0},{0.0,5e-14},{0.0,51.0},
        {std::numeric_limits<double>::infinity(),0.0},
        {-std::numeric_limits<double>::infinity(),0.0},
        {0.0,std::numeric_limits<double>::infinity()},
        {0.0,-std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()},
        {std::numeric_limits<double>::quiet_NaN(),0.0},
        {0.0,std::numeric_limits<double>::quiet_NaN()},
        {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}
    };
    for (auto z : zs) compare_case(z, 0.0, "mit_w", passed, failed, max_rel_error);
}

void test_erf_mit_list(int& passed, int& failed, double& max_rel_error) {
    // z list from MIT TEST_FADDEEVA erf block (41 entries)
    const std::complex<double> zs[] = {
        {1,2},{-1,2},{1,-2},{-1,-2},{9,-28},{21,-33},{1e3,1e3},{-3001,-1000},{1e160,-1e159},
        {5.1e-3,1e-8},{-4.9e-3,4.95e-3},{4.9e-3,0.5},{4.9e-4,-5.0},{-4.9e-5,-50.0},
        {5.1e-3,0.5},{5.1e-4,-5.0},{-5.1e-5,-50.0},{1e-6,2e-6},{0,2e-6},
        {0,2},{0,20},{0,200},
        {std::numeric_limits<double>::infinity(),0.0},
        {-std::numeric_limits<double>::infinity(),0.0},
        {0.0,std::numeric_limits<double>::infinity()},
        {0.0,-std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()},
        {std::numeric_limits<double>::quiet_NaN(),0.0},
        {0.0,std::numeric_limits<double>::quiet_NaN()},
        {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()},
        {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()},
        {1e-3,std::numeric_limits<double>::quiet_NaN()},
        {7e-2,7e-2},{7e-2,-7e-4},{-9e-2,7e-4},{-9e-2,9e-2},{-7e-4,9e-2},{7e-2,0.9e-2},{7e-2,1.1e-2}
    };
    for (auto z : zs) {
        std::complex<double> ref = Faddeeva::erf(z, 0.0);
        double mag = std::abs(ref);
        if (!std::isfinite(ref.real()) || !std::isfinite(ref.imag()) || mag > 1e150) continue;
        Complex<double> zt(z.real(), z.imag());
        auto wt = FaddeevaT::erf(zt, 0.0);
        std::complex<double> got(wt.real, wt.imag);
        double err = rel_err(ref, got);
        max_rel_error = std::max(max_rel_error, err);
        constexpr double tol = 1e-11;
        if (err < tol) ++passed;
        else {
            ++failed;
            std::cout << "FAIL mit_erf z=" << z
                      << " err=" << std::scientific << err
                      << " ref=" << ref << " got=" << got << "\n";
        }
    }
}

} // namespace

BOOST_AUTO_TEST_CASE(FaddeevaTemplatedImplementationTests) {
    int passed = 0, failed = 0;
    double max_rel_error = 0.0;

    // Basic sanity points
    std::vector<std::complex<double>> basic = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0},
        {2.0, 0.5}, {0.5, 2.0}, {10.0, 0.0}, {0.0, 10.0}
    };
    for (auto z : basic) compare_case(z, 0.0, "basic", passed, failed, max_rel_error);

    test_dense_grid(passed, failed, max_rel_error);
    test_gamma_relevant(passed, failed, max_rel_error);
    test_fuzz(passed, failed, max_rel_error);
    test_extremes(passed, failed, max_rel_error);
    test_w_mit_list(passed, failed, max_rel_error);
    test_erf_mit_list(passed, failed, max_rel_error);
    test_real_functions(passed, failed);
    test_derivative_jet(passed, failed);

    // Check results - should have no failures
    BOOST_CHECK_EQUAL(failed, 0);
}
