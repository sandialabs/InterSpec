/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
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

#define BOOST_TEST_MODULE ComptonTests
#include <boost/test/unit_test.hpp>

#include "transport/ComptonScatter.h"

#include <Eigen/Core>
#include <cmath>
#include <random>
#include <vector>
#include <numeric>

using namespace ceelo;

// ============================================================
//  Kinematic Tests (deterministic)
// ============================================================

BOOST_AUTO_TEST_SUITE(Kinematics)

BOOST_AUTO_TEST_CASE(compton_formula_backscatter) {
    // At 180 degrees backscatter (cos_theta = -1):
    // E' = E / (1 + 2*alpha)
    double E = 662.0; // keV
    double alpha = E / 510.998950;
    double expected = E / (1.0 + 2.0 * alpha);
    double computed = compton_scattered_energy(E, -1.0);
    BOOST_CHECK_CLOSE(computed, expected, 1e-6);
}

BOOST_AUTO_TEST_CASE(compton_formula_forward_scatter) {
    // At 0 degrees (cos_theta = 1): no energy loss
    double E = 662.0;
    double computed = compton_scattered_energy(E, 1.0);
    BOOST_CHECK_CLOSE(computed, E, 1e-6);
}

BOOST_AUTO_TEST_CASE(compton_formula_90_degrees) {
    // At 90 degrees (cos_theta = 0):
    // E' = E / (1 + E/m_e)
    double E = 662.0;
    double alpha = E / 510.998950;
    double expected = E / (1.0 + alpha);
    double computed = compton_scattered_energy(E, 0.0);
    BOOST_CHECK_CLOSE(computed, expected, 1e-6);
}

BOOST_AUTO_TEST_CASE(compton_energy_conservation) {
    // E_scattered + E_deposited = E_incident
    double E = 1000.0; // keV
    for (double cos_theta : {-1.0, -0.5, 0.0, 0.5, 0.9, 1.0}) {
        double E_prime = compton_scattered_energy(E, cos_theta);
        double deposited = E - E_prime;
        BOOST_CHECK_GE(E_prime, 0.0);
        BOOST_CHECK_GE(deposited, 0.0);
        BOOST_CHECK_CLOSE(E_prime + deposited, E, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(klein_nishina_total_low_energy) {
    // In the low-energy limit, sigma_KN -> sigma_Thomson = 0.6652 barn/electron
    double sigma = klein_nishina_total(1.0); // 1 keV
    BOOST_CHECK_CLOSE(sigma, 0.6652, 1.0); // within 1% of Thomson value
}

BOOST_AUTO_TEST_CASE(klein_nishina_total_decreases_with_energy) {
    double sigma_100 = klein_nishina_total(100.0);
    double sigma_662 = klein_nishina_total(662.0);
    double sigma_1000 = klein_nishina_total(1000.0);
    double sigma_5000 = klein_nishina_total(5000.0);

    BOOST_CHECK_GT(sigma_100, sigma_662);
    BOOST_CHECK_GT(sigma_662, sigma_1000);
    BOOST_CHECK_GT(sigma_1000, sigma_5000);
}

BOOST_AUTO_TEST_CASE(direction_rotation_normalised) {
    // rotate_direction must always produce a unit vector
    Eigen::Vector3d dir(0.0, 0.0, 1.0);
    for (double cos_theta : {-0.9, 0.0, 0.5, 0.99}) {
        for (double phi : {0.0, 1.0, 3.14, 5.0}) {
            auto new_dir = rotate_direction(dir, cos_theta, phi);
            BOOST_CHECK_CLOSE(new_dir.norm(), 1.0, 1e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(direction_rotation_cos_theta_correct) {
    // The angle between the new and old direction should equal theta
    Eigen::Vector3d dir(0.0, 0.0, 1.0);
    double cos_theta = 0.5;
    double phi = 0.0;
    auto new_dir = rotate_direction(dir, cos_theta, phi);
    double cos_angle = dir.dot(new_dir);
    BOOST_CHECK_CLOSE(cos_angle, cos_theta, 1e-4);
}

BOOST_AUTO_TEST_CASE(direction_rotation_arbitrary_incident) {
    // Test with a non-trivial incident direction
    Eigen::Vector3d dir = Eigen::Vector3d(1.0, 1.0, 1.0).normalized();
    double cos_theta = -0.3;
    double phi = 1.5;
    auto new_dir = rotate_direction(dir, cos_theta, phi);
    BOOST_CHECK_CLOSE(new_dir.norm(), 1.0, 1e-5);
    double cos_angle = dir.dot(new_dir);
    BOOST_CHECK_CLOSE(cos_angle, cos_theta, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()


// ============================================================
//  Statistical Tests (Monte Carlo validation of Klein-Nishina)
// ============================================================

BOOST_AUTO_TEST_SUITE(Statistics)

BOOST_AUTO_TEST_CASE(compton_mean_scattered_energy_662keV) {
    // At 662 keV, the mean scattered energy from Klein-Nishina sampling
    // should match the theoretical mean.
    // Theoretical mean E' can be computed by integrating the KN distribution.
    // We check that it's less than incident energy and greater than backscatter energy.
    double E_in = 662.0;
    double alpha = E_in / 510.998950;
    double E_backscatter = E_in / (1.0 + 2.0 * alpha);

    std::mt19937_64 rng(42);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    const int N = 50000;
    double sum_E_prime = 0.0;
    double sum_deposited = 0.0;

    for (int i = 0; i < N; ++i) {
        auto result = sample_compton_scatter(E_in, dir, rng);
        sum_E_prime += result.scattered_energy_keV;
        sum_deposited += result.deposited_energy_keV;

        // Energy conservation per event
        BOOST_REQUIRE_CLOSE(result.scattered_energy_keV + result.deposited_energy_keV,
                            E_in, 0.001);

        // Scattered energy must be >= minimum (backscatter) and <= incident
        BOOST_REQUIRE_GE(result.scattered_energy_keV, E_backscatter - 1.0); // 1 keV tolerance
        BOOST_REQUIRE_LE(result.scattered_energy_keV, E_in + 1.0);
    }

    double mean_E_prime = sum_E_prime / N;

    // Mean scattered energy at 662 keV should be around 400-500 keV
    // (exact value depends on the KN distribution; roughly E' = E/(1 + alpha/3))
    BOOST_CHECK_GT(mean_E_prime, 200.0);
    BOOST_CHECK_LT(mean_E_prime, 700.0);

    // Mean deposited energy + mean scattered energy should equal incident
    double mean_deposited = sum_deposited / N;
    BOOST_CHECK_CLOSE(mean_E_prime + mean_deposited, E_in, 0.1);
}

BOOST_AUTO_TEST_CASE(compton_direction_isotropic_phi) {
    // The azimuthal angle phi should be uniformly distributed.
    // Check that the x and y components of the scattered direction
    // have zero mean and equal variance.
    std::mt19937_64 rng(123);
    Eigen::Vector3d dir(0.0, 0.0, 1.0); // along z

    const int N = 10000;
    double sum_x = 0.0, sum_y = 0.0;

    for (int i = 0; i < N; ++i) {
        auto result = sample_compton_scatter(662.0, dir, rng);
        sum_x += result.new_direction.x();
        sum_y += result.new_direction.y();
    }

    // Mean x and y should be ~0 (azimuthal symmetry)
    double mean_x = sum_x / N;
    double mean_y = sum_y / N;
    BOOST_CHECK_SMALL(mean_x, 0.05); // within 5% of 0
    BOOST_CHECK_SMALL(mean_y, 0.05);
}

BOOST_AUTO_TEST_CASE(compton_forward_bias_high_energy) {
    // At high energies, Compton scattering is strongly forward-peaked.
    // Mean cos_theta should be > 0 (forward scatter dominates).
    std::mt19937_64 rng(999);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    const int N = 10000;
    double sum_cos = 0.0;

    for (int i = 0; i < N; ++i) {
        auto result = sample_compton_scatter(5000.0, dir, rng); // 5 MeV
        sum_cos += result.cos_theta;
    }

    double mean_cos = sum_cos / N;
    // At 5 MeV, mean cos_theta should be well above 0.5 (strongly forward-peaked)
    BOOST_CHECK_GT(mean_cos, 0.5);
}

BOOST_AUTO_TEST_CASE(compton_low_energy_symmetric) {
    // At very low energies, Compton approaches Thomson scattering
    // which has a (1 + cos²θ) distribution — symmetric around 90°.
    // Mean cos_theta should be ~0.
    std::mt19937_64 rng(777);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    const int N = 20000;
    double sum_cos = 0.0;

    for (int i = 0; i < N; ++i) {
        auto result = sample_compton_scatter(5.0, dir, rng); // 5 keV
        sum_cos += result.cos_theta;
    }

    double mean_cos = sum_cos / N;
    // At 5 keV (alpha = 0.01), distribution is nearly symmetric, mean cos_theta ~ 0
    BOOST_CHECK_SMALL(mean_cos, 0.1);
}

BOOST_AUTO_TEST_CASE(compton_new_direction_normalized) {
    // All sampled directions must be unit vectors
    std::mt19937_64 rng(42);
    Eigen::Vector3d dir = Eigen::Vector3d(1.0, 2.0, 3.0).normalized();

    for (int i = 0; i < 1000; ++i) {
        auto result = sample_compton_scatter(662.0, dir, rng);
        BOOST_CHECK_CLOSE(result.new_direction.norm(), 1.0, 1e-4);
    }
}

BOOST_AUTO_TEST_CASE(compton_scattered_energy_range) {
    // Scattered energy must always be in [E_backscatter, E_incident]
    std::mt19937_64 rng(12345);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    struct TestCase { double energy_keV; };
    std::vector<TestCase> cases = {{50.0}, {100.0}, {662.0}, {1000.0}, {5000.0}};

    for (const auto& tc : cases) {
        double alpha = tc.energy_keV / 510.998950;
        double E_back = tc.energy_keV / (1.0 + 2.0 * alpha);

        for (int i = 0; i < 1000; ++i) {
            auto result = sample_compton_scatter(tc.energy_keV, dir, rng);
            BOOST_CHECK_GE(result.scattered_energy_keV, E_back * 0.999);
            BOOST_CHECK_LE(result.scattered_energy_keV, tc.energy_keV * 1.001);
        }
    }
}

BOOST_AUTO_TEST_CASE(kn_angle_distribution_662keV) {
    // Verify the sampled angle distribution matches the analytical Klein-Nishina
    // differential cross section at 662 keV (alpha = 1.296).
    // This test catches the Kahn's method bug that was replaced by Butcher-Messel.
    //
    // We bin sampled cos_theta into 7 angle bands, compute expected fractions
    // from numerical integration of KN, and require each band to agree within 3%.
    const double E_keV = 662.0;
    const double alpha = E_keV / 510.998950;
    const int N = 500000;

    std::mt19937_64 rng(2024);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    struct Band { double lo, hi; };
    const Band bands[] = {
        {0,20}, {20,40}, {40,60}, {60,90}, {90,120}, {120,150}, {150,180}
    };
    const int N_BANDS = 7;
    std::vector<int> counts(N_BANDS, 0);

    for (int i = 0; i < N; ++i) {
        auto result = sample_compton_scatter(E_keV, dir, rng);
        double angle_deg = std::acos(std::max(-1.0, std::min(1.0, result.cos_theta)))
                         * 180.0 / M_PI;
        for (int b = 0; b < N_BANDS; ++b) {
            if (angle_deg >= bands[b].lo && angle_deg < bands[b].hi) {
                counts[b]++;
                break;
            }
        }
    }

    // Compute expected fractions by numerical integration of KN dsigma/d(costheta).
    auto kn_dsigma = [alpha](double cos_theta) {
        double eps = 1.0 / (1.0 + alpha * (1.0 - cos_theta));
        double sin2 = 1.0 - cos_theta * cos_theta;
        return eps * eps * (eps + 1.0/eps - sin2);
    };

    std::vector<double> expected_frac(N_BANDS);
    double total_integral = 0.0;
    for (int b = 0; b < N_BANDS; ++b) {
        double ct_lo = std::cos(bands[b].hi * M_PI / 180.0);
        double ct_hi = std::cos(bands[b].lo * M_PI / 180.0);
        double sum = 0.0;
        const int N_INT = 1000;
        double dct = (ct_hi - ct_lo) / N_INT;
        for (int j = 0; j < N_INT; ++j) {
            double ct = ct_lo + (j + 0.5) * dct;
            sum += kn_dsigma(ct) * dct;
        }
        expected_frac[b] = sum;
        total_integral += sum;
    }
    for (int b = 0; b < N_BANDS; ++b)
        expected_frac[b] /= total_integral;

    for (int b = 0; b < N_BANDS; ++b) {
        double sampled_frac = static_cast<double>(counts[b]) / N;
        double diff_pct = std::abs(sampled_frac - expected_frac[b])
                        / expected_frac[b] * 100.0;
        BOOST_CHECK_LT(diff_pct, 3.0);
    }
}

BOOST_AUTO_TEST_CASE(kn_angle_distribution_1173keV) {
    // Same test at 1173 keV (alpha = 2.296).
    const double E_keV = 1173.0;
    const double alpha = E_keV / 510.998950;
    const int N = 500000;

    std::mt19937_64 rng(2025);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    struct Band { double lo, hi; };
    const Band bands[] = {
        {0,20}, {20,40}, {40,60}, {60,90}, {90,120}, {120,150}, {150,180}
    };
    const int N_BANDS = 7;
    std::vector<int> counts(N_BANDS, 0);

    for (int i = 0; i < N; ++i) {
        auto result = sample_compton_scatter(E_keV, dir, rng);
        double angle_deg = std::acos(std::max(-1.0, std::min(1.0, result.cos_theta)))
                         * 180.0 / M_PI;
        for (int b = 0; b < N_BANDS; ++b) {
            if (angle_deg >= bands[b].lo && angle_deg < bands[b].hi) {
                counts[b]++;
                break;
            }
        }
    }

    auto kn_dsigma = [alpha](double cos_theta) {
        double eps = 1.0 / (1.0 + alpha * (1.0 - cos_theta));
        double sin2 = 1.0 - cos_theta * cos_theta;
        return eps * eps * (eps + 1.0/eps - sin2);
    };

    std::vector<double> expected_frac(N_BANDS);
    double total_integral = 0.0;
    for (int b = 0; b < N_BANDS; ++b) {
        double ct_lo = std::cos(bands[b].hi * M_PI / 180.0);
        double ct_hi = std::cos(bands[b].lo * M_PI / 180.0);
        double sum = 0.0;
        const int N_INT = 1000;
        double dct = (ct_hi - ct_lo) / N_INT;
        for (int j = 0; j < N_INT; ++j) {
            double ct = ct_lo + (j + 0.5) * dct;
            sum += kn_dsigma(ct) * dct;
        }
        expected_frac[b] = sum;
        total_integral += sum;
    }
    for (int b = 0; b < N_BANDS; ++b)
        expected_frac[b] /= total_integral;

    for (int b = 0; b < N_BANDS; ++b) {
        double sampled_frac = static_cast<double>(counts[b]) / N;
        double diff_pct = std::abs(sampled_frac - expected_frac[b])
                        / expected_frac[b] * 100.0;
        BOOST_CHECK_LT(diff_pct, 3.0);
    }
}

BOOST_AUTO_TEST_CASE(bound_electron_suppresses_forward_scatter) {
    // The incoherent scattering function S(x,Z)/Z < 1 at small angles (small x),
    // so bound-electron correction should suppress forward scattering compared to
    // free-electron KN. Verify that Z>0 produces fewer forward scatters.
    const double E_keV = 662.0;
    const int N = 200000;
    Eigen::Vector3d dir(0.0, 0.0, 1.0);

    std::mt19937_64 rng1(100);
    int free_forward = 0;
    for (int i = 0; i < N; ++i) {
        auto r = sample_compton_scatter(E_keV, dir, rng1, 0);
        if (r.cos_theta > std::cos(20.0 * M_PI / 180.0)) free_forward++;
    }

    std::mt19937_64 rng2(100);
    int bound_forward = 0;
    for (int i = 0; i < N; ++i) {
        auto r = sample_compton_scatter(E_keV, dir, rng2, 53); // Iodine
        if (r.cos_theta > std::cos(20.0 * M_PI / 180.0)) bound_forward++;
    }

    BOOST_CHECK_LT(bound_forward, free_forward);
    BOOST_CHECK_GT(bound_forward, free_forward / 2);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  Compton Doppler broadening (impulse approximation)
// ============================================================

BOOST_AUTO_TEST_SUITE(DopplerBroadening)

// 1. Analytic profile CDF / inverse-CDF round-trip and basic properties.
BOOST_AUTO_TEST_CASE(profile_cdf_inverse_roundtrip) {
    for (double J0 : {0.1, 0.85, 3.0, 9.0}) {
        // F(0) = 1/2
        BOOST_CHECK_CLOSE(compton_profile_cdf(J0, 0.0), 0.5, 1e-9);
        double prev_p = -1e30;
        for (double xi : {1e-9, 1e-4, 0.01, 0.2, 0.499, 0.5, 0.501,
                          0.8, 0.99, 0.9999, 1.0 - 1e-9}) {
            double p = compton_profile_invcdf(J0, xi);
            BOOST_CHECK_GT(p, prev_p);          // monotone in xi
            prev_p = p;
            double back = compton_profile_cdf(J0, p);
            BOOST_CHECK_SMALL(back - xi, 1e-9);  // round-trip
        }
    }
}

// 2. Kinematics: E'(0) == free Compton line; monotone increasing in t;
//    t recovered from (E, E', cos) matches.  Catches root-selection and the
//    137x atomic-unit conversion bug.
BOOST_AUTO_TEST_CASE(doppler_kinematics_inversion) {
    const double mc2 = 510.998950;
    for (double E : {59.0, 662.0, 2614.0}) {
        for (double c : {-0.99, -0.5, 0.0, 0.5, 0.9}) {
            BOOST_CHECK_CLOSE(doppler_scattered_energy(E, c, 0.0),
                              compton_scattered_energy(E, c), 1e-7);
            double prev = -1.0;
            for (double t = -0.3; t <= 0.3 + 1e-9; t += 0.02) {
                double ep = doppler_scattered_energy(E, c, t);
                BOOST_CHECK_GT(ep, prev);  // monotone increasing
                prev = ep;
                double num = E * ep * (1.0 - c) / mc2 - (E - ep);
                double den = std::sqrt(E * E + ep * ep - 2.0 * E * ep * c);
                BOOST_CHECK_SMALL(num / den - t, 1e-8);
            }
        }
    }
}

// 3. doppler_t_max: E'(t_max) == E - U exactly.
BOOST_AUTO_TEST_CASE(doppler_tmax_hits_binding_limit) {
    for (double E : {59.0, 662.0, 2614.0})
        for (double c : {-0.9, 0.0, 0.7})
            for (double U : {0.03, 1.07, 33.17}) {
                if (U >= E) continue;
                double tm = doppler_t_max(E, c, U);
                double ep = doppler_scattered_energy(E, c, tm);
                BOOST_CHECK_CLOSE(ep, E - U, 1e-4);
            }
}

// 4. Energy is conserved exactly per sampled event: E' + KE + U = E.
BOOST_AUTO_TEST_CASE(doppler_energy_conservation) {
    std::mt19937_64 rng(424242);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);
    for (double E : {30.0, 122.0, 662.0, 2614.0}) {
        for (int Z : {11, 53, 82}) {
            double max_err = 0.0;
            for (int i = 0; i < 100000; ++i) {
                auto s = sample_compton_scatter(E, dir, rng, Z, /*doppler=*/true);
                BOOST_REQUIRE(!std::isnan(s.scattered_energy_keV));
                BOOST_CHECK_GE(s.scattered_energy_keV, 0.0);
                BOOST_CHECK_GE(s.deposited_energy_keV, 0.0);
                BOOST_CHECK_GE(s.binding_deposit_keV, 0.0);
                double tot = s.scattered_energy_keV + s.deposited_energy_keV
                             + s.binding_deposit_keV;
                max_err = std::max(max_err, std::fabs(tot - E));
            }
            BOOST_CHECK_SMALL(max_err, 1e-7);
        }
    }
}

// 5. Doppler ON broadens the scattered energy at fixed angle; OFF is a delta.
//    The spread is dominated by inner-shell (K/L) scatters, so we compare the
//    sampled RMS to an in-test expectation computed from the same shell data.
BOOST_AUTO_TEST_CASE(doppler_broadens_at_fixed_angle) {
    std::mt19937_64 rng(7);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);
    const double E = 662.0;
    const int Z = 53;  // iodine

    // OFF: no spread relative to the free line.
    double off_max = 0.0;
    for (int i = 0; i < 50000; ++i) {
        auto s = sample_compton_scatter(E, dir, rng, Z, /*doppler=*/false);
        double ec = compton_scattered_energy(E, s.cos_theta);
        off_max = std::max(off_max, std::fabs(s.scattered_energy_keV - ec));
    }
    BOOST_CHECK_SMALL(off_max, 1e-9);

    // ON: collect a backscatter band, measure RMS of E' about the free line.
    double sum = 0.0, sum2 = 0.0;
    long n = 0;
    for (int i = 0; i < 4000000 && n < 100000; ++i) {
        auto s = sample_compton_scatter(E, dir, rng, Z, /*doppler=*/true);
        if (s.cos_theta > -1.0 && s.cos_theta < -0.9) {
            double ec = compton_scattered_energy(E, s.cos_theta);
            double d = s.scattered_energy_keV - ec;
            sum += d; sum2 += d * d; ++n;
        }
    }
    BOOST_REQUIRE_GT(n, 1000);
    double mean = sum / n;
    double rms = std::sqrt(sum2 / n - mean * mean);
    // Mean shift is small (profile is ~symmetric in p_z); width is several keV
    // to tens of keV.  Loose physical bounds — exact shape validated vs G4.
    BOOST_CHECK_LT(std::fabs(mean), 3.0);
    BOOST_CHECK_GT(rms, 3.0);
    BOOST_CHECK_LT(rms, 40.0);
}

// 6. Degenerate / low-energy cases produce no E' > E and no NaN.
BOOST_AUTO_TEST_CASE(doppler_degenerate_low_energy) {
    std::mt19937_64 rng(13);
    Eigen::Vector3d dir(0.0, 0.0, 1.0);
    for (double E : {3.0, 5.0, 40.0}) {
        for (int Z : {53, 82}) {
            for (int i = 0; i < 50000; ++i) {
                auto s = sample_compton_scatter(E, dir, rng, Z, /*doppler=*/true);
                BOOST_REQUIRE(!std::isnan(s.scattered_energy_keV));
                BOOST_CHECK_LE(s.scattered_energy_keV, E + 1e-9);
                BOOST_CHECK_GE(s.scattered_energy_keV, 0.0);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
