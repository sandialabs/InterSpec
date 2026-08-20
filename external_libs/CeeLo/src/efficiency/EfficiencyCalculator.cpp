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

#include "efficiency/EfficiencyCalculator.h"
#include "cascade/AngularSampling.h"
#include "cascade/LevelDag.h"
#include "cross_sections/CrossSectionData.h"
#include "export/Geant4Export.h"
#include "geometry/RayTrace.h"
#include "geometry/SourceGeometry.h"
#include "physics/ElectronCsda.h"
#include "transport/ComptonScatter.h"
#include "transport/PhotonTransport.h"
#include "transport/TransportUtils.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <iterator>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <future>
#include <array>
#include <random>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <ctime>
#endif

namespace ceelo {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kFourPi = 4.0 * kPi;
constexpr double kFepTolerance = 1.5; // keV tolerance for full-energy peak

// Below this kinetic energy an internal-conversion / Auger electron is treated
// as depositing locally (or dropped): its CSDA range is sub-mm and it cannot
// bridge any real source->detector gap, so transporting it buys nothing. Also
// the Seltzer-Berger brems tables floor at 200 keV, so no brems is lost below.
constexpr double kIcElectronMin_keV = 15.0;

/// Current efficiency estimates from a merged tally. IS variance:
/// var(eps) = (sum_w_sq/N - eps^2)/N; reduces to binomial eps(1-eps)/N for
/// unit weights. All zero when no weight has been accumulated yet.
struct TallyEstimates {
    double eps_fep = 0.0, sig_fep = 0.0;
    double eps_tot = 0.0, sig_tot = 0.0;
};

/// CPU time consumed by the CALLING thread so far, in seconds. Per-thread (not
/// process-wide) so the sum over worker threads stays correct when several
/// simulations share one process (future node-level parallelism) and under
/// machine load.
double thread_cpu_seconds() {
#if defined(_WIN32)
    FILETIME creation, exit_t, kernel, user;
    if (GetThreadTimes(GetCurrentThread(), &creation, &exit_t, &kernel, &user)) {
        const auto to_sec = [](const FILETIME& ft) {
            const uint64_t ticks = (uint64_t(ft.dwHighDateTime) << 32) | ft.dwLowDateTime;
            return double(ticks) * 1e-7;  // 100 ns ticks
        };
        return to_sec(kernel) + to_sec(user);
    }
    return 0.0;
#else
    timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0)
        return double(ts.tv_sec) + double(ts.tv_nsec) * 1e-9;
    return 0.0;
#endif
}

TallyEstimates compute_estimates(const EfficiencyCalculator::ThreadTally& t) {
    TallyEstimates est;
    if (t.sum_weights <= 0.0 || t.num_events == 0)
        return est;
    const double N = static_cast<double>(t.num_events);
    est.eps_fep = t.sum_fep_weights / t.sum_weights;
    est.sig_fep = std::sqrt(std::max(0.0,
        (t.sum_fep_w_sq / N - est.eps_fep * est.eps_fep) / N));
    est.eps_tot = t.sum_any_weights / t.sum_weights;
    est.sig_tot = std::sqrt(std::max(0.0,
        (t.sum_any_w_sq / N - est.eps_tot * est.eps_tot) / N));
    return est;
}

/// Shared simulation state for batched cooperative execution.
struct SimulationState {
    // Merged tally (protected by mutex)
    EfficiencyCalculator::ThreadTally merged;
    std::mutex mutex;

    // Termination flag (lock-free)
    std::atomic<bool> stop_flag{false};

    // Configuration
    TerminationCriteria termination;
    ProgressCallback progress_callback;
    std::chrono::steady_clock::time_point start_time;
    std::chrono::steady_clock::time_point last_callback_time;

    // Result
    StopReason stop_reason = StopReason::MaxEvents;

    // CPU time summed across worker threads (protected by mutex)
    double cpu_seconds = 0.0;

    void merge_and_check(EfficiencyCalculator::ThreadTally& batch, double batch_cpu_s) {
        std::lock_guard<std::mutex> lock(mutex);

        cpu_seconds += batch_cpu_s;

        merged.num_events += batch.num_events;
        merged.num_fep += batch.num_fep;
        merged.num_any += batch.num_any;
        merged.sum_weights += batch.sum_weights;
        merged.sum_fep_weights += batch.sum_fep_weights;
        merged.sum_any_weights += batch.sum_any_weights;
        merged.sum_fep_w_sq += batch.sum_fep_w_sq;
        merged.sum_any_w_sq += batch.sum_any_w_sq;
        merged.num_forced_absorption += batch.num_forced_absorption;
        for (size_t i = 0; i < merged.pulse_height_counts.size(); ++i) {
            merged.pulse_height_counts[i] += batch.pulse_height_counts[i];
            merged.pulse_height_counts_sq[i] += batch.pulse_height_counts_sq[i];
        }
        // Merge Marinelli diagnostics
        merged.dbg_initial_hit += batch.dbg_initial_hit;
        merged.dbg_miss_bounce_hit += batch.dbg_miss_bounce_hit;
        merged.dbg_reentry_hit += batch.dbg_reentry_hit;
        merged.dbg_secondary_hit += batch.dbg_secondary_hit;
        merged.dbg_pp_secondary += batch.dbg_pp_secondary;
        merged.dbg_dep_initial += batch.dbg_dep_initial;
        merged.dbg_dep_reentry += batch.dbg_dep_reentry;
        merged.dbg_dep_secondary += batch.dbg_dep_secondary;
        merged.dbg_dep_pp += batch.dbg_dep_pp;
        merged.dbg_pp_only_any += batch.dbg_pp_only_any;
        merged.dbg_pp_only_any_w += batch.dbg_pp_only_any_w;
        merged.dbg_n_escaped += batch.dbg_n_escaped;
        merged.dbg_n_can_reenter += batch.dbg_n_can_reenter;
        merged.dbg_n_wall_pass += batch.dbg_n_wall_pass;
        merged.dbg_n_wall_scatter += batch.dbg_n_wall_scatter;
        merged.dbg_n_wall_absorb += batch.dbg_n_wall_absorb;
        merged.dbg_n_water_survived += batch.dbg_n_water_survived;
        merged.dbg_n_trace_hit += batch.dbg_n_trace_hit;
        merged.dbg_n_trace_miss += batch.dbg_n_trace_miss;
        merged.dbg_src.merge(batch.dbg_src);

        // Check termination criteria
        auto now = std::chrono::steady_clock::now();
        auto elapsed = now - start_time;
        double elapsed_sec = std::chrono::duration<double>(elapsed).count();

        // Max events
        if (termination.max_events > 0 && merged.num_events >= termination.max_events) {
            stop_reason = StopReason::MaxEvents;
            stop_flag.store(true, std::memory_order_release);
        }

        // Wall time
        if (termination.max_wall_seconds > 0.0 && elapsed_sec >= termination.max_wall_seconds) {
            stop_reason = StopReason::WallTime;
            stop_flag.store(true, std::memory_order_release);
        }

        // CPU time (summed across worker threads; load-invariant budget cap)
        if (termination.max_cpu_seconds > 0.0 && cpu_seconds >= termination.max_cpu_seconds) {
            stop_reason = StopReason::CpuTime;
            stop_flag.store(true, std::memory_order_release);
        }

        // Current efficiency estimates (IS variance; shared by the precision
        // checks and the progress callback payload).
        const TallyEstimates est = compute_estimates(merged);

        // Precision checks (only after min_events)
        if (merged.num_events >= termination.min_events) {
            if (termination.target_fep_rel_precision > 0.0 && est.eps_fep > 0.0 &&
                est.sig_fep / est.eps_fep <= termination.target_fep_rel_precision) {
                stop_reason = StopReason::FepPrecision;
                stop_flag.store(true, std::memory_order_release);
            }

            if (termination.target_total_rel_precision > 0.0 && est.eps_tot > 0.0 &&
                est.sig_tot / est.eps_tot <= termination.target_total_rel_precision) {
                stop_reason = StopReason::TotalPrecision;
                stop_flag.store(true, std::memory_order_release);
            }
        }

        // Progress callback (throttled to at most once per second; the
        // guaranteed completion call is fired by compute() from the final
        // result, so no stop-triggered fire here).
        if (progress_callback) {
            double frac = 0.0;
            if (termination.max_events > 0)
                frac = std::max(frac, static_cast<double>(merged.num_events) /
                                      static_cast<double>(termination.max_events));
            if (termination.max_wall_seconds > 0.0)
                frac = std::max(frac, elapsed_sec / termination.max_wall_seconds);
            if (termination.max_cpu_seconds > 0.0)
                frac = std::max(frac, cpu_seconds / termination.max_cpu_seconds);
            // For precision-based, estimate fraction from current precision vs
            // target: precision scales as 1/sqrt(N), so frac ~ (target/current)^2.
            if (termination.target_fep_rel_precision > 0.0 && est.eps_fep > 0.0) {
                double rel = est.sig_fep / est.eps_fep;
                double ratio = termination.target_fep_rel_precision / std::max(rel, 1e-30);
                frac = std::max(frac, std::min(ratio * ratio, 1.0));
            }
            if (termination.target_total_rel_precision > 0.0 && est.eps_tot > 0.0) {
                double rel = est.sig_tot / est.eps_tot;
                double ratio = termination.target_total_rel_precision / std::max(rel, 1e-30);
                frac = std::max(frac, std::min(ratio * ratio, 1.0));
            }
            frac = std::min(frac, 1.0);

            auto since_last = std::chrono::duration<double>(now - last_callback_time).count();
            if (since_last >= 1.0) {
                last_callback_time = now;
                Progress p;
                p.num_events = merged.num_events;
                p.elapsed = elapsed;
                p.frac_complete = frac;
                p.fep_efficiency = est.eps_fep;
                p.fep_uncertainty = est.sig_fep;
                p.total_efficiency = est.eps_tot;
                p.total_uncertainty = est.sig_tot;
                progress_callback(p);
            }
        }
    }
};

} // anonymous namespace

EfficiencyCalculator::EfficiencyCalculator() = default;

void EfficiencyCalculator::set_detector(DetectorShape type, const Material* material,
                                        const std::vector<double>& dimensions) {
    geometry_.set_detector(type, material, dimensions);
}

void EfficiencyCalculator::set_bore_hole(double bore_radius, double bore_depth) {
    geometry_.set_bore_hole(bore_radius, bore_depth);
}

void EfficiencyCalculator::set_dead_layer(double front, double side, double back) {
    geometry_.set_dead_layer(front, side, back);
}

void EfficiencyCalculator::add_attenuator(const Material* material,
                                          double front_thickness, double side_thickness,
                                          double z_start, double z_end) {
    geometry_.add_attenuator(material, front_thickness, side_thickness, z_start, z_end);
}

void EfficiencyCalculator::add_collimator(const Material* material, double side_thickness,
                                          double z_start, double z_end) {
    geometry_.add_collimator(material, side_thickness, z_start, z_end);
}

void EfficiencyCalculator::set_point_source(const Eigen::Vector3d& position) {
    source_type_ = SourceType::Point;
    source_position_ = position;
    source_geometry_.configure_point(position);
}

void EfficiencyCalculator::set_cylindrical_source(
    const Eigen::Vector3d& center,
    double radius, double half_length,
    const Eigen::Matrix3d& rotation,
    double inner_radius)
{
    source_type_ = SourceType::Cylindrical;
    cyl_src_.center = center;
    cyl_src_.radius = radius;
    cyl_src_.inner_radius = inner_radius;
    cyl_src_.half_length = half_length;
    cyl_src_.rotation = rotation;
    source_position_ = center;
    source_geometry_.configure_cylindrical(center, radius, half_length, rotation,
                                           inner_radius);
}

void EfficiencyCalculator::set_spherical_source(
    const Eigen::Vector3d& center, double radius,
    const Eigen::Matrix3d& rotation, double inner_radius)
{
    source_type_ = SourceType::Spherical;
    sph_src_.center = center;
    sph_src_.inner_radius = inner_radius;
    sph_src_.radius = radius;
    sph_src_.rotation = rotation;
    source_position_ = center;
    source_geometry_.configure_spherical(center, inner_radius, radius, rotation);
}

void EfficiencyCalculator::set_rectangular_source(
    const Eigen::Vector3d& center,
    const Eigen::Vector3d& half_dims,
    const Eigen::Matrix3d& rotation,
    const Eigen::Vector3d& inner_half_dims)
{
    // Inner void must be all-zero (solid) or strictly inside the outer box.
    assert((inner_half_dims.array() == 0.0).all()
           || ((inner_half_dims.array() >= 0.0).all()
               && (inner_half_dims.array() < half_dims.array()).all()));
    source_type_ = SourceType::Rectangular;
    rect_src_.center = center;
    rect_src_.half_dims = half_dims;
    rect_src_.inner_half_dims = inner_half_dims;
    rect_src_.rotation = rotation;
    source_position_ = center;
    source_geometry_.configure_rectangular(center, half_dims, rotation,
                                           inner_half_dims);
}

void EfficiencyCalculator::set_marinelli_beaker(
    double well_inner_radius, double well_depth,
    double outer_radius, double fill_height,
    double endcap_to_beaker,
    const Material* sample_material,
    const Material* beaker_material, double beaker_thickness)
{
    source_type_ = SourceType::Marinelli;
    marinelli_cfg_.well_inner_radius = well_inner_radius;
    marinelli_cfg_.well_depth = well_depth;
    marinelli_cfg_.outer_radius = outer_radius;
    marinelli_cfg_.fill_height = fill_height;
    marinelli_cfg_.endcap_to_beaker = endcap_to_beaker;
    marinelli_cfg_.beaker_material = beaker_material;
    marinelli_cfg_.beaker_thickness = beaker_thickness;

    // Compute absolute z-coordinates from detector geometry
    double z_det_min = geometry_.outer_z_extent().first;
    double z_bk = z_det_min - endcap_to_beaker;
    double z_we = z_det_min + well_depth;
    double z_bot = z_bk - fill_height;

    // Source position = center of the L-shaped volume (approximate, for cone direction)
    source_position_ = Eigen::Vector3d(0.0, 0.0, (z_bot + z_we) / 2.0);

    // Cone sampling doesn't make sense for Marinelli (source surrounds detector)
    use_cone_sampling_ = false;

    // Configure source geometry with absolute z-coordinates
    source_geometry_.configure_marinelli(well_inner_radius, outer_radius,
                                         z_bk, z_we, z_bot);
    source_geometry_.set_source_material(sample_material);
    if (beaker_material && beaker_thickness > 0.0) {
        source_geometry_.add_shield(beaker_material, beaker_thickness);
    }
}

void EfficiencyCalculator::set_exponential_depth_distribution(double relaxation_length) {
    assert(relaxation_length > 0.0 && "Relaxation length must be positive");
    depth_distribution_ = DepthDistribution::Exponential;
    relaxation_length_ = relaxation_length;
}

void EfficiencyCalculator::set_uniform_depth_distribution() {
    depth_distribution_ = DepthDistribution::Uniform;
}

void EfficiencyCalculator::enable_source_electron_transport(bool on) {
    source_electron_transport_ = on;
    source_geometry_.set_source_electron_transport(on);
}

double EfficiencyCalculator::sample_exponential_depth(double total_depth, double U) const {
    double ratio = total_depth / relaxation_length_;
    if (ratio < 1e-3) {
        // D << lambda: nearly uniform, fall back to uniform sampling
        return total_depth * U;
    }
    double exp_neg_ratio = std::exp(-ratio);
    if (exp_neg_ratio < 1e-15) {
        // D >> lambda: truncation negligible, sample from untruncated exponential
        // Clamp to [0, D] for safety
        return std::min(-relaxation_length_ * std::log(1.0 - U), total_depth);
    }
    return -relaxation_length_ * std::log(1.0 - U * (1.0 - exp_neg_ratio));
}

void EfficiencyCalculator::set_source_material(const Material* mat) {
    source_geometry_.set_source_material(mat);
}

void EfficiencyCalculator::add_source_shield(const Material* mat, double thickness) {
    source_geometry_.add_shield(mat, thickness);
}

void EfficiencyCalculator::add_source_shield(const Material* mat, double t_radial, double t_end) {
    source_geometry_.add_shield(mat, t_radial, t_end);
}

void EfficiencyCalculator::add_source_shield(const Material* mat, double t_x, double t_y, double t_z) {
    source_geometry_.add_shield(mat, t_x, t_y, t_z);
}

// --- Position biasing ---

void EfficiencyCalculator::set_position_bias(const PositionBiasConfig& config) {
    position_bias_ = config;
    position_bias_enabled_ = true;
    position_bias_manual_ = true;
}

void EfficiencyCalculator::enable_position_bias() {
    position_bias_enabled_ = true;
    position_bias_manual_ = false;
}

void EfficiencyCalculator::disable_position_bias() {
    position_bias_enabled_ = false;
    position_bias_manual_ = false;
    position_bias_ = PositionBiasConfig{};
}

namespace {

/// Sample from a truncated exponential on [0, extent] with scale lambda.
/// Returns {sampled_value, weight} where weight = p_uniform / p_bias.
/// p_uniform = 1/extent, p_bias = (1/lambda)*exp(-x/lambda) / (1-exp(-extent/lambda)).
std::pair<double, double> sample_truncated_exp_bias(
    double extent, double lambda, double U)
{
    double ratio = extent / lambda;
    if (ratio < 1e-4) {
        // Lambda >> extent: bias is negligible, sample uniformly
        return {extent * U, 1.0};
    }
    double exp_neg_ratio = std::exp(-ratio);
    double x = -lambda * std::log(1.0 - U * (1.0 - exp_neg_ratio));
    double w = (1.0 - exp_neg_ratio) * lambda / (extent * std::exp(-x / lambda));
    return {x, w};
}

/// Sample from a truncated exponential on [0, extent] with bias scale lambda_b,
/// when the physical distribution is also exponential with scale lambda_p.
/// Returns {sampled_value, weight} where weight = p_physical / p_bias.
std::pair<double, double> sample_truncated_exp_bias_on_exp(
    double extent, double lambda_b, double lambda_p, double U)
{
    double ratio_b = extent / lambda_b;
    double exp_neg_b = std::exp(-ratio_b);
    double x = -lambda_b * std::log(1.0 - U * (1.0 - exp_neg_b));

    // weight = p_phys(x) / p_bias(x)
    // = (lb/lp) * (1-exp(-D/lb))/(1-exp(-D/lp)) * exp(x*(1/lb - 1/lp))
    double ratio_p = extent / lambda_p;
    double exp_neg_p = std::exp(-ratio_p);
    double w = (lambda_b / lambda_p)
               * (1.0 - exp_neg_b) / (1.0 - exp_neg_p)
               * std::exp(x * (1.0/lambda_b - 1.0/lambda_p));
    return {x, w};
}

} // anonymous namespace

std::pair<Eigen::Vector3d, double> EfficiencyCalculator::sample_biased_position(
    std::mt19937_64& rng,
    const PositionBiasConfig& bias) const
{
    // Point, Marinelli, Spherical: no position biasing (sphere has no preferred
    // lateral/depth axis — sample its volume distribution unbiased).
    if (source_type_ == SourceType::Point || source_type_ == SourceType::Marinelli
        || source_type_ == SourceType::Spherical) {
        return {sample_source_position(rng), 1.0};
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double w_lat = 1.0;
    double w_depth = 1.0;

    if (source_type_ == SourceType::Cylindrical) {
        double R = cyl_src_.radius;
        double H = cyl_src_.half_length;
        double D = 2.0 * H;

        // --- Lateral biasing ---
        double r;
        if (bias.solid_angle_h_cm > 0.0 && R > 0.0) {
            // Solid-angle-matched sampling: p(r) ∝ r/(r²+h²)
            // CDF inverse: r = h * sqrt(exp(U * ln(base)) - 1)
            double h = bias.solid_angle_h_cm;
            double log_base = bias.log_solid_angle_base;
            if (log_base <= 0.0) {
                // Fallback if not precomputed
                log_base = std::log(1.0 + (R * R) / (h * h));
            }
            double U = uniform(rng);
            r = h * std::sqrt(std::exp(U * log_base) - 1.0);
            // Weight = (r²+h²) * ln(1+R²/h²) / R²
            w_lat = (r * r + h * h) * log_base / (R * R);
        } else if (bias.lambda_lateral_cm > 0.0 && R > 0.0) {
            // Exponential-in-r² biasing (fallback for smaller sources)
            double beta = (R * R) / (2.0 * bias.lambda_lateral_cm * bias.lambda_lateral_cm);
            if (beta < 1e-4) {
                r = R * std::sqrt(uniform(rng));
            } else {
                double exp_neg_beta = std::exp(-beta);
                double u = -std::log(1.0 - uniform(rng) * (1.0 - exp_neg_beta)) / beta;
                r = R * std::sqrt(u);
                w_lat = (1.0 - exp_neg_beta) / (beta * std::exp(-beta * u));
            }
        } else {
            r = R * std::sqrt(uniform(rng));
        }
        double phi = kTwoPi * uniform(rng);

        // --- Depth biasing (z_local) ---
        // Determine depth direction from detector-frame projection
        Eigen::Vector3d depth_dir = cyl_src_.rotation * Eigen::Vector3d(0.0, 0.0, 1.0);
        double dz = depth_dir.z();  // local-z contribution to detector depth

        double z_local;
        bool apply_depth_bias = bias.lambda_depth_cm > 0.0 && std::abs(dz) > 0.1;

        if (apply_depth_bias) {
            double lambda_z_eff = bias.lambda_depth_cm / std::abs(dz);

            // d=0 at the detector-facing surface
            // If dz > 0: face is at z_local=+H, d = H - z_local
            // If dz < 0: face is at z_local=-H, d = z_local + H
            if (depth_distribution_ == DepthDistribution::Exponential) {
                auto [d, w] = sample_truncated_exp_bias_on_exp(
                    D, lambda_z_eff, relaxation_length_, uniform(rng));
                w_depth = w;
                z_local = (dz > 0) ? (H - d) : (-H + d);
            } else {
                auto [d, w] = sample_truncated_exp_bias(D, lambda_z_eff, uniform(rng));
                w_depth = w;
                z_local = (dz > 0) ? (H - d) : (-H + d);
            }
        } else if (depth_distribution_ == DepthDistribution::Exponential) {
            double d = sample_exponential_depth(D, uniform(rng));
            z_local = H - d;
        } else {
            z_local = D * (uniform(rng) - 0.5);
        }

        Eigen::Vector3d local_pos(r * std::cos(phi), r * std::sin(phi), z_local);
        Eigen::Vector3d pos = cyl_src_.rotation.transpose() * local_pos + cyl_src_.center;
        double w = std::min(w_lat * w_depth, 100.0);
        return {pos, w};
    }

    if (source_type_ == SourceType::Rectangular) {
        Eigen::Vector3d half = rect_src_.half_dims;

        // --- Lateral biasing (|x| and |y| independently) ---
        double x_local, y_local;
        double w_x = 1.0, w_y = 1.0;

        if (bias.lambda_lateral_cm > 0.0) {
            if (half.x() > 0.0) {
                auto [ax, w] = sample_truncated_exp_bias(half.x(), bias.lambda_lateral_cm, uniform(rng));
                x_local = (uniform(rng) < 0.5) ? ax : -ax;
                w_x = w;
            } else {
                x_local = 0.0;
            }
            if (half.y() > 0.0) {
                auto [ay, w] = sample_truncated_exp_bias(half.y(), bias.lambda_lateral_cm, uniform(rng));
                y_local = (uniform(rng) < 0.5) ? ay : -ay;
                w_y = w;
            } else {
                y_local = 0.0;
            }
            w_lat = w_x * w_y;
        } else {
            x_local = half.x() * (2.0 * uniform(rng) - 1.0);
            y_local = half.y() * (2.0 * uniform(rng) - 1.0);
        }

        // --- Depth biasing (z_local) ---
        Eigen::Vector3d depth_dir = rect_src_.rotation * Eigen::Vector3d(0.0, 0.0, 1.0);
        double dz = depth_dir.z();
        double D = 2.0 * half.z();

        double z_local;
        bool apply_depth_bias = bias.lambda_depth_cm > 0.0 && std::abs(dz) > 0.1;

        if (apply_depth_bias) {
            double lambda_z_eff = bias.lambda_depth_cm / std::abs(dz);

            if (depth_distribution_ == DepthDistribution::Exponential) {
                auto [d, w] = sample_truncated_exp_bias_on_exp(
                    D, lambda_z_eff, relaxation_length_, uniform(rng));
                w_depth = w;
                z_local = (dz > 0) ? (half.z() - d) : (-half.z() + d);
            } else {
                auto [d, w] = sample_truncated_exp_bias(D, lambda_z_eff, uniform(rng));
                w_depth = w;
                z_local = (dz > 0) ? (half.z() - d) : (-half.z() + d);
            }
        } else if (depth_distribution_ == DepthDistribution::Exponential) {
            double d = sample_exponential_depth(D, uniform(rng));
            z_local = half.z() - d;
        } else {
            z_local = half.z() * (2.0 * uniform(rng) - 1.0);
        }

        Eigen::Vector3d local_pos(x_local, y_local, z_local);
        Eigen::Vector3d pos = rect_src_.rotation.transpose() * local_pos + rect_src_.center;
        double w = std::min(w_lat * w_depth, 100.0);
        return {pos, w};
    }

    return {sample_source_position(rng), 1.0};
}

PositionBiasConfig EfficiencyCalculator::compute_auto_bias_params(double energy_keV) const {
    PositionBiasConfig bias;

    if (source_type_ == SourceType::Point || source_type_ == SourceType::Marinelli
        || source_type_ == SourceType::Spherical) {
        return bias;  // No biasing for these types (sphere has no preferred axis)
    }

    // Source-to-detector distance (approximate: from source center to detector front face)
    double z_det_front = geometry_.outer_z_extent().first;
    double distance = std::abs(source_position_.z() - z_det_front);

    // Lateral bias: scale with source-detector distance
    // Solid angle falls as ~1/(r² + h²), so characteristic scale is h
    double lateral_scale = distance / std::sqrt(2.0);

    // Adjust for offset: if source center is far off-axis, lateral biasing is less helpful
    double offset = std::sqrt(source_position_.x() * source_position_.x()
                            + source_position_.y() * source_position_.y());
    if (offset > 2.0 * lateral_scale) {
        lateral_scale = offset;  // Reduce aggressiveness
    }

    // Only enable lateral bias if source is significantly larger than the bias scale
    double source_lateral_extent = 0.0;
    if (source_type_ == SourceType::Cylindrical) {
        source_lateral_extent = cyl_src_.radius;
    } else if (source_type_ == SourceType::Rectangular) {
        source_lateral_extent = std::max(rect_src_.half_dims.x(), rect_src_.half_dims.y());
    }

    if (source_lateral_extent > 2.0 * lateral_scale
        && source_type_ == SourceType::Cylindrical) {
        // For cylindrical sources with R >> h, use solid-angle-matched sampling:
        // p(r) ∝ r/(r²+h²) naturally matches the detector solid angle falloff.
        bias.solid_angle_h_cm = distance;
        double R = cyl_src_.radius;
        double h = distance;
        bias.log_solid_angle_base = std::log(1.0 + (R * R) / (h * h));
    } else if (source_lateral_extent > 2.0 * lateral_scale) {
        // Rectangular or moderate-size cylindrical: exponential lateral bias
        constexpr double beta_max = 200.0;
        double beta = source_lateral_extent * source_lateral_extent
                     / (2.0 * lateral_scale * lateral_scale);
        if (beta > beta_max) {
            lateral_scale = source_lateral_extent / std::sqrt(2.0 * beta_max);
        }
        bias.lambda_lateral_cm = lateral_scale;
    }

    // Depth bias: scale with medium attenuation length
    const Material* medium = source_geometry_.source_material();
    if (medium) {
        double energy_MeV = energy_keV * 1e-3;
        MacroscopicXS xs = medium->macroscopic_xs(energy_MeV);
        double mu = xs.mu_total();
        if (mu > 0.0) {
            double lambda_atten = 1.0 / mu;
            if (depth_distribution_ == DepthDistribution::Exponential) {
                // Combined effective depth: 1/(1/lambda_relax + mu)
                bias.lambda_depth_cm = 1.0 / (1.0 / relaxation_length_ + mu);
            } else {
                bias.lambda_depth_cm = lambda_atten;
            }

            // Only enable if source depth is significantly larger than bias scale
            double source_depth = 0.0;
            if (source_type_ == SourceType::Cylindrical)
                source_depth = 2.0 * cyl_src_.half_length;
            else if (source_type_ == SourceType::Rectangular)
                source_depth = 2.0 * rect_src_.half_dims.z();

            if (source_depth < 2.0 * bias.lambda_depth_cm) {
                bias.lambda_depth_cm = 0.0;  // Not worth biasing
            }
        }
    }

    return bias;
}

// --- Source position sampling ---

Eigen::Vector3d EfficiencyCalculator::sample_source_position(std::mt19937_64& rng) const {
    if (source_type_ == SourceType::Point) {
        return source_position_;
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    if (source_type_ == SourceType::Cylindrical) {
        // Annular (uniform in area between inner bore and outer radius);
        // inner_radius == 0 reduces to r = radius·sqrt(U).
        double r_in2  = cyl_src_.inner_radius * cyl_src_.inner_radius;
        double r_out2 = cyl_src_.radius * cyl_src_.radius;
        double r   = std::sqrt(r_in2 + (r_out2 - r_in2) * uniform(rng));
        double phi = kTwoPi * uniform(rng);

        double z_local;
        if (depth_distribution_ == DepthDistribution::Exponential) {
            // Depth d=0 at detector-facing surface (z_local = +half_length),
            // increasing away from detector (toward z_local = -half_length).
            // With identity rotation and center_z < 0, z_local = +half_length
            // maps to the most positive z in detector frame (nearest z=0).
            double D = 2.0 * cyl_src_.half_length;
            double d = sample_exponential_depth(D, uniform(rng));
            z_local = cyl_src_.half_length - d;
        } else {
            z_local = 2.0 * cyl_src_.half_length * (uniform(rng) - 0.5);
        }

        Eigen::Vector3d local_pos(r * std::cos(phi), r * std::sin(phi), z_local);
        return cyl_src_.rotation.transpose() * local_pos + cyl_src_.center;
    }

    if (source_type_ == SourceType::Rectangular) {
        // Hollow box shell: rejection-sample the outer-box proposal, rejecting
        // draws inside the inner void. Acceptance = 1 - Vi/Vo (weighted by the
        // depth distribution when exponential); the loop is a few RNG draws,
        // cheap vs transport. Solid boxes (inner all-zero) take exactly one
        // pass — bit-identical to the pre-hollow sampling.
        const Eigen::Vector3d& inner = rect_src_.inner_half_dims;
        const bool hollow = inner.minCoeff() > 0.0;

        Eigen::Vector3d local_pos;
        do {
            double z_local;
            if (depth_distribution_ == DepthDistribution::Exponential) {
                double D = 2.0 * rect_src_.half_dims.z();
                double d = sample_exponential_depth(D, uniform(rng));
                z_local = rect_src_.half_dims.z() - d;
            } else {
                z_local = rect_src_.half_dims.z() * (2.0 * uniform(rng) - 1.0);
            }

            local_pos = Eigen::Vector3d(
                rect_src_.half_dims.x() * (2.0 * uniform(rng) - 1.0),
                rect_src_.half_dims.y() * (2.0 * uniform(rng) - 1.0),
                z_local
            );
        } while (hollow
                 && std::abs(local_pos.x()) < inner.x()
                 && std::abs(local_pos.y()) < inner.y()
                 && std::abs(local_pos.z()) < inner.z());
        return rect_src_.rotation.transpose() * local_pos + rect_src_.center;
    }

    if (source_type_ == SourceType::Spherical) {
        double r_in = sph_src_.inner_radius;
        double r_out = sph_src_.radius;
        double r;
        if (depth_distribution_ == DepthDistribution::Exponential) {
            // Surface-inward exponential contamination: depth d from the outer
            // surface, r = r_out − d, clamped to the shell. Thin-skin (constant
            // shell-area) approximation — the curvature r² weight is neglected,
            // appropriate for surface contamination where the active skin ≪ r_out.
            double D = r_out - r_in;
            double d = sample_exponential_depth(D, uniform(rng));
            r = r_out - d;
        } else {
            // Uniform in the spherical shell volume: r = cbrt(r_in³+(r_out³−r_in³)U).
            double r_in3 = r_in * r_in * r_in;
            double r_out3 = r_out * r_out * r_out;
            r = std::cbrt(r_in3 + (r_out3 - r_in3) * uniform(rng));
        }
        // Isotropic direction → uniform over the sphere surface at radius r.
        double cos_t = 2.0 * uniform(rng) - 1.0;
        double sin_t = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
        double phi = kTwoPi * uniform(rng);
        Eigen::Vector3d local_pos(r * sin_t * std::cos(phi),
                                  r * sin_t * std::sin(phi),
                                  r * cos_t);
        return sph_src_.rotation.transpose() * local_pos + sph_src_.center;
    }

    if (source_type_ == SourceType::Marinelli) {
        const auto& m = marinelli_cfg_;
        double ring_depth = m.endcap_to_beaker + m.well_depth;
        double or2 = m.outer_radius * m.outer_radius;
        double wr2 = m.well_inner_radius * m.well_inner_radius;
        double V_disk = kPi * or2 * m.fill_height;
        double V_ring = kPi * (or2 - wr2) * ring_depth;
        double V_total = V_disk + V_ring;

        // Absolute z-coordinates
        double z_det_min = geometry_.outer_z_extent().first;
        double z_bk = z_det_min - m.endcap_to_beaker;

        if (uniform(rng) < V_disk / V_total) {
            // Disk: z ∈ [z_bk - fill_height, z_bk], r ∈ [0, outer_r]
            double r = m.outer_radius * std::sqrt(uniform(rng));
            double phi = kTwoPi * uniform(rng);
            double z = z_bk - m.fill_height * uniform(rng);
            return {r * std::cos(phi), r * std::sin(phi), z};
        } else {
            // Ring: z ∈ [z_bk, z_bk + ring_depth], r ∈ [well_r, outer_r]
            double r = std::sqrt(wr2 + (or2 - wr2) * uniform(rng));
            double phi = kTwoPi * uniform(rng);
            double z = z_bk + ring_depth * uniform(rng);
            return {r * std::cos(phi), r * std::sin(phi), z};
        }
    }

    return source_position_;
}

// --- Cone sampling ---

double EfficiencyCalculator::compute_raw_cone_half_angle(const Eigen::Vector3d& from_pos) const {
    double r = geometry_.outer_bounding_radius();
    auto [z_min, z_max] = geometry_.outer_z_extent();

    double sz = from_pos.z();
    if (sz >= z_min && sz <= z_max) {
        return kPi;
    }
    double dz_front = std::abs(sz - z_min);
    double dz_back  = std::abs(sz - z_max);
    if (std::min(dz_front, dz_back) < 1e-6) {
        return kPi;
    }

    // Cone center direction: from source toward origin (detector front face center).
    // Work with unnormalized vectors to avoid sqrt; use cos(angle) = A·B/(|A||B|).
    double cx = -from_pos.x(), cy = -from_pos.y(), cz = -from_pos.z();
    double c_len_sq = cx * cx + cy * cy + cz * cz;

    // Lateral direction from axis toward source.
    double dx = std::sqrt(from_pos.x() * from_pos.x()
                        + from_pos.y() * from_pos.y());
    double ux, uy;
    if (dx > 1e-10) {
        ux = from_pos.x() / dx;
        uy = from_pos.y() / dx;
    } else {
        ux = 1.0; uy = 0.0;
    }

    // Check 4 critical candidates: far rim and near rim at both z-ends.
    // The maximum angular deviation is always at one of these points.
    double min_cos = 1.0;
    const double zs[] = {z_min, z_max};
    for (double z : zs) {
        // Far rim (away from source laterally) and near rim (toward source)
        for (int sign = -1; sign <= 1; sign += 2) {
            double px = sign * ux * r - from_pos.x();
            double py = sign * uy * r - from_pos.y();
            double pz = z - from_pos.z();
            double dot = cx * px + cy * py + cz * pz;
            double p_len_sq = px * px + py * py + pz * pz;
            // cos(angle) = dot / sqrt(c_len_sq * p_len_sq)
            double cos_a = dot / std::sqrt(c_len_sq * p_len_sq);
            min_cos = std::min(min_cos, cos_a);
        }
    }

    // Also check perpendicular rim points (important for sources near the axis)
    for (double z : zs) {
        double px = -uy * r - from_pos.x();
        double py =  ux * r - from_pos.y();
        double pz = z - from_pos.z();
        double dot = cx * px + cy * py + cz * pz;
        double p_len_sq = px * px + py * py + pz * pz;
        double cos_a = dot / std::sqrt(c_len_sq * p_len_sq);
        min_cos = std::min(min_cos, cos_a);
    }

    min_cos = std::clamp(min_cos, -1.0, 1.0);
    return std::min(std::acos(min_cos), kPi);
}

double EfficiencyCalculator::compute_cone_half_angle(const Eigen::Vector3d& from_pos) const {
    double raw = compute_raw_cone_half_angle(from_pos);
    return std::min(raw * 1.025, kPi);
}

double EfficiencyCalculator::compute_worst_case_cone() const {
    // Marinelli source surrounds detector — always isotropic
    if (source_type_ == SourceType::Marinelli) {
        return kPi;
    }

    std::vector<Eigen::Vector3d> candidates;

    if (source_type_ == SourceType::Cylindrical) {
        const auto& c = cyl_src_;
        const auto& R = c.rotation.transpose();
        // The worst-case cone angle occurs not just at r=0 and r=R, but at
        // intermediate radii where off-axis geometry creates wider angles to
        // the detector's far rim. Sample multiple radii to capture this.
        const double phis[] = {0.0, kPi / 2.0, kPi, 3.0 * kPi / 2.0};
        const double zs[] = {-c.half_length, 0.0, c.half_length};

        // Intermediate radii: the worst-case often occurs at r ~ few × R_det —
        // AND at interior radii 0 < r < R_det. For a source plane close to the
        // detector face, a position at r ≈ 0.5·R_det needs a cone reaching the
        // near rim at >120° from its (tilted) toward-origin axis, far wider
        // than either the r=0 or the r=R_det candidate (found July 2026 by the
        // det-response campaign S9: the {0, ≥R_det} candidate set clipped the
        // detector for a 10 cm trace disk 1 cm from a 3"x3" NaI face, biasing
        // eps_tot −9.5%; interior radii close the hole).
        double R_det = geometry_.outer_bounding_radius();
        std::vector<double> radii = {0.0};
        for (double scale : {0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0}) {
            double r = scale * R_det;
            if (r < c.radius) radii.push_back(r);
        }
        radii.push_back(c.radius);

        for (double z : zs) {
            for (double r_cand : radii) {
                if (r_cand < 1e-6) {
                    candidates.push_back(R * Eigen::Vector3d(0.0, 0.0, z) + c.center);
                } else {
                    for (double phi : phis) {
                        Eigen::Vector3d local(r_cand * std::cos(phi), r_cand * std::sin(phi), z);
                        candidates.push_back(R * local + c.center);
                    }
                }
            }
        }
    } else if (source_type_ == SourceType::Rectangular) {
        const auto& r = rect_src_;
        const auto& R = r.rotation.transpose();
        // Corners, face centers, and center — PLUS interior scalings of each,
        // for the same reason as the cylindrical interior radii (positions
        // between the axis and the detector rim need the widest cones).
        std::vector<Eigen::Vector3d> locals;
        for (int ix = -1; ix <= 1; ix += 2)
            for (int iy = -1; iy <= 1; iy += 2)
                for (int iz = -1; iz <= 1; iz += 2)
                    locals.emplace_back(ix * r.half_dims.x(),
                                        iy * r.half_dims.y(),
                                        iz * r.half_dims.z());
        locals.emplace_back( r.half_dims.x(), 0.0, 0.0);
        locals.emplace_back(-r.half_dims.x(), 0.0, 0.0);
        locals.emplace_back(0.0,  r.half_dims.y(), 0.0);
        locals.emplace_back(0.0, -r.half_dims.y(), 0.0);
        locals.emplace_back(0.0, 0.0,  r.half_dims.z());
        locals.emplace_back(0.0, 0.0, -r.half_dims.z());
        for (const auto& l : locals)
            for (double f : {0.25, 0.5, 0.75, 1.0})
                candidates.push_back(R * (f * l) + r.center);
        candidates.push_back(r.center);
    } else if (source_type_ == SourceType::Spherical) {
        const auto& s = sph_src_;
        const double R = s.radius;
        candidates.push_back(s.center);
        // 6 axis + 8 diagonal points, at the surface AND at interior radii
        // (a sphere is rotation-invariant, so no orientation sweep is needed;
        // interior points needed for the same reason as the cylinder).
        const double d = R / std::sqrt(3.0);
        for (double f : {0.25, 0.5, 0.75, 1.0}) {
            for (int i = -1; i <= 1; i += 2) {
                candidates.push_back(s.center + f * Eigen::Vector3d(i * R, 0.0, 0.0));
                candidates.push_back(s.center + f * Eigen::Vector3d(0.0, i * R, 0.0));
                candidates.push_back(s.center + f * Eigen::Vector3d(0.0, 0.0, i * R));
            }
            for (int ix = -1; ix <= 1; ix += 2)
                for (int iy = -1; iy <= 1; iy += 2)
                    for (int iz = -1; iz <= 1; iz += 2)
                        candidates.push_back(
                            s.center + f * Eigen::Vector3d(ix * d, iy * d, iz * d));
        }
    }

    double max_raw = 0.0;
    for (const auto& pos : candidates) {
        max_raw = std::max(max_raw, compute_raw_cone_half_angle(pos));
    }

    return std::min(max_raw * 1.025, kPi);
}

Eigen::Vector3d EfficiencyCalculator::sample_cone_direction(
    const Eigen::Vector3d& from_pos,
    double cos_theta_max, std::mt19937_64& rng) const
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double cos_theta = 1.0 - uniform(rng) * (1.0 - cos_theta_max);
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double phi = kTwoPi * uniform(rng);

    Eigen::Vector3d cone_dir(
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta
    );

    Eigen::Vector3d to_detector = -from_pos;
    to_detector.normalize();

    Eigen::Vector3d z_axis(0.0, 0.0, 1.0);
    double cos_rot = z_axis.dot(to_detector);

    if (cos_rot > 0.9999) {
        return cone_dir;
    }
    if (cos_rot < -0.9999) {
        return Eigen::Vector3d(-cone_dir.x(), cone_dir.y(), -cone_dir.z());
    }

    Eigen::Vector3d axis = z_axis.cross(to_detector);
    axis.normalize();
    double sin_rot = std::sqrt(1.0 - cos_rot * cos_rot);

    Eigen::Vector3d result = cone_dir * cos_rot
                           + axis.cross(cone_dir) * sin_rot
                           + axis * (axis.dot(cone_dir)) * (1.0 - cos_rot);
    result.normalize();
    return result;
}

std::pair<Eigen::Vector3d, double> EfficiencyCalculator::sample_surface_direction(
    const Eigen::Vector3d& src_pos, std::mt19937_64& rng) const
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    auto [z_min, z_max] = geometry_.outer_z_extent();
    double L = z_max - z_min;

    Eigen::Vector3d target;
    Eigen::Vector3d outward_normal;

    if (geometry_.shape() == DetectorShape::Cylinder) {
        double R = geometry_.outer_bounding_radius();
        double A_front = kPi * R * R;
        double A_side  = kTwoPi * R * L;
        double A_total = 2.0 * A_front + A_side;

        double u = uniform(rng) * A_total;
        if (u < A_front) {
            // Front disk (z = z_min)
            double r = R * std::sqrt(uniform(rng));
            double phi = kTwoPi * uniform(rng);
            target = {r * std::cos(phi), r * std::sin(phi), z_min};
            outward_normal = {0.0, 0.0, -1.0};
        } else if (u < 2.0 * A_front) {
            // Back disk (z = z_max)
            double r = R * std::sqrt(uniform(rng));
            double phi = kTwoPi * uniform(rng);
            target = {r * std::cos(phi), r * std::sin(phi), z_max};
            outward_normal = {0.0, 0.0, 1.0};
        } else {
            // Side cylinder
            double phi = kTwoPi * uniform(rng);
            double z = z_min + uniform(rng) * L;
            double cp = std::cos(phi), sp = std::sin(phi);
            target = {R * cp, R * sp, z};
            outward_normal = {cp, sp, 0.0};
        }

        Eigen::Vector3d to_src = src_pos - target;
        double dist = to_src.norm();
        double cos_inc = to_src.dot(outward_normal) / dist;

        if (cos_inc <= 0.0) {
            return {{0.0, 0.0, 0.0}, 0.0};
        }

        Eigen::Vector3d direction = (target - src_pos) / dist;
        double weight = A_total * cos_inc / (kFourPi * dist * dist);
        return {direction, weight};
    }

    // Box detector
    double hx = geometry_.detector_half_x();
    double hy = geometry_.detector_half_y();
    // Accumulate outer half-dims (same as ray_hits_outer_boundary)
    if (geometry_.has_dead_layer()) {
        auto dl = geometry_.dead_layer();
        if (dl) {
            hx += dl->side;
            hy += dl->side;
        }
    }
    for (const auto& att : geometry_.attenuators()) {
        hx += att.side_thickness;
        hy += att.side_thickness;
    }

    double A_front = 4.0 * hx * hy;              // front and back each
    double A_lr    = 2.0 * hy * L;                // left/right each
    double A_tb    = 2.0 * hx * L;                // top/bottom each
    double A_total = 2.0 * A_front + 2.0 * A_lr + 2.0 * A_tb;

    double u = uniform(rng) * A_total;
    if (u < A_front) {
        // Front (z = z_min)
        double x = hx * (2.0 * uniform(rng) - 1.0);
        double y = hy * (2.0 * uniform(rng) - 1.0);
        target = {x, y, z_min};
        outward_normal = {0.0, 0.0, -1.0};
    } else if (u < 2.0 * A_front) {
        // Back (z = z_max)
        double x = hx * (2.0 * uniform(rng) - 1.0);
        double y = hy * (2.0 * uniform(rng) - 1.0);
        target = {x, y, z_max};
        outward_normal = {0.0, 0.0, 1.0};
    } else if (u < 2.0 * A_front + A_lr) {
        // Right (+x face)
        double y = hy * (2.0 * uniform(rng) - 1.0);
        double z = z_min + uniform(rng) * L;
        target = {hx, y, z};
        outward_normal = {1.0, 0.0, 0.0};
    } else if (u < 2.0 * A_front + 2.0 * A_lr) {
        // Left (-x face)
        double y = hy * (2.0 * uniform(rng) - 1.0);
        double z = z_min + uniform(rng) * L;
        target = {-hx, y, z};
        outward_normal = {-1.0, 0.0, 0.0};
    } else if (u < 2.0 * A_front + 2.0 * A_lr + A_tb) {
        // Top (+y face)
        double x = hx * (2.0 * uniform(rng) - 1.0);
        double z = z_min + uniform(rng) * L;
        target = {x, hy, z};
        outward_normal = {0.0, 1.0, 0.0};
    } else {
        // Bottom (-y face)
        double x = hx * (2.0 * uniform(rng) - 1.0);
        double z = z_min + uniform(rng) * L;
        target = {x, -hy, z};
        outward_normal = {0.0, -1.0, 0.0};
    }

    Eigen::Vector3d to_src = src_pos - target;
    double dist = to_src.norm();
    double cos_inc = to_src.dot(outward_normal) / dist;

    if (cos_inc <= 0.0) {
        return {{0.0, 0.0, 0.0}, 0.0};
    }

    Eigen::Vector3d direction = (target - src_pos) / dist;
    double weight = A_total * cos_inc / (kFourPi * dist * dist);
    return {direction, weight};
}

// --- Thread-level simulation (one batch) ---

EfficiencyCalculator::ThreadTally EfficiencyCalculator::simulate_thread(
    double energy_keV,
    uint64_t num_events,
    uint64_t seed,
    const std::vector<float>& energy_bin_edges) const
{
    ThreadTally tally;
    tally.num_events = num_events;

    if (energy_bin_edges.size() >= 2) {
        tally.pulse_height_counts.resize(energy_bin_edges.size() - 1, 0.0);
        tally.pulse_height_counts_sq.resize(energy_bin_edges.size() - 1, 0.0);
    }

    const uint64_t seed_mixed = seed * 2654435761ULL;
    const std::array<uint_least32_t, 6> seed_data{{
      static_cast<uint_least32_t>(seed), static_cast<uint_least32_t>(seed >> 32),
      static_cast<uint_least32_t>(seed_mixed), static_cast<uint_least32_t>(seed_mixed >> 32),
      static_cast<uint_least32_t>(seed ^ (seed >> 16)),
      static_cast<uint_least32_t>((seed >> 48) | (seed << 16))
    }};
    std::seed_seq seq( seed_data.begin(), seed_data.end() );
    std::mt19937_64 rng(seq);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double energy_MeV = energy_keV * 1e-3;

    // Pre-compute air attenuation coefficient (excluding Rayleigh) once per thread.
    double mu_air_no_rs = 0.0;
    if (air_attenuation_ == AirAttenuation::AnalyticNoScatter) {
        static const Material air_mat = make_Air();
        MacroscopicXS air_xs = air_mat.macroscopic_xs(energy_MeV);
        mu_air_no_rs = air_xs.mu_pe + air_xs.mu_cs + air_xs.mu_pp;
    }

    // Pre-compute forced-first-interaction transport config (variance
    // reduction). The history weight 1 - exp(-tau) is returned in
    // TransportResult::weight and multiplied into the event score weight.
    // Gated off for the Marinelli primary path: its re-entry physics needs
    // the un-collided pass-through branch (and the event-level weight could
    // not represent per-branch forcing).
    TransportConfig forced_transport_config = transport_config_;
    forced_transport_config.force_first_interaction = true;
    const bool force_interaction = active_biasing_.force_detector_interaction;
    const bool force_src_primary =
        force_interaction && source_type_ != SourceType::Marinelli;

    // Pre-compute position bias parameters (constant across all events in this batch).
    PositionBiasConfig active_bias;
    if (position_bias_enabled_) {
        if (position_bias_manual_) {
            active_bias = position_bias_;
        } else {
            active_bias = compute_auto_bias_params(energy_keV);
        }
    }

    // Pre-compute cone parameters.
    // For point sources: cone from the fixed source position (constant).
    // For extended sources: per-position cone if beneficial, else worst-case.
    double pre_cone_half_angle = kPi;
    double pre_cos_theta_max   = -1.0;
    double pre_event_weight    = 1.0;
    bool per_position_cone = false;
    if (use_cone_sampling_) {
        if (source_type_ == SourceType::Point) {
            pre_cone_half_angle = compute_cone_half_angle(source_position_);
        } else {
            // For solid-angle-matched sampling, per-position cones are
            // essential: the radial IS weight q(r) ∝ r/(r²+h²) is designed
            // to cancel the per-position cone weight Ω(r)/(4π), producing
            // constant total weight per hit event. A fixed cone would break
            // this cancellation and create enormous weight variance.
            if (active_bias.solid_angle_h_cm > 0.0) {
                per_position_cone = true;
            } else {
                // For other extended sources, check if per-position cones
                // are worthwhile by comparing worst vs best solid angles.
                double worst = compute_worst_case_cone();
                double best = compute_cone_half_angle(source_position_);
                double omega_worst = kTwoPi * (1.0 - std::cos(worst));
                double omega_best  = kTwoPi * (1.0 - std::cos(best));
                if (omega_worst > 2.0 * omega_best) {
                    per_position_cone = true;
                } else {
                    pre_cone_half_angle = worst;
                }
            }
        }
        pre_cos_theta_max = std::cos(pre_cone_half_angle);
        double solid_angle = kTwoPi * (1.0 - pre_cos_theta_max);
        pre_event_weight = solid_angle / kFourPi;
    }

    std::vector<PathSegment> event_segments;
    for (uint64_t event = 0; event < num_events; ++event) {
        double position_weight = 1.0;
        Eigen::Vector3d src_pos;
        if (position_bias_enabled_) {
            auto [pos, w] = sample_biased_position(rng, active_bias);
            src_pos = pos;
            position_weight = w;
        } else {
            src_pos = sample_source_position(rng);
        }

        // ================================================================
        // Full mode with source shielding: isotropic / mixture emission.
        //
        // Default (alpha = 0): emit isotropically, transport through source
        // material/shielding, then check if the photon hits the detector.
        // Photons that exit the source heading away from the detector
        // genuinely miss. Pure cone sampling would be biased here: photons
        // emitted at wide angles can scatter in the source material toward
        // the detector.
        //
        // Mixture angular biasing (alpha in (0,1)): emit from
        //   q(dir) = alpha*(1/4pi) + (1-alpha)*1[dir in cone]/Omega_c,
        // weight w = (1/4pi)/q(dir). Unbiased for any downstream physics
        // including source scatter, because q > 0 over all 4pi and the
        // weight depends only on the initial emission direction. Bounded:
        // w <= 1/alpha. Direct-component statistics improve by up to
        // Omega_c/4pi; scattered-in photons keep alpha-fraction coverage.
        // ================================================================
        if (source_geometry_.has_source_effects() && !fep_only_mode_) {
            // -------- Two-stream stratified estimator (BiasingConfig::two_stream)
            // Deterministic round-robin allocation over 20-event cycles;
            // direct events carry weight factor 1/f, scatter events
            // 1/(1-f), so the single-pool estimator sum(w_i I_i)/N is
            // unbiased for every tally (FEP, total, spectrum). The variance
            // reported by the standard formula is slightly conservative for
            // this stratified allocation (it includes a between-strata term
            // that deterministic allocation does not have).
            const bool two_stream = active_biasing_.two_stream &&
                source_type_ != SourceType::Marinelli;
            bool stream_direct = false;
            double w_stream = 1.0;
            if (two_stream) {
                constexpr int kStreamSlots = 20;
                const int d_slots = std::clamp(
                    static_cast<int>(std::lround(
                        active_biasing_.two_stream_direct_fraction *
                        kStreamSlots)),
                    1, kStreamSlots - 1);
                stream_direct =
                    (event % kStreamSlots) < static_cast<uint64_t>(d_slots);
                w_stream = stream_direct
                    ? static_cast<double>(kStreamSlots) / d_slots
                    : static_cast<double>(kStreamSlots) /
                          (kStreamSlots - d_slots);
            }

            if (stream_direct) {
                // ---- Direct stream: the unscattered-primary (u) term ----
                // Emit in the detector-subtending cone and weight by
                // (Omega/4pi) * T, where T = exp(-sum mu_tot*l) is the
                // deterministic probability that the primary crosses the
                // source geometry with zero interactions (same unmasked
                // mu_total as the analog sampler, so T is exact).
                // Conditional on zero interactions the photon keeps its
                // direction and energy and produces no secondaries, so
                // plain detector transport estimates the u-term of every
                // tally. Directions outside the cone provably score zero
                // (the cone bounds the full outer detector envelope with a
                // 1.025 margin). The primary is always the sole depositing
                // branch here, so forced first interaction is always valid.
                tally.sum_weights += 1.0;
                tally.dbg_src.n_u++;

                Eigen::Vector3d direction;
                double omega_frac = 1.0;
                bool coned = false;
                if (use_cone_sampling_) {
                    const double cone_half = compute_cone_half_angle(src_pos);
                    if (cone_half < kPi - 0.01) {
                        const double cos_max = std::cos(cone_half);
                        direction = sample_cone_direction(src_pos, cos_max, rng);
                        omega_frac = 0.5 * (1.0 - cos_max);
                        coned = true;
                    }
                }
                if (!coned) {
                    double ct = 2.0 * uniform(rng) - 1.0;
                    double st = std::sqrt(1.0 - ct * ct);
                    double ph = kTwoPi * uniform(rng);
                    direction = Eigen::Vector3d(
                        st * std::cos(ph), st * std::sin(ph), ct);
                }

                double src_path_cm = 0.0;
                const double T = source_geometry_.no_interaction_probability(
                    src_pos, direction, energy_keV, &src_path_cm);
                double w = position_weight * w_stream * omega_frac * T;
                if (w <= 0.0) continue;

                geometry_.trace_ray(src_pos, direction, event_segments);
                if (event_segments.empty()) continue;
                const PathSegment* first_seg = nullptr;
                for (const auto& seg : event_segments) {
                    if (seg.material) { first_seg = &seg; break; }
                }
                if (!first_seg) continue;

                if (mu_air_no_rs > 0.0) {
                    // Air gap = distance to detector entry minus the path
                    // inside the source geometry along this ray.
                    double air_dist =
                        std::max(first_seg->t_start, 0.0) - src_path_cm;
                    if (air_dist > 0.0)
                        w *= std::exp(-mu_air_no_rs * air_dist);
                }

                Eigen::Vector3d entry_point =
                    src_pos + direction * std::max(first_seg->t_start, 0.0);
                auto tr = transport_photon(
                    entry_point, direction, energy_keV, geometry_,
                    force_interaction ? forced_transport_config
                                      : transport_config_,
                    rng);
                w *= tr.weight;
                if (tr.forced_absorption) tally.num_forced_absorption++;

                if (tr.any_interaction_in_scoring) {
                    tally.num_any++;
                    tally.sum_any_weights += w;
                    tally.sum_any_w_sq += w * w;
                    tally.dbg_src.any_w_u += w;
                    tally.dbg_src.any_w2_u += w * w;
                }
                if (std::abs(tr.energy_deposited_scoring - energy_keV) <
                    kFepTolerance) {
                    tally.num_fep++;
                    tally.sum_fep_weights += w;
                    tally.sum_fep_w_sq += w * w;
                    tally.dbg_src.fep_w_u += w;
                    tally.dbg_src.fep_w2_u += w * w;
                }
                if (!tally.pulse_height_counts.empty() &&
                    tr.energy_deposited_scoring > 0.0) {
                    auto it = std::upper_bound(
                        energy_bin_edges.begin(), energy_bin_edges.end(),
                        static_cast<float>(tr.energy_deposited_scoring));
                    if (it != energy_bin_edges.begin()) {
                        size_t bin = static_cast<size_t>(
                            it - energy_bin_edges.begin()) - 1;
                        if (bin < tally.pulse_height_counts.size()) {
                            tally.pulse_height_counts[bin] += w;
                            tally.pulse_height_counts_sq[bin] += w * w;
                        }
                    }
                }
                continue;
            }

            Eigen::Vector3d direction;
            double w_mix = 1.0;
            // Scatter stream emits purely isotropically (it scores the
            // s-term under the analog angular distribution); otherwise the
            // legacy mixture emission applies.
            const double alpha =
                two_stream ? 0.0 : active_biasing_.mixture_cone_alpha;
            bool mixture_done = false;

            if (alpha > 0.0 && alpha <= 1.0) {
                const double cone_half = compute_cone_half_angle(src_pos);
                // Degenerate cone (source beside/inside the detector
                // envelope): fall back to pure isotropic with weight 1.
                if (cone_half < kPi - 0.01) {
                    const double cos_max = std::cos(cone_half);
                    const double omega_c = kTwoPi * (1.0 - cos_max);
                    bool in_cone;
                    if (uniform(rng) < alpha) {
                        // Isotropic branch; membership test must use the
                        // same axis the cone sampler rotates onto:
                        // normalize(-src_pos), toward the detector origin.
                        double ct = 2.0 * uniform(rng) - 1.0;
                        double st = std::sqrt(1.0 - ct * ct);
                        double ph = kTwoPi * uniform(rng);
                        direction = Eigen::Vector3d(
                            st * std::cos(ph), st * std::sin(ph), ct);
                        const Eigen::Vector3d axis = -src_pos.normalized();
                        in_cone = direction.dot(axis) >= cos_max;
                    } else {
                        direction = sample_cone_direction(src_pos, cos_max, rng);
                        in_cone = true;
                    }
                    w_mix = in_cone
                        ? omega_c / (alpha * omega_c + (1.0 - alpha) * kFourPi)
                        : 1.0 / alpha;
                    mixture_done = true;
                }
            }

            if (!mixture_done) {
                double cos_theta = 2.0 * uniform(rng) - 1.0;
                double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                double phi = kTwoPi * uniform(rng);
                direction = Eigen::Vector3d(
                    sin_theta * std::cos(phi),
                    sin_theta * std::sin(phi),
                    cos_theta);
            }

            tally.sum_weights += 1.0;

            // Transport through source material + shielding
            auto src_result = source_geometry_.transport_source_photon(
                src_pos, direction, energy_keV, rng);

            // Save PP annihilation gammas from source material/shielding
            // (511 keV). These have already been transported through the
            // remaining source geometry by the recursive transport call and
            // are ready for trace_ray to the detector. All source shapes
            // produce them; the detector-transport lambda below is safe for
            // non-Marinelli sources (compute_marinelli_reentry reports no
            // re-entry, so the photon gets plain trace + crystal transport).
            std::vector<SourceGeometry::SourceSecondaryPhoton> source_pp =
                std::move(src_result.secondaries);

            // Classify the event by the primary's source-transport outcome
            // (u = zero interactions, s = interacted) for the variance
            // decomposition diagnostics. Counted before any early-out so the
            // class event counts sum to the number of source-effect events.
            // (In two-stream mode the direct stream counts n_u; killed
            // scatter-stream events are counted by neither class.)
            const bool src_interacted = src_result.interacted;
            if (src_interacted)    tally.dbg_src.n_s++;
            else if (!two_stream)  tally.dbg_src.n_u++;

            // Scatter stream scores ONLY the s-term: events whose primary
            // crossed the source geometry without interacting belong to the
            // direct stream's term. Killing them here is both required for
            // unbiasedness of the stream split and a per-event time win
            // (skips ray trace + detector transport).
            if (two_stream && !src_interacted) continue;

            // Absorbed-primary escaped-electron channel: the primary died in
            // the source geometry and no secondary photon survived, but
            // escaped recoil electrons can still reach the crystal (e.g.
            // Compton recoil escapes, then the degraded photon is
            // photo-absorbed in the shield). Tracked separately so the
            // channel's contribution is measurable (it was dropped entirely
            // before June 2026).
            const bool electron_only_event = !src_result.survived &&
                source_pp.empty() && !src_result.source_electrons.empty();
            if (electron_only_event) {
                tally.dbg_src.n_e_only++;
                tally.dbg_src.n_e_only_electrons +=
                    src_result.source_electrons.size();
                for (const auto& e : src_result.source_electrons)
                    tally.dbg_src.e_only_energy_keV += e.energy_keV;
            }

            if (!src_result.survived && source_pp.empty() &&
                src_result.source_electrons.empty()) continue;

            // Forced first interaction is only valid when the primary is the
            // sole depositing branch of the event: the forced-collision
            // weight multiplies the whole event score, and the FEP/total
            // indicators are non-linear in per-branch deposits. Decide per
            // event from the source-transport outcome — unbiased, since the
            // estimator choice conditions only on already-sampled source
            // randomness, never on detector randomness.
            const bool force_this_event = force_src_primary &&
                source_pp.empty() && src_result.source_electrons.empty();

            double total_dep_scoring = 0.0;
            bool had_any_scoring = false;
            // One weight per event: position bias x mixture/stream factors x
            // Compton-angle bias from source transport (1.0 when disabled).
            double src_air_w = position_weight * w_mix * w_stream *
                               src_result.bias_weight;
            bool air_applied = false;

            // Lambda: full Marinelli detector processing for a photon exiting
            // source water. Includes miss-bounce, crystal transport, post-crystal
            // re-entry bounce loop, and escaped secondary processing.
            // Used for both the primary photon and PP annihilation gammas.
            // Wall transport result: Passed = no interaction, Scattered = new
            // direction/energy from Compton/Rayleigh, Absorbed = PE/PP in wall.
            enum class WallResult { Passed, Scattered, Absorbed };

            // Full MC transport through beaker wall for re-entry.
            // On Compton/Rayleigh scatter, updates dir/eng in place.
            auto wall_transport = [&](double wall_path,
                Eigen::Vector3d& w_dir, double& w_eng) -> WallResult
            {
                if (wall_path <= 0.0 || source_geometry_.shields().empty())
                    return WallResult::Passed;
                const auto& wall_mat = *source_geometry_.shields()[0].material;
                double w_eng_MeV = w_eng * 1e-3;
                auto wall_xs = wall_mat.macroscopic_xs(w_eng_MeV);
                double mu_tot = wall_xs.mu_total();
                if (mu_tot <= 0.0) return WallResult::Passed;
                double s = -std::log(uniform(rng)) / mu_tot;
                if (s >= wall_path) return WallResult::Passed;
                auto type = select_interaction(wall_xs, uniform(rng));
                if (type == InteractionType::Photoelectric ||
                    type == InteractionType::PairProduction)
                    return WallResult::Absorbed;
                if (type == InteractionType::Compton) {
                    int Z = wall_mat.select_element(w_eng_MeV, 1, uniform(rng));
                    auto scatter = sample_compton_scatter(w_eng, w_dir, rng, Z);
                    w_eng = scatter.scattered_energy_keV;
                    w_dir = scatter.new_direction;
                    if (w_eng < 1.0) return WallResult::Absorbed;
                } else { // Rayleigh
                    int Z = wall_mat.select_element(w_eng_MeV, 2, uniform(rng));
                    double cos_theta = sample_rayleigh_cos_theta(Z, w_eng, rng);
                    double phi = 2 * M_PI * uniform(rng);
                    w_dir = rotate_direction(w_dir, cos_theta, phi);
                }
                return WallResult::Scattered;
            };

            // Per-event diagnostic accumulators for re-entry tracking
            double evt_dep_initial = 0.0;
            double evt_dep_reentry = 0.0;
            double evt_dep_secondary = 0.0;
            bool evt_had_initial_hit = false;
            bool evt_had_miss_bounce = false;
            uint64_t evt_reentry_hits = 0;
            uint64_t evt_secondary_hits = 0;

            auto marinelli_detector_transport = [&](
                Eigen::Vector3d pos, Eigen::Vector3d dir, double eng_keV,
                bool force_first = false)
            {
                auto segments = geometry_.trace_ray(pos, dir);

                // Miss-bounce: if photon exits water but misses crystal,
                // try re-entering water on another surface.
                bool used_miss_bounce = false;
                if (segments.empty()) {
                    constexpr int kMaxMissRetries = 5;
                    for (int retry = 0; retry < kMaxMissRetries && eng_keV > 1.0; ++retry) {
                        auto reentry = source_geometry_.compute_marinelli_reentry(pos, dir);
                        if (!reentry.can_reenter) break;
                        auto wr_wall = wall_transport(reentry.wall_path, dir, eng_keV);
                        if (wr_wall == WallResult::Absorbed) break;
                        if (wr_wall == WallResult::Scattered) {
                            pos = reentry.water_entry_pos;
                            continue;  // New direction — re-evaluate re-entry
                        }
                        auto wr = source_geometry_.transport_source_photon(
                            reentry.water_entry_pos, dir, eng_keV, rng);
                        if (!wr.survived) break;
                        segments = geometry_.trace_ray(wr.position, wr.direction);
                        if (!segments.empty()) {
                            pos = wr.position;
                            dir = wr.direction;
                            eng_keV = wr.energy_keV;
                            used_miss_bounce = true;
                            break;
                        }
                        pos = wr.position;
                        dir = wr.direction;
                        eng_keV = wr.energy_keV;
                    }
                }

                const PathSegment* first_seg = nullptr;
                if (!segments.empty()) {
                    for (const auto& seg : segments) {
                        if (seg.material) { first_seg = &seg; break; }
                    }
                }
                if (!first_seg) return;

                // Diagnostic: track initial hit type
                if (used_miss_bounce) evt_had_miss_bounce = true;
                else evt_had_initial_hit = true;

                Eigen::Vector3d entry_point = pos +
                    dir * std::max(first_seg->t_start, 0.0);

                if (mu_air_no_rs > 0.0 && !air_applied) {
                    double air_dist = std::max(first_seg->t_start, 0.0);
                    src_air_w *= std::exp(-mu_air_no_rs * air_dist);
                    air_applied = true;
                }

                auto tr = transport_photon(
                    entry_point, dir, eng_keV, geometry_,
                    force_first ? forced_transport_config : transport_config_,
                    rng);
                if (force_first) {
                    // One weight per event: the forced-collision weight of the
                    // primary's initial detector transport scales the whole
                    // event score. Only valid when the primary is the sole
                    // depositing branch — the caller guarantees this by
                    // passing force_first only when the event has no source
                    // secondaries/electrons (force_this_event).
                    src_air_w *= tr.weight;
                }

                if (tr.forced_absorption) tally.num_forced_absorption++;
                total_dep_scoring += tr.energy_deposited_scoring;
                evt_dep_initial += tr.energy_deposited_scoring;
                if (tr.any_interaction_in_scoring) had_any_scoring = true;

                // Post-crystal re-entry bounce loop.
                constexpr int kMaxReentries = 10;
                Eigen::Vector3d bounce_pos, bounce_dir;
                double bounce_energy;
                bool have_bounce = false;

                auto escaped_pool = std::move(tr.escaped_secondaries);

                if (tr.escaped && tr.exit_energy_keV > 1.0) {
                    bounce_pos = tr.exit_position;
                    bounce_dir = tr.exit_direction;
                    bounce_energy = tr.exit_energy_keV;
                    have_bounce = true;
                    tally.dbg_n_escaped++;
                }

                for (int re = 0; re < kMaxReentries && have_bounce; ++re) {
                    auto reentry = source_geometry_.compute_marinelli_reentry(
                        bounce_pos, bounce_dir);
                    if (!reentry.can_reenter) break;
                    tally.dbg_n_can_reenter++;

                    auto bw = wall_transport(reentry.wall_path, bounce_dir, bounce_energy);
                    if (bw == WallResult::Absorbed) { tally.dbg_n_wall_absorb++; break; }
                    if (bw == WallResult::Scattered) {
                        tally.dbg_n_wall_scatter++;
                        bounce_pos = reentry.water_entry_pos;
                        continue;  // New direction — re-evaluate re-entry
                    }
                    tally.dbg_n_wall_pass++;

                    auto wr = source_geometry_.transport_source_photon(
                        reentry.water_entry_pos, bounce_dir,
                        bounce_energy, rng);
                    if (!wr.survived) break;
                    tally.dbg_n_water_survived++;

                    auto re_segs = geometry_.trace_ray(wr.position, wr.direction);
                    if (re_segs.empty()) {
                        tally.dbg_n_trace_miss++;
                        bounce_pos = wr.position;
                        bounce_dir = wr.direction;
                        bounce_energy = wr.energy_keV;
                        continue;
                    }
                    tally.dbg_n_trace_hit++;

                    const PathSegment* re_first = nullptr;
                    for (const auto& seg : re_segs) {
                        if (seg.material) { re_first = &seg; break; }
                    }
                    if (!re_first) break;

                    Eigen::Vector3d re_entry = wr.position +
                        wr.direction * std::max(re_first->t_start, 0.0);

                    auto re_tr = transport_photon(re_entry, wr.direction, wr.energy_keV,
                                                  geometry_, transport_config_, rng);

                    total_dep_scoring += re_tr.energy_deposited_scoring;
                    evt_dep_reentry += re_tr.energy_deposited_scoring;
                    evt_reentry_hits++;
                    if (re_tr.any_interaction_in_scoring) had_any_scoring = true;

                    for (auto& ep : re_tr.escaped_secondaries)
                        escaped_pool.push_back(ep);

                    if (!re_tr.escaped || re_tr.exit_energy_keV <= 1.0) break;

                    bounce_pos = re_tr.exit_position;
                    bounce_dir = re_tr.exit_direction;
                    bounce_energy = re_tr.exit_energy_keV;
                }

                // Escaped secondary photons (fluorescence, brems, 511 keV
                // from crystal PP) through the water-bounce-crystal loop.
                constexpr int kMaxSecondaryBounces = 3;
                for (size_t si = 0; si < escaped_pool.size(); ++si) {
                    auto& sec = escaped_pool[si];
                    Eigen::Vector3d sp = sec.position, sd = sec.direction;
                    double se = sec.energy_keV;

                    for (int sb = 0; sb < kMaxSecondaryBounces; ++sb) {
                        auto reentry = source_geometry_.compute_marinelli_reentry(sp, sd);
                        if (!reentry.can_reenter) break;

                        auto sw = wall_transport(reentry.wall_path, sd, se);
                        if (sw == WallResult::Absorbed) break;
                        if (sw == WallResult::Scattered) {
                            sp = reentry.water_entry_pos;
                            continue;  // New direction — re-evaluate re-entry
                        }

                        auto wr = source_geometry_.transport_source_photon(
                            reentry.water_entry_pos, sd, se, rng);
                        if (!wr.survived) break;

                        auto s_segs = geometry_.trace_ray(wr.position, wr.direction);
                        if (s_segs.empty()) {
                            sp = wr.position; sd = wr.direction; se = wr.energy_keV;
                            continue;
                        }

                        const PathSegment* sf = nullptr;
                        for (const auto& seg : s_segs) {
                            if (seg.material) { sf = &seg; break; }
                        }
                        if (!sf) break;

                        Eigen::Vector3d s_entry = wr.position +
                            wr.direction * std::max(sf->t_start, 0.0);
                        auto s_tr = transport_photon(s_entry, wr.direction, wr.energy_keV,
                                                     geometry_, transport_config_, rng);

                        total_dep_scoring += s_tr.energy_deposited_scoring;
                        evt_dep_secondary += s_tr.energy_deposited_scoring;
                        evt_secondary_hits++;
                        if (s_tr.any_interaction_in_scoring) had_any_scoring = true;
                        break;  // Secondary deposited — done with this one
                    }
                }
            };

            // ------ Process primary photon if it survived source transport ------
            if (src_result.survived) {
                marinelli_detector_transport(
                    src_result.position, src_result.direction, src_result.energy_keV,
                    /*force_first=*/force_this_event);
            }

            // ------ Process 511 keV annihilation gammas from PP in source water ------
            // Each has already been transported through remaining source material;
            // now give it the full detector processing (miss-bounce, crystal
            // transport, post-crystal re-entry, escaped secondaries).
            double dep_before_pp = total_dep_scoring;
            const bool had_any_before_pp = had_any_scoring;
            for (auto& pp : source_pp) {
                tally.dbg_pp_secondary++;
                marinelli_detector_transport(pp.position, pp.direction, pp.energy_keV);
            }
            tally.dbg_dep_pp += (total_dep_scoring - dep_before_pp);
            if (had_any_scoring && !had_any_before_pp) {
                tally.dbg_pp_only_any++;
                tally.dbg_pp_only_any_w += src_air_w;
            }

            // ------ Process source electrons from Compton recoil in source material ------
            if (!src_result.source_electrons.empty()) {
                for (const auto& src_e : src_result.source_electrons) {
                    auto e_result = ElectronCsda::instance().deposited_in_scoring(
                        src_e.position, src_e.direction, src_e.energy_keV,
                        geometry_, rng,
                        transport_config_.disable_moliere,
                        transport_config_.disable_brems);

                    total_dep_scoring += e_result.deposited_scoring_keV;
                    if (e_result.deposited_scoring_keV > 0.0) had_any_scoring = true;

                    // Transport bremsstrahlung photons emitted by the electron
                    for (const auto& bp : e_result.brems_photons) {
                        auto br = transport_photon(
                            bp.position, bp.direction, bp.energy_keV,
                            geometry_, transport_config_, rng);
                        total_dep_scoring += br.energy_deposited_scoring;
                        if (br.any_interaction_in_scoring) had_any_scoring = true;
                    }
                }
            }

            // Accumulate per-event diagnostics into tally
            if (evt_had_initial_hit) tally.dbg_initial_hit++;
            if (evt_had_miss_bounce) tally.dbg_miss_bounce_hit++;
            tally.dbg_reentry_hit += evt_reentry_hits;
            tally.dbg_secondary_hit += evt_secondary_hits;
            tally.dbg_dep_initial += evt_dep_initial;
            tally.dbg_dep_reentry += evt_dep_reentry;
            tally.dbg_dep_secondary += evt_dep_secondary;

            if (had_any_scoring) {
                tally.num_any++;
                tally.sum_any_weights += src_air_w;
                tally.sum_any_w_sq += src_air_w * src_air_w;
                if (src_interacted) {
                    tally.dbg_src.any_w_s += src_air_w;
                    tally.dbg_src.any_w2_s += src_air_w * src_air_w;
                } else {
                    tally.dbg_src.any_w_u += src_air_w;
                    tally.dbg_src.any_w2_u += src_air_w * src_air_w;
                }
                if (electron_only_event) {
                    tally.dbg_src.e_only_any_w += src_air_w;
                    tally.dbg_src.e_only_any_w2 += src_air_w * src_air_w;
                }
            }

            if (std::abs(total_dep_scoring - energy_keV) < kFepTolerance) {
                tally.num_fep++;
                tally.sum_fep_weights += src_air_w;
                tally.sum_fep_w_sq += src_air_w * src_air_w;
                if (src_interacted) {
                    tally.dbg_src.fep_w_s += src_air_w;
                    tally.dbg_src.fep_w2_s += src_air_w * src_air_w;
                } else {
                    tally.dbg_src.fep_w_u += src_air_w;
                    tally.dbg_src.fep_w2_u += src_air_w * src_air_w;
                }
                if (electron_only_event)
                    tally.dbg_src.e_only_fep_w += src_air_w;
            }

            if (!tally.pulse_height_counts.empty() && total_dep_scoring > 0.0) {
                double dep = total_dep_scoring;
                auto it = std::upper_bound(energy_bin_edges.begin(),
                                           energy_bin_edges.end(),
                                           static_cast<float>(dep));
                if (it != energy_bin_edges.begin()) {
                    size_t bin = static_cast<size_t>(
                        it - energy_bin_edges.begin()) - 1;
                    if (bin < tally.pulse_height_counts.size()) {
                        tally.pulse_height_counts[bin] += src_air_w;
                        tally.pulse_height_counts_sq[bin] += src_air_w * src_air_w;
                    }
                }
            }
            continue;
        }

        // ================================================================
        // Standard path: direction sampling (surface, cone, or isotropic).
        // ================================================================
        double event_weight = 1.0;
        Eigen::Vector3d direction;

        if (active_bias.use_surface_sampling) {
            // Full-surface sampling: target a random point on the detector's
            // outer surface (front + side + back). Eliminates per-position
            // cone computation and Rodrigues rotation.
            auto [dir, w] = sample_surface_direction(src_pos, rng);
            if (w <= 0.0) {
                // Surface faces away from source — skip cheaply
                tally.sum_weights += 1.0;
                continue;
            }
            direction = dir;
            event_weight = w;
        } else if (use_cone_sampling_) {
            double cone_half_angle, cos_theta_max;
            if (per_position_cone) {
                cone_half_angle = compute_cone_half_angle(src_pos);
                cos_theta_max = std::cos(cone_half_angle);
                double solid_angle = kTwoPi * (1.0 - cos_theta_max);
                event_weight = solid_angle / kFourPi;
            } else {
                cone_half_angle = pre_cone_half_angle;
                cos_theta_max   = pre_cos_theta_max;
                event_weight    = pre_event_weight;
            }

            if (cone_half_angle < kPi - 0.01) {
                direction = sample_cone_direction(src_pos, cos_theta_max, rng);
            } else {
                double cos_theta = 2.0 * uniform(rng) - 1.0;
                double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                double phi = kTwoPi * uniform(rng);
                direction = Eigen::Vector3d(
                    sin_theta * std::cos(phi),
                    sin_theta * std::sin(phi),
                    cos_theta);
                event_weight = 1.0;
            }
        } else {
            double cos_theta = 2.0 * uniform(rng) - 1.0;
            double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
            double phi = kTwoPi * uniform(rng);
            direction = Eigen::Vector3d(
                sin_theta * std::cos(phi),
                sin_theta * std::sin(phi),
                cos_theta);
        }

        tally.sum_weights += 1.0;

        geometry_.trace_ray(src_pos, direction, event_segments);

        if (event_segments.empty()) {
            continue;
        }

        const PathSegment* first_seg = nullptr;
        for (const auto& seg : event_segments) {
            if (seg.material) {
                first_seg = &seg;
                break;
            }
        }

        if (!first_seg) {
            continue;
        }

        double total_weight = event_weight * position_weight;

        // FEP-only mode with source effects: stochastic Rayleigh
        if (source_geometry_.has_source_effects()) {
            auto src_result = source_geometry_.compute_transmission_fep_only(
                src_pos, direction, energy_keV * 1e-3, rng);
            total_weight *= src_result.weight;

            if (src_result.direction != direction) {
                direction = src_result.direction;
                geometry_.trace_ray(src_pos, direction, event_segments);
                if (event_segments.empty()) continue;
                first_seg = nullptr;
                for (const auto& seg : event_segments) {
                    if (seg.material) { first_seg = &seg; break; }
                }
                if (!first_seg) continue;
            }

            // Air attenuation: distance from source exit to detector entry
            if (mu_air_no_rs > 0.0) {
                Eigen::Vector3d det_entry = src_pos + direction * std::max(first_seg->t_start, 0.0);
                double air_dist = (det_entry - src_result.exit_position).norm();
                total_weight *= std::exp(-mu_air_no_rs * air_dist);
            }
        } else if (mu_air_no_rs > 0.0) {
            // No source effects: air gap = distance from source to detector entry
            double air_dist = std::max(first_seg->t_start, 0.0);
            total_weight *= std::exp(-mu_air_no_rs * air_dist);
        }

        Eigen::Vector3d entry_point = src_pos + direction * std::max(first_seg->t_start, 0.0);

        if (fep_only_mode_) {
            // FEP-only: stochastic Rayleigh in detector non-scoring segments,
            // then full MC in scoring volume with boundary-exit kill.
            auto rayleigh_result = compute_transmission_stochastic_rayleigh(
                src_pos, direction, energy_keV, geometry_, rng);

            if (!rayleigh_result.hit_scoring) {
                continue;
            }
            total_weight *= rayleigh_result.weight;

            TransportConfig fep_config = transport_config_;
            fep_config.fep_only_mode = true;
            fep_config.force_first_interaction = force_interaction;
            auto transport_result = transport_photon(
                rayleigh_result.position, rayleigh_result.direction, energy_keV,
                geometry_, fep_config, rng);
            total_weight *= transport_result.weight;
            if (transport_result.forced_absorption) tally.num_forced_absorption++;

            if (std::abs(transport_result.energy_deposited_scoring - energy_keV) < kFepTolerance) {
                tally.num_fep++;
                tally.sum_fep_weights += total_weight;
                tally.sum_fep_w_sq += total_weight * total_weight;
            }
        } else {
            // Full mode without source effects
            auto transport_result = transport_photon(
                entry_point, direction, energy_keV, geometry_,
                force_interaction ? forced_transport_config : transport_config_,
                rng);
            total_weight *= transport_result.weight;
            if (transport_result.forced_absorption) tally.num_forced_absorption++;

            if (transport_result.any_interaction_in_scoring) {
                tally.num_any++;
                tally.sum_any_weights += total_weight;
                tally.sum_any_w_sq += total_weight * total_weight;
            }

            if (std::abs(transport_result.energy_deposited_scoring - energy_keV) < kFepTolerance) {
                tally.num_fep++;
                tally.sum_fep_weights += total_weight;
                tally.sum_fep_w_sq += total_weight * total_weight;
            }

            if (!tally.pulse_height_counts.empty() &&
                transport_result.energy_deposited_scoring > 0.0) {
                double dep = transport_result.energy_deposited_scoring;
                auto it = std::upper_bound(energy_bin_edges.begin(),
                                           energy_bin_edges.end(),
                                           static_cast<float>(dep));
                if (it != energy_bin_edges.begin()) {
                    size_t bin = static_cast<size_t>(it - energy_bin_edges.begin()) - 1;
                    if (bin < tally.pulse_height_counts.size()) {
                        tally.pulse_height_counts[bin] += total_weight;
                        tally.pulse_height_counts_sq[bin] += total_weight * total_weight;
                    }
                }
            }
        }
    }

    return tally;
}

// ===================== Correlated-cascade MC (compute_cascade) ================

namespace {

/// Sample an isotropic unit direction.
inline Eigen::Vector3d sample_isotropic_dir(
        std::uniform_real_distribution<double>& u, std::mt19937_64& rng) {
    const double ct = 2.0 * u(rng) - 1.0;
    const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
    const double ph = kTwoPi * u(rng);
    return Eigen::Vector3d(st * std::cos(ph), st * std::sin(ph), ct);
}

/// Find an angular-correlation link between members a and b (either direction).
/// Returns true and fills a2/a4 if a correlated link exists.
inline bool cascade_corr_link(const DecayCascade& dc, std::size_t a, std::size_t b,
                              double& a2, double& a4) {
    for (const auto& c : dc.members[a].coincident)
        if (c.partner == b && c.has_correlation) { a2 = c.a2; a4 = c.a4; return true; }
    for (const auto& c : dc.members[b].coincident)
        if (c.partner == a && c.has_correlation) { a2 = c.a2; a4 = c.a4; return true; }
    return false;
}

/// Sample the x-ray emitted when a vacancy in shell `f` is filled: returns the
/// line energy (keV), or 0 for an Auger transition (no x-ray). Works for both
/// FluorescenceData (K) and LSubshellFluor (one L subshell) -- same member layout.
template <class FD>
inline double sample_vacancy_xray(const FD* f,
                                  std::uniform_real_distribution<double>& u,
                                  std::mt19937_64& rng) {
    if (!f || f->num_lines == 0) return 0.0;
    if (u(rng) >= f->fluorescence_yield) return 0.0;  // Auger, no x-ray
    double xi = u(rng), cum = 0.0;
    for (int li = 0; li < f->num_lines; ++li) {
        cum += f->line_probability[li];
        if (xi <= cum) return f->line_energy_keV[li];
    }
    return f->line_energy_keV[0];
}

/// P(member m emitted | member pi emitted) within one decay branch, from the
/// pairwise coincidence data. A direct link pi->m gives P(m|pi); a reverse link
/// m->pi gives P(pi|m), inverted by Bayes using the marginal intensities (this
/// reproduces the transition-group result P(xray|gamma) = intensity_xray); with
/// no link it falls back to the marginal (conditional-independence given pi).
double cascade_partner_prob(const DecayCascade& dc, std::size_t pi, std::size_t m) {
    const CascadeMember& P = dc.members[pi];
    const CascadeMember& M = dc.members[m];
    for (const auto& c : P.coincident)
        if (c.partner == m)
            return std::clamp(c.prob, 0.0, 1.0);
    for (const auto& c : M.coincident)
        if (c.partner == pi) {
            if (P.intensity <= 0.0) return 0.0;
            return std::clamp(c.prob * M.intensity / P.intensity, 0.0, 1.0);
        }
    return std::clamp(M.intensity, 0.0, 1.0);
}

/// Sample the coherent, marginal-preserving fallback forest. The forest is
/// precomputed once per branch; one uniform is consumed per member.
void cascade_sample_realization(
                                std::uniform_real_distribution<double>& u,
                                std::mt19937_64& rng,
                                const std::vector<CascadeFallbackNode>& forest,
                                int forced_annihilation_state,
                                std::vector<char>& emitted) {
    sample_cascade_fallback_forest(
        forest, forced_annihilation_state, [&]() { return u(rng); }, emitted);
}

// --- beta+ / annihilation helpers (cascade 511 summing) --------------------

/// Sample a positron kinetic energy (keV) from the allowed beta+ spectrum with
/// endpoint `endpoint_keV`: N(T) ~ p*E_tot*(E0-T)^2 (the Fermi/Coulomb function
/// is omitted -- a few-percent shape tweak, negligible for the annihilation site
/// / escape). Rejection against a numerically-found envelope max.
double sample_beta_plus_KE(double endpoint_keV,
                           std::uniform_real_distribution<double>& u,
                           std::mt19937_64& rng) {
    constexpr double mc2 = 510.998950;
    const double E0 = endpoint_keV;
    if (E0 <= 0.0) return 0.0;
    auto shape = [&](double T) {
        const double p = std::sqrt(T * (T + 2.0 * mc2));   // momentum*c
        return p * (T + mc2) * (E0 - T) * (E0 - T);
    };
    double fmax = 0.0;
    for (int i = 1; i < 40; ++i) fmax = std::max(fmax, shape(E0 * i / 40.0));
    if (fmax <= 0.0) return 0.5 * E0;
    for (int tries = 0; tries < 100; ++tries) {
        const double T = E0 * u(rng);
        if (u(rng) * fmax <= shape(T)) return T;
    }
    return 0.4 * E0;
}

/// Total positron-electron annihilation cross-section per electron (Heitler), in
/// cm^2, for a positron of kinetic energy T_keV (Geant4 G4eplusAnnihilation form).
double annih_xsec_per_electron(double T_keV) {
    constexpr double mc2 = 510.998950;
    const double gamma = 1.0 + T_keV / mc2;
    if (gamma <= 1.0 + 1.0e-6) return 0.0;
    const double g2 = gamma * gamma;
    const double s = std::sqrt(g2 - 1.0);
    constexpr double pi_re2 = 2.4946e-25;  // pi * r_e^2 (cm^2)
    return (pi_re2 / (gamma + 1.0)) *
           (((g2 + 4.0 * gamma + 1.0) / (g2 - 1.0)) * std::log(gamma + s)
            - (gamma + 3.0) / s);
}

/// Probability a positron of initial KE T_keV annihilates IN FLIGHT before
/// stopping in `mat`: the macroscopic annihilation rate integrated over the
/// slowing-down against the collision stopping power, P = int n_e sigma_a / S dE.
double inflight_annih_prob(const Material& mat, double T_keV) {
    if (T_keV <= 1.0) return 0.0;
    constexpr double NA = 6.02214076e23;
    const double rho = mat.density();
    double ne = 0.0;  // electron density (1/cm^3)
    for (const auto& c : mat.composition())
        ne += rho * NA * c.mass_fraction * c.Z / ElectronCsda::atomic_weight(c.Z);
    auto S_keV_cm = [&](double E_keV) {        // positron collision stopping (keV/cm)
        double s = 0.0;
        for (const auto& c : mat.composition())
            s += c.mass_fraction *
                 ElectronCsda::stopping_power_MeV_cm2_g(
                     c.Z, ElectronCsda::atomic_weight(c.Z), E_keV,
                     /*is_positron=*/true);
        return s * rho * 1000.0;
    };
    constexpr int N = 20;
    double P = 0.0;
    for (int i = 0; i < N; ++i) {
        const double E = T_keV * (i + 0.5) / N;
        const double S = S_keV_cm(E);
        if (S > 0.0) P += (ne * annih_xsec_per_electron(E) / S) * (T_keV / N);
    }
    return std::min(P, 0.5);
}

// Absolute daughter-path probabilities.  LevelDag::pass()/n_subset_joint()
// average over LevelScheme::feeding after normalising it, so they are only the
// right primitive for metadata-free legacy cascades.  A SandiaDecay
// WeakOutcomeLaw already gives an absolute, branch-conditional probability for
// each EC/beta+ outcome and (when matched) its exact fed level.
std::vector<double> weak_unknown_start_weights(const DecayCascade& dc,
                                               const LevelDag& dag)
{
    std::vector<double> w(dag.n, 0.0), known(dag.n, 0.0);
    double unknown_mass = 0.0;
    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
        if (o.fed_level >= 0 && o.fed_level < dag.n)
            known[o.fed_level] += o.selected_mass;
        else
            unknown_mass += o.selected_mass;
    }
    if (!(unknown_mass > 0.0)) return w;
    double residual = 0.0;
    for (int l = 0; l < dag.n; ++l) {
        w[l] = std::max(0.0, dc.level_scheme.levels[l].feeding - known[l]);
        residual += w[l];
    }
    // The residual graph can itself be inconsistent with the categorical mass.
    // Preserve its relative level distribution, but never assign more than the
    // available unknown probability; the unassigned remainder lands in a
    // terminal/direct state and emits no scheme transition.
    const double scale = residual > unknown_mass ? unknown_mass / residual : 1.0;
    for (double& x : w) x *= scale / unknown_mass;
    return w;  // conditional start mass per unknown outcome; sum <= 1
}

double weak_outcome_reach(const DecayCascade& dc, const LevelDag& dag,
                          const WeakOutcome& o,
                          const std::vector<int>& subset)
{
    // Selecting a categorical outcome is certain even when its fed level is
    // unknown or terminal.  Residual graph weights below describe only the
    // probability of entering a transition path, not the outcome's existence
    // (notably its EC vacancy).
    if (subset.empty()) return 1.0;
    if (o.fed_level >= 0 && o.fed_level < dag.n)
        return dag.subset_from_level(o.fed_level, subset);
    const std::vector<double> w = weak_unknown_start_weights(dc, dag);
    double p = 0.0;
    for (int l = 0; l < dag.n; ++l)
        if (w[l] > 0.0) p += w[l] * dag.subset_from_level(l, subset);
    return p;
}

double weak_subset_probability(const DecayCascade& dc, const LevelDag& dag,
                               const std::vector<int>& subset,
                               WeakOutcomeKind* only_kind = nullptr)
{
    if (!dag.valid) return 0.0;
    if (!dc.weak_outcome_law.usable())
        return dc.level_scheme.entry_probability * dag.n_subset_joint(subset);

    double p = 0.0;
    for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
        if (only_kind && o.kind != *only_kind) continue;
        const double reach = weak_outcome_reach(dc, dag, o, subset);
        p += o.selected_mass * reach;
    }
    return std::clamp(p, 0.0, 1.0);
}

double weak_pass_probability(const DecayCascade& dc, const LevelDag& dag, int t,
                             WeakOutcomeKind* only_kind = nullptr)
{
    return weak_subset_probability(dc, dag, std::vector<int>{t}, only_kind);
}

// P(outcome o | all transitions in subset were passed).  Photon/IC outcome
// factors cancel when they are common to every weak mode, which is the case for
// conditioning on emitted nuclear gammas.  Unknown-level outcomes retain the
// audited legacy feeding distribution.
std::vector<double> weak_outcome_posterior(const DecayCascade& dc,
                                           const LevelDag& dag,
                                           const std::vector<int>& subset)
{
    std::vector<double> post;
    if (!dc.weak_outcome_law.usable() || !dag.valid) return post;
    post.resize(dc.weak_outcome_law.outcomes.size(), 0.0);
    double den = 0.0;
    for (std::size_t k = 0; k < post.size(); ++k) {
        const WeakOutcome& o = dc.weak_outcome_law.outcomes[k];
        const double reach = weak_outcome_reach(dc, dag, o, subset);
        post[k] = o.selected_mass * reach;
        den += post[k];
    }
    if (den > 0.0)
        for (double& x : post) x /= den;
    else
        post.clear();
    return post;
}

int residual_transition_of(const DecayCascade& dc, int member)
{
    for (int r = 0; r < static_cast<int>(dc.residual_transitions.size()); ++r)
        if (dc.residual_transitions[static_cast<std::size_t>(r)].gamma_member == member)
            return r;
    return -1;
}

// P(member m emitted | the conditioned transition subset, weak outcome k).
// The outcome is drawn once by the Conditional estimator, so every nuclear
// partner and conversion-vacancy draw must use the same row rather than the
// outcome-averaged marginal.  Otherwise mutually exclusive EC/beta+ feed paths
// are silently recombined in a single simulated history.
std::vector<std::vector<double>> weak_outcome_partner_probabilities(
    const DecayCascade& dc, const LevelDag& dag,
    const std::vector<int>& conditioned_subset, int reference_member = -1)
{
    std::vector<std::vector<double>> rows;
    if (!dc.weak_outcome_law.usable() || !dag.valid) return rows;
    rows.assign(dc.weak_outcome_law.outcomes.size(),
                std::vector<double>(dc.members.size(), 0.0));
    for (std::size_t k = 0; k < rows.size(); ++k) {
        const WeakOutcome& o = dc.weak_outcome_law.outcomes[k];
        const double den = weak_outcome_reach(dc, dag, o, conditioned_subset);
        if (!(den > 0.0)) continue;
        for (std::size_t m = 0; m < dc.members.size(); ++m) {
            if (dc.members[m].type == CascadeParticleType::Annih511) {
                rows[k][m] = o.kind == WeakOutcomeKind::Positron ? 1.0 : 0.0;
                continue;
            }
            double p = 0.0;
            bool matched = false;
            for (int t = 0; t < static_cast<int>(dag.ts.size()); ++t) {
                if (dag.ts[t].gamma_member != static_cast<int>(m)) continue;
                matched = true;
                std::vector<int> joint = conditioned_subset;
                joint.push_back(t);
                p += weak_outcome_reach(dc, dag, o, joint) / den *
                     dag.ts[t].p_gamma;
            }
            if (matched)
                rows[k][m] = std::clamp(p, 0.0, 1.0);
            else if (const int r = residual_transition_of(
                         dc, static_cast<int>(m)); r >= 0)
                rows[k][m] = std::clamp(
                    dc.residual_transitions[static_cast<std::size_t>(r)].p_gamma,
                    0.0, 1.0);
            else
                // Flat legacy x-ray members have no fed-level metadata.  They
                // retain their historical pairwise conditional; atomic
                // vacancies constructed by this change never use this path.
                rows[k][m] = reference_member >= 0
                    ? cascade_partner_prob(dc,
                        static_cast<std::size_t>(reference_member), m)
                    : std::clamp(dc.members[m].intensity, 0.0, 1.0);
        }
    }
    return rows;
}

} // namespace

double EfficiencyCalculator::transport_into_detector(
    const Eigen::Vector3d& start, const Eigen::Vector3d& direction,
    double energy_keV, std::vector<PathSegment>& seg_buffer,
    std::mt19937_64& rng) const
{
    geometry_.trace_ray(start, direction, seg_buffer);
    const PathSegment* first = nullptr;
    for (const auto& s : seg_buffer)
        if (s.material) { first = &s; break; }
    if (!first)
        return 0.0;  // ray misses the detector geometry
    const Eigen::Vector3d entry = start + direction * std::max(first->t_start, 0.0);
    const TransportResult tr = transport_photon(entry, direction, energy_keV,
                                                geometry_, transport_config_, rng);
    return tr.energy_deposited_scoring;
}

double EfficiencyCalculator::transport_cascade_member(
    const Eigen::Vector3d& vertex, const Eigen::Vector3d& direction,
    double energy_keV, std::vector<PathSegment>& seg_buffer,
    std::mt19937_64& rng) const
{
    // No source material/shield: emit straight into the detector.
    if (!source_geometry_.has_source_effects())
        return transport_into_detector(vertex, direction, energy_keV, seg_buffer, rng);

    // Transport through the source material + shielding first (self-attenuation
    // and Compton/Rayleigh scatter), then sum the survivor and any source-
    // generated secondary photons (511s from pair production, bremsstrahlung).
    const SourceGeometry::SourceFullTransportResult src =
        source_geometry_.transport_source_photon(vertex, direction, energy_keV, rng);
    double dep = 0.0;
    if (src.survived)
        dep += transport_into_detector(src.position, src.direction,
                                       src.energy_keV, seg_buffer, rng);
    for (const SourceGeometry::SourceSecondaryPhoton& s : src.secondaries)
        dep += transport_into_detector(s.position, s.direction,
                                       s.energy_keV, seg_buffer, rng);
    return dep;
}

// LevelDag (the shared visit/reach DAG DP) now lives in cascade/LevelDag.h so the
// analytic summing path (compute_cascade_analytic) can reuse the exact combinatorics
// the MC estimators below rely on. Included at the top of this file.

std::vector<double> EfficiencyCalculator::cascade_level_pmate(
    const DecayCascade& dc, std::size_t primary_index)
{
    const LevelDag dag(dc.level_scheme);
    if (!dag.valid) return {};

    // A valid-scheme primary can be either a matched DAG transition or one of
    // the categorical unmatched residuals.  Only flat x-rays/legacy members
    // fall back to the pairwise construction.
    const int ip = dag.transition_of(static_cast<int>(primary_index));
    const int ir = residual_transition_of(dc, static_cast<int>(primary_index));
    if (ip < 0 && ir < 0) return {};
    const double pass_p = ip >= 0 ? weak_pass_probability(dc, dag, ip) : 1.0;
    if (pass_p <= 0.0) return {};
    std::vector<double> post = ip >= 0
        ? weak_outcome_posterior(dc, dag, {ip}) : std::vector<double>{};
    if (ip < 0 && dc.weak_outcome_law.usable()) {
        post.reserve(dc.weak_outcome_law.outcomes.size());
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
            post.push_back(o.selected_mass);
    }
    double post_pos = 0.0;
    for (std::size_t k = 0; k < post.size(); ++k)
        if (dc.weak_outcome_law.outcomes[k].kind == WeakOutcomeKind::Positron)
            post_pos += post[k];

    // pmate[m] = P(member m gamma emitted | primary gamma emitted). The primary's
    // own p_gamma cancels between the joint's numerator and pass_p, so
    // conditioning on "primary gamma emitted" (vs. just "transition ip passed")
    // is exact.
    std::vector<double> pmate(dc.members.size(), 0.0);
    for (std::size_t m = 0; m < dc.members.size(); ++m) {
        if (m == primary_index) continue;
        double num = 0.0;
        bool matched = false;
        for (int i = 0; i < static_cast<int>(dag.ts.size()); ++i)
            if (dag.ts[i].gamma_member == static_cast<int>(m)) {
                matched = true;
                num += (ip >= 0
                            ? weak_subset_probability(dc, dag, {ip, i}) / pass_p
                            : weak_subset_probability(dc, dag, {i}))
                     * dag.ts[i].p_gamma;
            }
        // x-ray / annihilation members are not level-scheme gammas; they are
        // emitted independently of the gamma path (as FullRealization treats
        // them), so keep their pairwise/marginal probability.
        if (matched)
            pmate[m] = std::clamp(num, 0.0, 1.0);
        else if (const int r = residual_transition_of(
                     dc, static_cast<int>(m)); r >= 0)
            pmate[m] = std::clamp(
                dc.residual_transitions[static_cast<std::size_t>(r)].p_gamma,
                0.0, 1.0);
        else if (dc.members[m].type == CascadeParticleType::Annih511 && !post.empty())
            pmate[m] = std::clamp(post_pos, 0.0, 1.0);
        else
            pmate[m] = cascade_partner_prob(dc, primary_index, m);
    }
    return pmate;
}

std::vector<EfficiencyCalculator::CascadeVacancyDraw>
EfficiencyCalculator::cascade_level_vacancies(const DecayCascade& dc,
                                              std::size_t primary_index)
{
    const LevelDag dag(dc.level_scheme);
    if (!dag.valid) return {};
    const int ip = dag.transition_of(static_cast<int>(primary_index));
    const int ir = residual_transition_of(dc, static_cast<int>(primary_index));
    if (ip < 0 && ir < 0) return {};
    const double pass_p = ip >= 0 ? weak_pass_probability(dc, dag, ip) : 1.0;
    if (pass_p <= 0.0) return {};

    std::vector<CascadeVacancyDraw> out;
    // Electron-capture vacancy at the directly-fed level, coincident with the
    // primary: P(EC-K | primary emitted) = Σ_L feed[L]·feed_ecK[L]·R_ip[L] / pass_p
    // (the primary's own emit/convert factor cancels, so this is independent of
    // whether the primary converts). Same for EC-L.
    double ecK = 0.0, ecL = 0.0;
    std::vector<double> post = ip >= 0
        ? weak_outcome_posterior(dc, dag, {ip}) : std::vector<double>{};
    if (ip < 0 && dc.weak_outcome_law.usable()) {
        post.reserve(dc.weak_outcome_law.outcomes.size());
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
            post.push_back(o.selected_mass);
    }
    if (!post.empty()) {
        for (std::size_t k = 0; k < post.size(); ++k) {
            const WeakOutcome& o = dc.weak_outcome_law.outcomes[k];
            if (o.kind != WeakOutcomeKind::ElectronCapture) continue;
            ecK += post[k] * o.ec_K;
            ecL += post[k] * o.ec_L;
        }
    } else {
        if (ip >= 0) {
            for (int L = 0; L < dag.n; ++L) {
                const double f = dag.feed[L] * dag.R[ip][L];
                ecK += f * dc.level_scheme.levels[L].feed_ecK;
                ecL += f * dc.level_scheme.levels[L].feed_ecL;
            }
            ecK /= dag.pass(ip);
            ecL /= dag.pass(ip);
        } else {
            const double ep = dc.level_scheme.entry_probability;
            for (int L = 0; L < dag.n; ++L) {
                const double f = ep * dag.feed[L];
                ecK += f * dc.level_scheme.levels[L].feed_ecK;
                ecL += f * dc.level_scheme.levels[L].feed_ecL;
            }
            if (dag.n > 0) {
                ecK += (1.0 - ep) * dc.level_scheme.levels[0].feed_ecK;
                ecL += (1.0 - ep) * dc.level_scheme.levels[0].feed_ecL;
            }
        }
    }
    if (ecK > 1e-9) out.push_back({std::clamp(ecK, 0.0, 1.0), false, -1,
                                   0.0, 0, -1, {}});
    if (ecL > 1e-9) out.push_back({std::clamp(ecL, 0.0, 1.0), true, -1,
                                   0.0, 0, -1, {}});

    // Internal-conversion vacancies from the OTHER transitions on the cascade
    // path (a transition either emits its gamma or converts). Conditioning on the
    // primary gamma emitted: P(IC-K from j | primary) = joint(ip,j)·p_icK[j]/pass_p
    // (j == ip excluded: the primary emitted its gamma, so it did not convert).
    for (int j = 0; j < static_cast<int>(dag.ts.size()); ++j) {
        if (j == ip) continue;
        const double base = ip >= 0
            ? weak_subset_probability(dc, dag, {ip, j}) / pass_p
            : weak_subset_probability(dc, dag, {j});
        if (base <= 1e-12) continue;
        // IC transition energy, taken from the transition itself: an E0 has no
        // photon member, and recovering the energy through gamma_member gave 0,
        // which the `trans_keV > 0` guard downstream then read as "no electron".
        const double e_tr = dag.ts[j].gamma_keV;
        const double icK = base * dag.ts[j].p_icK;
        auto outcome_shell_prob = [&](double shell_p) {
            std::vector<double> p;
            if (!post.empty()) {
                p.resize(dc.weak_outcome_law.outcomes.size(), 0.0);
                for (std::size_t k = 0; k < p.size(); ++k) {
                    const WeakOutcome& o = dc.weak_outcome_law.outcomes[k];
                    const std::vector<int> conditioned = ip >= 0
                        ? std::vector<int>{ip} : std::vector<int>{};
                    const double den = weak_outcome_reach(
                        dc, dag, o, conditioned);
                    std::vector<int> joint = conditioned;
                    joint.push_back(j);
                    p[k] = den > 0.0
                        ? std::clamp(weak_outcome_reach(dc, dag, o, joint) /
                                     den * shell_p, 0.0, 1.0)
                        : 0.0;
                }
            }
            return p;
        };
        if (icK > 1e-9) {
            CascadeVacancyDraw v{std::clamp(icK, 0.0, 1.0), false, -1,
                                 e_tr, j + 1, dag.ts[j].gamma_member, {}};
            v.weak_outcome_prob = outcome_shell_prob(dag.ts[j].p_icK);
            out.push_back(std::move(v));
        }
        for (int s = 0; s < 3; ++s) {
            const double icL = base * dag.ts[j].p_icL[s];
            if (icL > 1e-9) {
                CascadeVacancyDraw v{std::clamp(icL, 0.0, 1.0), true, s,
                                     e_tr, j + 1, dag.ts[j].gamma_member, {}};
                v.weak_outcome_prob = outcome_shell_prob(dag.ts[j].p_icL[s]);
                out.push_back(std::move(v));
            }
        }
        const double unresolved = base * dag.ts[j].p_icUnresolved;
        if (unresolved > 1e-9) {
            CascadeVacancyDraw v{std::clamp(unresolved, 0.0, 1.0), false, -1,
                                 e_tr, j + 1, dag.ts[j].gamma_member, {}};
            v.weak_outcome_prob =
                outcome_shell_prob(dag.ts[j].p_icUnresolved);
            v.produces_vacancy = false;
            out.push_back(std::move(v));
        }
    }

    // Independent categorical residual conversions.  Probabilities remain
    // absolute here; Conditional's common group sampler divides by the
    // complementary residual-gamma probability when that member did not fire.
    // When the conditioned primary itself is residual its gamma outcome is
    // forced, so its own conversion alternatives are omitted entirely.
    const int residual_group_base = static_cast<int>(dag.ts.size()) + 1;
    for (int rix = 0; rix < static_cast<int>(dc.residual_transitions.size()); ++rix) {
        if (rix == ir) continue;
        const ResidualTransition& r =
            dc.residual_transitions[static_cast<std::size_t>(rix)];
        const int group = residual_group_base + rix;
        auto append = [&](double p, bool is_L, int subshell,
                          bool produces_vacancy) {
            if (p <= 1e-9) return;
            CascadeVacancyDraw v{std::clamp(p, 0.0, 1.0), is_L, subshell,
                                 r.transition_energy_keV, group,
                                 r.gamma_member, {}};
            if (!post.empty())
                v.weak_outcome_prob.assign(post.size(), std::clamp(p, 0.0, 1.0));
            v.produces_vacancy = produces_vacancy;
            out.push_back(std::move(v));
        };
        append(r.p_icK, false, -1, true);
        append(r.p_icL1, true, 0, true);
        append(r.p_icL2, true, 1, true);
        append(r.p_icL3, true, 2, true);
        append(r.p_icOuter + r.p_icUnresolved, false, -1, false);
    }
    return out;
}

std::vector<EfficiencyCalculator::CascadeSumPairChannel>
EfficiencyCalculator::cascade_sum_pair_channels(
    const std::vector<DecayCascade>& cascades, std::size_t primary_branch,
    std::size_t primary_index, double peak_keV, double tol_keV)
{
    validate_cascade_member_references(
        cascades, "EfficiencyCalculator::cascade_sum_pair_channels");
    if (primary_branch >= cascades.size()) return {};
    const DecayCascade& dcA = cascades[primary_branch];
    const LevelDag dagA(dcA.level_scheme);
    if (!dagA.valid) return {};
    const int ipA = dagA.transition_of(static_cast<int>(primary_index));
    const int irA = residual_transition_of(dcA, static_cast<int>(primary_index));
    if (ipA < 0 && irA < 0) return {};
    const double passA = ipA >= 0
        ? weak_pass_probability(dcA, dagA, ipA) : 1.0;
    const double pg_p = ipA >= 0 ? dagA.ts[ipA].p_gamma
        : dcA.residual_transitions[static_cast<std::size_t>(irA)].p_gamma;
    const double bwA = dcA.branch_weight;
    const double denom = bwA * passA * pg_p;  // P(primary gamma emitted per decay)
    if (denom <= 0.0) return {};

    // Both pair photons must exceed 2·tol (and >= 5 keV) so that when the
    // cone-biased photon a misses the detector the window is unreachable -- this
    // is what makes cone-biasing photon a unbiased (see the block comment above
    // compute_cascade). Sub-5-keV sum feeders are also the x-ray regime (handled
    // by the vacancy model, not gamma pairs).
    const double min_e = std::max(2.0 * tol_keV, 5.0);

    std::vector<CascadeSumPairChannel> out;
    for (std::size_t B = 0; B < cascades.size(); ++B) {
        const DecayCascade& dcB = cascades[B];
        const LevelDag dagB(dcB.level_scheme);
        if (!dagB.valid) continue;
        const double bwB = dcB.branch_weight;
        const int nts = static_cast<int>(dagB.ts.size());
        for (int i = 0; i < nts; ++i) {
            if (dagB.ts[i].p_gamma <= 0.0) continue;
            const int mi = dagB.ts[i].gamma_member;
            if (mi < 0 || static_cast<std::size_t>(mi) >= dcB.members.size())
                continue;
            if (B == primary_branch && mi == static_cast<int>(primary_index)) continue;
            const double Ei = dcB.members[mi].energy_keV;
            for (int j = i + 1; j < nts; ++j) {
                if (dagB.ts[j].p_gamma <= 0.0) continue;
                const int mj = dagB.ts[j].gamma_member;
                if (mj < 0 || static_cast<std::size_t>(mj) >= dcB.members.size())
                    continue;
                if (B == primary_branch && mj == static_cast<int>(primary_index)) continue;
                const double Ej = dcB.members[mj].energy_keV;
                if (std::abs(Ei + Ej - peak_keV) > tol_keV) continue;
                if (std::min(Ei, Ej) < min_e) continue;
                const int excluded_gamma_t = B == primary_branch ? ipA : -1;
                const bool excluded_residual =
                    B == primary_branch && irA >= 0;
                auto event_probability = [&](std::vector<int> subset,
                                             int ic_transition = -1) {
                    double p = weak_subset_probability(dcB, dagB, subset);
                    if (excluded_gamma_t >= 0 &&
                        ic_transition != excluded_gamma_t) {
                        subset.push_back(excluded_gamma_t);
                        p -= weak_subset_probability(dcB, dagB, subset) * pg_p;
                    }
                    if (excluded_residual)
                        p *= 1.0 - pg_p;
                    return std::max(0.0, p);
                };
                auto outcome_event_probability = [&](const WeakOutcome& o,
                                                     std::vector<int> subset,
                                                     int ic_transition = -1) {
                    double p = weak_outcome_reach(dcB, dagB, o, subset);
                    if (excluded_gamma_t >= 0 &&
                        ic_transition != excluded_gamma_t) {
                        subset.push_back(excluded_gamma_t);
                        p -= weak_outcome_reach(dcB, dagB, o, subset) * pg_p;
                    }
                    if (excluded_residual)
                        p *= 1.0 - pg_p;
                    return std::max(0.0, p);
                };
                const double pj = event_probability({i, j});
                if (pj <= 0.0) continue;
                const double pgi = dagB.ts[i].p_gamma, pgj = dagB.ts[j].p_gamma;
                // w_k per primary emission (Eq. in the block comment): P(a&b gammas
                // emitted per decay) minus the triple overlap the primary
                // conditional stream already samples, over P(primary emitted).
                const double num = bwB * pj * pgi * pgj;
                const double w = num / denom;
                if (w < 1e-9) continue;

                CascadeSumPairChannel ch;
                ch.e_a_keV = Ei; ch.e_b_keV = Ej;
                ch.member_a = mi; ch.member_b = mj;
                ch.branch = B;
                ch.weight = w;
                if (dcB.weak_outcome_law.usable()) {
                    ch.weak_outcome_posterior.assign(
                        dcB.weak_outcome_law.outcomes.size(), 0.0);
                    double ps = 0.0;
                    for (std::size_t k = 0;
                         k < dcB.weak_outcome_law.outcomes.size(); ++k) {
                        const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                        ch.weak_outcome_posterior[k] = o.selected_mass *
                            outcome_event_probability(o, {i, j});
                        ps += ch.weak_outcome_posterior[k];
                    }
                    if (ps > 0.0)
                        for (double& p : ch.weak_outcome_posterior) p /= ps;
                    else
                        ch.weak_outcome_posterior.clear();
                    ch.weak_outcome_partner_prob =
                        weak_outcome_partner_probabilities(dcB, dagB, {i, j}, mi);
                }

                // Pair-conditioned partner emission probabilities P(m | a&b).
                ch.partner_prob.assign(dcB.members.size(), 0.0);
                for (std::size_t m = 0; m < dcB.members.size(); ++m) {
                    if (static_cast<int>(m) == mi || static_cast<int>(m) == mj) continue;
                    const int tm = dagB.transition_of(static_cast<int>(m));
                    if (tm == excluded_gamma_t) {
                        ch.partner_prob[m] = 0.0;
                        for (auto& row : ch.weak_outcome_partner_prob)
                            if (m < row.size()) row[m] = 0.0;
                    } else if (tm >= 0) {
                        ch.partner_prob[m] = std::clamp(
                            event_probability({i, j, tm})
                                * dagB.ts[tm].p_gamma / pj,
                            0.0, 1.0);
                        for (std::size_t k = 0;
                             k < ch.weak_outcome_partner_prob.size(); ++k) {
                            const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                            const double den_o = outcome_event_probability(o, {i, j});
                            ch.weak_outcome_partner_prob[k][m] = den_o > 0.0
                                ? std::clamp(outcome_event_probability(
                                      o, {i, j, tm}) * dagB.ts[tm].p_gamma / den_o,
                                      0.0, 1.0)
                                : 0.0;
                        }
                    } else if (dcB.members[m].type == CascadeParticleType::Annih511 &&
                             dcB.weak_outcome_law.usable()) {
                        double pp = 0.0;
                        for (std::size_t k = 0;
                             k < ch.weak_outcome_posterior.size(); ++k)
                            if (dcB.weak_outcome_law.outcomes[k].kind ==
                                WeakOutcomeKind::Positron)
                                pp += ch.weak_outcome_posterior[k];
                        ch.partner_prob[m] = std::clamp(pp, 0.0, 1.0);
                    } else  // legacy x-ray / annihilation member
                        ch.partner_prob[m] = cascade_partner_prob(
                            dcB, static_cast<std::size_t>(mi), m);
                }

                // Vacancies coincident with the pair. The EC feed sits upstream of
                // both, so it conditions on the higher (a_hi) transition alone.
                {
                    double ecK = 0.0, ecL = 0.0;
                    const std::vector<double>& post = ch.weak_outcome_posterior;
                    if (!post.empty()) {
                        for (std::size_t k = 0; k < post.size(); ++k) {
                            const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                            if (o.kind != WeakOutcomeKind::ElectronCapture) continue;
                            ecK += post[k] * o.ec_K;
                            ecL += post[k] * o.ec_L;
                        }
                    } else {
                        const int a_hi = (dagB.ts[i].from >= dagB.ts[j].from) ? i : j;
                        const double pass_hi = dagB.pass(a_hi);
                        if (pass_hi > 0.0) {
                            for (int L = 0; L < dagB.n; ++L) {
                                const double f = dagB.feed[L] * dagB.R[a_hi][L];
                                ecK += f * dcB.level_scheme.levels[L].feed_ecK;
                                ecL += f * dcB.level_scheme.levels[L].feed_ecL;
                            }
                            ecK /= pass_hi;
                            ecL /= pass_hi;
                        }
                    }
                    if (ecK > 1e-9) ch.vacancies.push_back({
                        std::clamp(ecK, 0.0, 1.0), false, -1, 0.0, 0, -1, {}});
                    if (ecL > 1e-9) ch.vacancies.push_back({
                        std::clamp(ecL, 0.0, 1.0), true, -1, 0.0, 0, -1, {}});
                }
                // IC vacancies from third transitions coincident with the pair.
                for (int c = 0; c < nts; ++c) {
                    if (c == i || c == j) continue;
                    const double base = event_probability({i, j, c}, c) / pj;
                    if (base <= 1e-12) continue;
                    // From the transition, not its photon member -- an E0 has none.
                    const double e_tr = dagB.ts[c].gamma_keV;
                    const double icK = base * dagB.ts[c].p_icK;
                    auto outcome_shell_prob = [&](double shell_p) {
                        std::vector<double> p;
                        if (!ch.weak_outcome_posterior.empty()) {
                            p.resize(dcB.weak_outcome_law.outcomes.size(), 0.0);
                            for (std::size_t k = 0; k < p.size(); ++k) {
                                const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                                const double den = outcome_event_probability(o, {i, j});
                                p[k] = den > 0.0
                                    ? std::clamp(outcome_event_probability(
                                          o, {i, j, c}, c) / den * shell_p,
                                          0.0, 1.0)
                                    : 0.0;
                            }
                        }
                        return p;
                    };
                    if (icK > 1e-9) {
                        CascadeVacancyDraw v{std::clamp(icK, 0.0, 1.0),
                            false, -1, e_tr, c + 1, dagB.ts[c].gamma_member, {}};
                        v.weak_outcome_prob = outcome_shell_prob(dagB.ts[c].p_icK);
                        ch.vacancies.push_back(std::move(v));
                    }
                    for (int s = 0; s < 3; ++s) {
                        const double icL = base * dagB.ts[c].p_icL[s];
                        if (icL > 1e-9) {
                            CascadeVacancyDraw v{std::clamp(icL, 0.0, 1.0),
                                true, s, e_tr, c + 1, dagB.ts[c].gamma_member, {}};
                            v.weak_outcome_prob = outcome_shell_prob(dagB.ts[c].p_icL[s]);
                            ch.vacancies.push_back(std::move(v));
                        }
                    }
                    const double unresolved =
                        base * dagB.ts[c].p_icUnresolved;
                    if (unresolved > 1e-9) {
                        CascadeVacancyDraw v{std::clamp(unresolved, 0.0, 1.0),
                            false, -1, e_tr, c + 1,
                            dagB.ts[c].gamma_member, {}};
                        v.weak_outcome_prob =
                            outcome_shell_prob(dagB.ts[c].p_icUnresolved);
                        v.produces_vacancy = false;
                        ch.vacancies.push_back(std::move(v));
                    }
                }
                // Residual conversions are independent of this DAG pair.  If
                // the residual is the excluded primary, the channel is already
                // conditioned on its gamma not firing, so expose the conversion
                // split conditional on that complement and remove the gamma gate.
                for (int rix = 0;
                     rix < static_cast<int>(dcB.residual_transitions.size()); ++rix) {
                    const ResidualTransition& r = dcB.residual_transitions[
                        static_cast<std::size_t>(rix)];
                    const bool conditioned_no_gamma =
                        B == primary_branch && rix == irA;
                    const double scale = conditioned_no_gamma
                        ? 1.0 / std::max(1e-300, 1.0 - r.p_gamma) : 1.0;
                    const int gm = conditioned_no_gamma ? -1 : r.gamma_member;
                    const int group = nts + 1 + rix;
                    auto append_residual = [&](double p, bool is_L, int subshell,
                                               bool produces_vacancy) {
                        p *= scale;
                        if (p <= 1e-9) return;
                        CascadeVacancyDraw v{std::clamp(p, 0.0, 1.0), is_L,
                            subshell, r.transition_energy_keV, group, gm, {}};
                        if (!ch.weak_outcome_posterior.empty())
                            v.weak_outcome_prob.assign(
                                ch.weak_outcome_posterior.size(), v.prob);
                        v.produces_vacancy = produces_vacancy;
                        ch.vacancies.push_back(std::move(v));
                    };
                    append_residual(r.p_icK, false, -1, true);
                    append_residual(r.p_icL1, true, 0, true);
                    append_residual(r.p_icL2, true, 1, true);
                    append_residual(r.p_icL3, true, 2, true);
                    append_residual(r.p_icOuter + r.p_icUnresolved,
                                    false, -1, false);
                }
                out.push_back(std::move(ch));
            }
        }
    }

    // Pairs containing at least one categorical residual gamma.  The residual
    // occurrence is independent of the matched DAG, but remains mutually
    // exclusive with conversion of that same transition.  This is separate from
    // the DAG-DAG loop above to keep its long-established conditioning intact.
    for (std::size_t B = 0; B < cascades.size(); ++B) {
        const DecayCascade& dcB = cascades[B];
        const LevelDag dagB(dcB.level_scheme);
        if (!dagB.valid || dcB.residual_transitions.empty()) continue;
        const int nts = static_cast<int>(dagB.ts.size());
        const double bwB = dcB.branch_weight;

        auto add_pair = [&](int ri, int tj, int rj) {
            const ResidualTransition& a = dcB.residual_transitions[
                static_cast<std::size_t>(ri)];
            const bool b_is_residual = rj >= 0;
            const ResidualTransition* rb = b_is_residual
                ? &dcB.residual_transitions[static_cast<std::size_t>(rj)] : nullptr;
            const int ma = a.gamma_member;
            const int mb = b_is_residual ? rb->gamma_member : dagB.ts[tj].gamma_member;
            const double pga = a.p_gamma;
            const double pgb = b_is_residual ? rb->p_gamma : dagB.ts[tj].p_gamma;
            if (ma < 0 || mb < 0 || pga <= 0.0 || pgb <= 0.0
                || static_cast<std::size_t>(ma) >= dcB.members.size()
                || static_cast<std::size_t>(mb) >= dcB.members.size()) return;
            if (B == primary_branch &&
                (ma == static_cast<int>(primary_index)
                 || mb == static_cast<int>(primary_index))) return;
            const double Ea = dcB.members[static_cast<std::size_t>(ma)].energy_keV;
            const double Eb = dcB.members[static_cast<std::size_t>(mb)].energy_keV;
            if (std::abs(Ea + Eb - peak_keV) > tol_keV
                || std::min(Ea, Eb) < min_e) return;

            std::vector<int> subset;
            if (tj >= 0) subset.push_back(tj);
            const bool exclude_here = B == primary_branch;
            auto path_event = [&](std::vector<int> s, int converting = -1) {
                double p = s.empty()
                    ? 1.0 : weak_subset_probability(dcB, dagB, s);
                if (exclude_here && ipA >= 0 && converting != ipA) {
                    s.push_back(ipA);
                    p -= weak_subset_probability(dcB, dagB, s) * pg_p;
                }
                if (exclude_here && irA >= 0)
                    p *= 1.0 - pg_p;
                return std::max(0.0, p);
            };
            auto outcome_path_event = [&](const WeakOutcome& o,
                                          std::vector<int> s,
                                          int converting = -1) {
                double p = s.empty()
                    ? 1.0 : weak_outcome_reach(dcB, dagB, o, s);
                if (exclude_here && ipA >= 0 && converting != ipA) {
                    s.push_back(ipA);
                    p -= weak_outcome_reach(dcB, dagB, o, s) * pg_p;
                }
                if (exclude_here && irA >= 0)
                    p *= 1.0 - pg_p;
                return std::max(0.0, p);
            };
            const double pair_path = path_event(subset);
            const double pair_emit = pair_path * pga * pgb;
            if (pair_emit <= 1e-12) return;

            CascadeSumPairChannel ch;
            ch.e_a_keV = Ea;
            ch.e_b_keV = Eb;
            ch.member_a = ma;
            ch.member_b = mb;
            ch.branch = B;
            ch.weight = bwB * pair_emit / denom;
            if (ch.weight < 1e-9) return;

            if (dcB.weak_outcome_law.usable()) {
                ch.weak_outcome_posterior.assign(
                    dcB.weak_outcome_law.outcomes.size(), 0.0);
                double ps = 0.0;
                for (std::size_t k = 0;
                     k < dcB.weak_outcome_law.outcomes.size(); ++k) {
                    const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                    ch.weak_outcome_posterior[k] =
                        o.selected_mass * outcome_path_event(o, subset);
                    ps += ch.weak_outcome_posterior[k];
                }
                if (ps > 0.0)
                    for (double& p : ch.weak_outcome_posterior) p /= ps;
                else
                    ch.weak_outcome_posterior.clear();
                ch.weak_outcome_partner_prob =
                    weak_outcome_partner_probabilities(dcB, dagB, subset, ma);
            }

            ch.partner_prob.assign(dcB.members.size(), 0.0);
            for (std::size_t m = 0; m < dcB.members.size(); ++m) {
                if (static_cast<int>(m) == ma || static_cast<int>(m) == mb) continue;
                const int tm = dagB.transition_of(static_cast<int>(m));
                const int rm = residual_transition_of(dcB, static_cast<int>(m));
                if (B == primary_branch && m == primary_index) {
                    ch.partner_prob[m] = 0.0;
                    for (auto& row : ch.weak_outcome_partner_prob)
                        if (m < row.size()) row[m] = 0.0;
                } else if (tm >= 0) {
                    std::vector<int> joint = subset;
                    joint.push_back(tm);
                    ch.partner_prob[m] = std::clamp(
                        path_event(joint) * dagB.ts[tm].p_gamma / pair_path,
                        0.0, 1.0);
                    for (std::size_t k = 0;
                         k < ch.weak_outcome_partner_prob.size(); ++k) {
                        const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                        const double den_o = outcome_path_event(o, subset);
                        ch.weak_outcome_partner_prob[k][m] = den_o > 0.0
                            ? std::clamp(outcome_path_event(o, joint)
                                  * dagB.ts[tm].p_gamma / den_o, 0.0, 1.0)
                            : 0.0;
                    }
                } else if (rm >= 0) {
                    ch.partner_prob[m] = std::clamp(
                        dcB.residual_transitions[static_cast<std::size_t>(rm)].p_gamma,
                        0.0, 1.0);
                } else if (dcB.members[m].type == CascadeParticleType::Annih511
                           && !ch.weak_outcome_posterior.empty()) {
                    double pp = 0.0;
                    for (std::size_t k = 0;
                         k < ch.weak_outcome_posterior.size(); ++k)
                        if (dcB.weak_outcome_law.outcomes[k].kind ==
                            WeakOutcomeKind::Positron)
                            pp += ch.weak_outcome_posterior[k];
                    ch.partner_prob[m] = std::clamp(pp, 0.0, 1.0);
                } else {
                    ch.partner_prob[m] = cascade_partner_prob(
                        dcB, static_cast<std::size_t>(ma), m);
                }
            }

            if (ch.weak_outcome_posterior.empty()) {
                double ecK = 0.0, ecL = 0.0;
                if (tj >= 0) {
                    const double pass = dagB.pass(tj);
                    if (pass > 0.0) {
                        for (int L = 0; L < dagB.n; ++L) {
                            const double f = dagB.feed[L] * dagB.R[tj][L];
                            ecK += f * dcB.level_scheme.levels[L].feed_ecK;
                            ecL += f * dcB.level_scheme.levels[L].feed_ecL;
                        }
                        ecK /= pass;
                        ecL /= pass;
                    }
                } else {
                    const double ep = dcB.level_scheme.entry_probability;
                    for (int L = 0; L < dagB.n; ++L) {
                        ecK += ep * dagB.feed[L]
                            * dcB.level_scheme.levels[L].feed_ecK;
                        ecL += ep * dagB.feed[L]
                            * dcB.level_scheme.levels[L].feed_ecL;
                    }
                    if (dagB.n > 0) {
                        ecK += (1.0 - ep) * dcB.level_scheme.levels[0].feed_ecK;
                        ecL += (1.0 - ep) * dcB.level_scheme.levels[0].feed_ecL;
                    }
                }
                if (ecK > 1e-9) ch.vacancies.push_back({
                    std::clamp(ecK, 0.0, 1.0), false, -1, 0.0, 0, -1, {}});
                if (ecL > 1e-9) ch.vacancies.push_back({
                    std::clamp(ecL, 0.0, 1.0), true, -1, 0.0, 0, -1, {}});
            }

            // Matched-DAG conversion outcomes conditional on the pair.
            for (int c = 0; c < nts; ++c) {
                if (c == tj) continue;
                std::vector<int> joint = subset;
                joint.push_back(c);
                const double base = path_event(joint, c) / pair_path;
                if (base <= 1e-12) continue;
                const LevelDag::Tr& trc = dagB.ts[c];
                const int group = c + 1;
                auto append = [&](double p, bool is_L, int subshell,
                                  bool produces_vacancy) {
                    p *= base;
                    if (p <= 1e-9) return;
                    CascadeVacancyDraw v{std::clamp(p, 0.0, 1.0), is_L,
                        subshell, trc.gamma_keV, group, trc.gamma_member, {}};
                    if (!ch.weak_outcome_posterior.empty()) {
                        v.weak_outcome_prob.resize(
                            dcB.weak_outcome_law.outcomes.size(), 0.0);
                        for (std::size_t k = 0;
                             k < v.weak_outcome_prob.size(); ++k) {
                            const WeakOutcome& o = dcB.weak_outcome_law.outcomes[k];
                            const double den_o = outcome_path_event(o, subset);
                            v.weak_outcome_prob[k] = den_o > 0.0
                                ? std::clamp(outcome_path_event(o, joint, c)
                                      / den_o * (p / base), 0.0, 1.0)
                                : 0.0;
                        }
                    }
                    v.produces_vacancy = produces_vacancy;
                    ch.vacancies.push_back(std::move(v));
                };
                append(trc.p_icK, false, -1, true);
                append(trc.p_icL[0], true, 0, true);
                append(trc.p_icL[1], true, 1, true);
                append(trc.p_icL[2], true, 2, true);
                append(trc.p_icUnresolved, false, -1, false);
            }

            // Residual conversion outcomes, excluding the two forced gammas.
            for (int rix = 0;
                 rix < static_cast<int>(dcB.residual_transitions.size()); ++rix) {
                if (rix == ri || rix == rj) continue;
                const ResidualTransition& r = dcB.residual_transitions[
                    static_cast<std::size_t>(rix)];
                const bool conditioned_no_gamma =
                    B == primary_branch && rix == irA;
                const double scale = conditioned_no_gamma
                    ? 1.0 / std::max(1e-300, 1.0 - r.p_gamma) : 1.0;
                const int gm = conditioned_no_gamma ? -1 : r.gamma_member;
                const int group = nts + 1 + rix;
                auto append = [&](double p, bool is_L, int subshell,
                                  bool produces_vacancy) {
                    p *= scale;
                    if (p <= 1e-9) return;
                    CascadeVacancyDraw v{std::clamp(p, 0.0, 1.0), is_L,
                        subshell, r.transition_energy_keV, group, gm, {}};
                    if (!ch.weak_outcome_posterior.empty())
                        v.weak_outcome_prob.assign(
                            ch.weak_outcome_posterior.size(), v.prob);
                    v.produces_vacancy = produces_vacancy;
                    ch.vacancies.push_back(std::move(v));
                };
                append(r.p_icK, false, -1, true);
                append(r.p_icL1, true, 0, true);
                append(r.p_icL2, true, 1, true);
                append(r.p_icL3, true, 2, true);
                append(r.p_icOuter + r.p_icUnresolved, false, -1, false);
            }
            out.push_back(std::move(ch));
        };

        for (int ri = 0;
             ri < static_cast<int>(dcB.residual_transitions.size()); ++ri) {
            for (int t = 0; t < nts; ++t)
                if (dagB.ts[t].p_gamma > 0.0 && dagB.ts[t].gamma_member >= 0)
                    add_pair(ri, t, -1);
            for (int rj = ri + 1;
                 rj < static_cast<int>(dcB.residual_transitions.size()); ++rj)
                add_pair(ri, -1, rj);
        }
    }
    return out;
}

double EfficiencyCalculator::emit_ic_electron(
    double trans_keV, double binding_keV, bool enabled,
    const Eigen::Vector3d& vertex, std::vector<PathSegment>& seg_buffer,
    std::uniform_real_distribution<double>& u, std::mt19937_64& rng) const
{
    if (!enabled) return 0.0;                       // OFF: no RNG draw
    double ke = trans_keV - binding_keV;
    if (ke <= kIcElectronMin_keV) return 0.0;       // sub-threshold: local/ignore
    static const Material kAir = make_Air();
    Eigen::Vector3d pos = vertex;
    Eigen::Vector3d dir = sample_isotropic_dir(u, rng);
    double dep = 0.0;
    // (a) source self-attenuation escape walk (only with source material).
    const Material* mat = source_geometry_.source_material();
    if (source_geometry_.has_source_effects() && mat && mat->density() > 0.0) {
        const double range_cm =
            ElectronCsda::instance().range_gcm2_material(*mat, ke) / mat->density();
        const double d_stay = source_geometry_.min_distance_to_boundary(
            vertex, /*include_shields=*/true);
        if (range_cm < d_stay) return 0.0;          // contained: cannot reach crystal
        const ElectronSourceWalkResult w =
            ElectronCsda::instance().walk_in_source_geometry(
                source_geometry_, *mat, vertex, dir, ke, rng);
        for (const auto& bp : w.brems_photons)      // brems still counts
            dep += transport_cascade_member(bp.position, bp.direction,
                                            bp.energy_keV, seg_buffer, rng);
        if (!w.escaped) return dep;                 // stopped in source: brems only
        pos = w.exit_position; dir = w.exit_direction; ke = w.exit_KE_keV;
    }
    // (b) air-gap gate + degradation (air is not in the MC geometry, so the
    // electron would otherwise cross it unattenuated). The gap is the pure vacuum
    // distance from `pos` to the first detector material.
    geometry_.trace_ray(pos, dir, seg_buffer);
    const PathSegment* first = nullptr;
    for (const auto& s : seg_buffer) if (s.material) { first = &s; break; }
    if (!first) return dep;                         // misses the detector
    const double d_gap = std::max(first->t_start, 0.0);
    if (d_gap > 0.0) {
        ke = ElectronCsda::instance().residual_energy_keV(
            kAir, ke, d_gap * kAir.density());
        if (ke <= kIcElectronMin_keV) return dep;   // absorbed in the air gap
    }
    // (c) deposit in the crystal (electron path; positron term unused).
    const ElectronDepositResult er =
        ElectronCsda::instance().deposited_in_scoring(
            pos, dir, ke, geometry_, rng,
            transport_config_.disable_moliere, transport_config_.disable_brems);
    dep += er.deposited_scoring_keV;
    for (const auto& bp : er.brems_photons)
        dep += transport_into_detector(bp.position, bp.direction,
                                       bp.energy_keV, seg_buffer, rng);
    return dep;
}

double EfficiencyCalculator::emit_vacancy_xray_deposit(
    int daughter_Z, bool is_L, int l_subshell, const Eigen::Vector3d& vertex,
    std::vector<PathSegment>& seg_buffer,
    std::uniform_real_distribution<double>& u, std::mt19937_64& rng,
    bool ic_electrons) const
{
    if (daughter_Z <= 0) return 0.0;
    const FluorescenceData* fl = CrossSectionData::instance().fluorescence(daughter_Z);
    const LFluorescenceData* flL = CrossSectionData::instance().l_fluorescence(daughter_Z);

    // L-vacancy relaxation with Coster-Kronig transfer. MIRRORS emit_l_xray in
    // cascade_full_thread -- keep in sync (FR is the reference, not modified).
    auto emit_l = [&](int s) -> double {
        if (!flL) return 0.0;
        if (s < 0) { const double r = u(rng); s = (r < 0.25) ? 0 : (r < 0.5 ? 1 : 2); }
        const float* f = flL->coster_kronig;  // {f12, f13, f23}
        for (int guard = 0; guard < 4; ++guard) {
            const LSubshellFluor& sub = flL->sub[s];
            const double r = u(rng);
            if (r < sub.fluorescence_yield) {           // radiate from sub[s]
                double xi = u(rng), cum = 0.0;
                double ex = sub.num_lines ? sub.line_energy_keV[0] : 0.0;
                for (int li = 0; li < sub.num_lines; ++li) {
                    cum += sub.line_probability[li];
                    if (xi <= cum) { ex = sub.line_energy_keV[li]; break; }
                }
                if (ex <= 0.0) return 0.0;
                return transport_cascade_member(
                    vertex, sample_isotropic_dir(u, rng), ex, seg_buffer, rng);
            }
            const double w = sub.fluorescence_yield;
            if (s == 0) {                                // L1
                if (r < w + f[1]) { s = 2; continue; }   // f13: L1 -> L3
                if (r < w + f[1] + f[0]) { s = 1; continue; }  // f12: L1 -> L2
            } else if (s == 1) {                         // L2
                if (r < w + f[2]) { s = 2; continue; }   // f23: L2 -> L3
            }
            return 0.0;  // normal Auger / L3 non-radiative
        }
        return 0.0;
    };

    if (is_L) return emit_l(l_subshell);
    const double ex = sample_vacancy_xray(fl, u, rng);
    if (ex <= 0.0)  // K-Auger: deposit the K binding locally (was dropped).
        return ic_electrons && fl
            ? emit_ic_electron(fl->k_edge_keV, 0.0, ic_electrons,
                               vertex, seg_buffer, u, rng)
            : 0.0;
    double D = transport_cascade_member(
        vertex, sample_isotropic_dir(u, rng), ex, seg_buffer, rng);
    if (fl && ex < 0.93 * fl->k_edge_keV) D += emit_l(-1);  // Kalpha -> L vacancy
    return D;
}

EfficiencyCalculator::CascadePeakTally EfficiencyCalculator::cascade_peak_thread(
    const std::vector<DecayCascade>& cascades, std::size_t primary_branch,
    std::size_t primary_index, double peak_keV, double tol_keV,
    double cos_theta_max, bool coned,
    const std::vector<CascadeVacancyDraw>& prim_vacancies,
    const std::vector<CascadeSumPairChannel>& channels,
    uint64_t num_events, uint64_t seed, bool ic_electrons,
    const std::function<void(const CascadePeakTally&)>& progress_flush) const
{
    CascadePeakTally tally;
    // Progress flush cadence: clock-check every 256 events, flush the local
    // tally at >= kFlushSec per thread. Reads only local counters + the
    // clock — never the RNG — so results are bit-identical either way.
    constexpr uint64_t kCheckMask = 0xFF;
    constexpr double kFlushSec = 0.25;
    auto last_flush = std::chrono::steady_clock::now();
    const uint64_t seed_mixed = seed * 2654435761ULL;
    const std::array<uint_least32_t, 6> seed_data{{
      static_cast<uint_least32_t>(seed), static_cast<uint_least32_t>(seed >> 32),
      static_cast<uint_least32_t>(seed_mixed), static_cast<uint_least32_t>(seed_mixed >> 32),
      static_cast<uint_least32_t>(seed ^ (seed >> 16)),
      static_cast<uint_least32_t>((seed >> 48) | (seed << 16))
    }};
    std::seed_seq seq( seed_data.begin(), seed_data.end() );
    std::mt19937_64 rng(seq);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    const DecayCascade& dc = cascades[primary_branch];
    const CascadeMember& primary = dc.members[primary_index];

    // Daughter K/L binding energies for the IC-electron KE (E_IC = E_trans - E_bind).
    const FluorescenceData* fl_ic =
        ic_electrons ? CrossSectionData::instance().fluorescence(dc.daughter_Z) : nullptr;
    const LFluorescenceData* flL_ic =
        ic_electrons ? CrossSectionData::instance().l_fluorescence(dc.daughter_Z) : nullptr;
    const double k_bind = fl_ic ? static_cast<double>(fl_ic->k_edge_keV) : 0.0;
    const double l_bind = flL_ic ? static_cast<double>(flL_ic->l3_edge_keV) : 0.0;

    // Conditional emission probability of each partner, given the primary fired.
    // Prefer the level-scheme visit/reach DP (exact joints; p=0 for transitions
    // mutually exclusive with the primary); fall back to the pairwise/marginal
    // construction only when there is no usable daughter level scheme.
    const bool has_scheme = dc.level_scheme.valid;
    std::vector<double> pmate = cascade_level_pmate(dc, primary_index);
    if (pmate.empty()) {
        pmate.assign(dc.members.size(), 0.0);
        for (std::size_t m = 0; m < dc.members.size(); ++m)
            if (m != primary_index)
                pmate[m] = cascade_partner_prob(dc, primary_index, m);
    }
    const LevelDag primary_dag(dc.level_scheme);
    const int primary_transition = primary_dag.valid
        ? primary_dag.transition_of(static_cast<int>(primary_index)) : -1;
    const int primary_residual = primary_dag.valid
        ? residual_transition_of(dc, static_cast<int>(primary_index)) : -1;
    std::vector<double> primary_weak_post = primary_transition >= 0
        ? weak_outcome_posterior(dc, primary_dag, {primary_transition})
        : std::vector<double>{};
    const std::vector<std::vector<double>> primary_weak_partner =
        primary_transition >= 0
            ? weak_outcome_partner_probabilities(
                  dc, primary_dag, {primary_transition},
                  static_cast<int>(primary_index))
            : primary_residual >= 0
            ? weak_outcome_partner_probabilities(
                  dc, primary_dag, {}, static_cast<int>(primary_index))
            : std::vector<std::vector<double>>{};
    if (primary_weak_post.empty() && dc.weak_outcome_law.usable()) {
        // A rejected/metadata-free branch has no topology with which to refine
        // the mode posterior.  Its branch-level categorical law is still the
        // physically bounded fallback and must replace independent EC/beta+
        // Bernoulli draws.
        for (const WeakOutcome& o : dc.weak_outcome_law.outcomes)
            primary_weak_post.push_back(o.selected_mass);
    }

    std::vector<PathSegment> seg_buffer;

    for (uint64_t e = 0; e < num_events; ++e) {
        ++tally.n;
        // One decay = one location: sample the shared vertex per history (a
        // constant for point sources, sampled over the volume for extended).
        // The sum-fed channels reuse this vertex -- an independent draw from the
        // same source distribution, so the estimator mean is unbiased (the shared
        // vertex only correlates the streams, affecting variance, not the mean).
        const Eigen::Vector3d vertex = sample_source_position(rng);

        // Primary gamma: cone-biased toward the detector (else isotropic).
        const Eigen::Vector3d pdir = coned
            ? sample_cone_direction(vertex, cos_theta_max, rng)
            : sample_isotropic_dir(uniform, rng);
        const double dep_primary =
            transport_cascade_member(vertex, pdir, primary.energy_keV, seg_buffer, rng);

        double summed = dep_primary;
        int selected_post = -1;
        if (!primary_weak_post.empty()) {
            const double uw = uniform(rng);
            double aw = 0.0;
            for (std::size_t k = 0; k < primary_weak_post.size(); ++k) {
                aw += primary_weak_post[k];
                if (uw <= aw) { selected_post = static_cast<int>(k); break; }
            }
            if (selected_post < 0) selected_post =
                static_cast<int>(primary_weak_post.size()) - 1;
        }
        std::vector<char> primary_partner_fired(dc.members.size(), 0);
        for (std::size_t m = 0; m < dc.members.size(); ++m) {
            if (m == primary_index) continue;
            double partner_p = pmate[m];
            if (selected_post >= 0 &&
                selected_post < static_cast<int>(primary_weak_partner.size()) &&
                m < primary_weak_partner[static_cast<std::size_t>(selected_post)].size())
                partner_p = primary_weak_partner[static_cast<std::size_t>(selected_post)][m];
            if (!(partner_p > 0.0)) continue;
            bool fires = false;
            if (dc.members[m].type == CascadeParticleType::Annih511 &&
                selected_post >= 0) {
                fires = dc.weak_outcome_law.outcomes[selected_post].kind ==
                    WeakOutcomeKind::Positron;
            } else {
                fires = uniform(rng) < partner_p;
            }
            if (!fires)
                continue;
            primary_partner_fired[m] = 1;
            // Partner direction: relative to the primary via W(theta) when the
            // pair carries an angular correlation, else isotropic. The primary's
            // cone-bias weight already corrects its own pdf; the partner is drawn
            // from its physical pdf W/4pi (which integrates to 1), so no extra
            // weight is needed.
            double a2 = 0.0, a4 = 0.0;
            Eigen::Vector3d d;
            if (cascade_corr_link(dc, primary_index, m, a2, a4)) {
                const double ct = sample_cos_theta_W(a2, a4, uniform, rng);
                d = direction_at_angle(pdir, ct, kTwoPi * uniform(rng));
            } else {
                d = sample_isotropic_dir(uniform, rng);
            }
            if (dc.members[m].type == CascadeParticleType::Annih511) {
                // Back-to-back annihilation pair from the (point) vertex.
                summed += transport_cascade_member(vertex,  d, 510.998950, seg_buffer, rng);
                summed += transport_cascade_member(vertex, -d, 510.998950, seg_buffer, rng);
            } else {
                summed += transport_cascade_member(vertex, d, dc.members[m].energy_keV,
                                                   seg_buffer, rng);
            }
        }

        // Vacancy x-rays coincident with the primary (summed out of the window,
        // like FullRealization; they do NOT affect the primary-alone tally). For a
        // level-scheme nuclide these come from prim_vacancies (EC feed + IC of the
        // passed transitions); otherwise from the independent k_vacancies model,
        // conditioned on the primary being emitted.
        if (has_scheme) {
            // The EC shell belongs to the same posterior weak-mode draw as the
            // annihilation decision above.  It is therefore impossible to emit
            // both in one conditioned history.
            if (selected_post >= 0) {
                const WeakOutcome& o = dc.weak_outcome_law.outcomes[selected_post];
                if (o.kind == WeakOutcomeKind::ElectronCapture) {
                    const double us = uniform(rng);
                    if (us < o.ec_K)
                        summed += emit_vacancy_xray_deposit(dc.daughter_Z, false,
                            -1, vertex, seg_buffer, uniform, rng, ic_electrons);
                    else if (us < o.ec_K + o.ec_L)
                        summed += emit_vacancy_xray_deposit(dc.daughter_Z, true,
                            -1, vertex, seg_buffer, uniform, rng, ic_electrons);
                }
            }
            std::vector<int> sampled_groups;
            for (const CascadeVacancyDraw& v : prim_vacancies) {
                if (selected_post >= 0 && v.exclusive_group == 0) continue;
                const CascadeVacancyDraw* fired = nullptr;
                if (v.exclusive_group < 0) {
                    if (v.prob > 0.0 && uniform(rng) < v.prob) fired = &v;
                } else if (std::find(sampled_groups.begin(), sampled_groups.end(),
                                     v.exclusive_group) == sampled_groups.end()) {
                    sampled_groups.push_back(v.exclusive_group);
                    if (v.gamma_member >= 0 &&
                        v.gamma_member < static_cast<int>(primary_partner_fired.size()) &&
                        primary_partner_fired[v.gamma_member])
                        continue;
                    double gate_scale = 1.0;
                    if (v.gamma_member >= 0 &&
                        v.gamma_member < static_cast<int>(pmate.size())) {
                        double pg = pmate[v.gamma_member];
                        if (selected_post >= 0 &&
                            selected_post < static_cast<int>(primary_weak_partner.size()) &&
                            v.gamma_member < static_cast<int>(
                                primary_weak_partner[static_cast<std::size_t>(selected_post)].size()))
                            pg = primary_weak_partner[static_cast<std::size_t>(selected_post)]
                                                     [static_cast<std::size_t>(v.gamma_member)];
                        const double room = 1.0 - pg;
                        if (!(room > 0.0)) continue;
                        gate_scale = 1.0 / room;
                    }
                    const double uv = uniform(rng);
                    double av = 0.0;
                    for (const CascadeVacancyDraw& q : prim_vacancies)
                        if (q.exclusive_group == v.exclusive_group) {
                            double pq = q.prob;
                            if (selected_post >= 0 &&
                                selected_post < static_cast<int>(q.weak_outcome_prob.size()))
                                pq = q.weak_outcome_prob[static_cast<std::size_t>(selected_post)];
                            av += pq * gate_scale;
                            if (!fired && uv < av) fired = &q;
                        }
                }
                if (fired) {
                    // IC conversion electron (v.trans_keV == 0 for an EC vacancy).
                    // A shell-unresolved E0 has binding=0 and no relaxation.
                    if (fired->trans_keV > 0.0)
                        summed += emit_ic_electron(fired->trans_keV,
                                      fired->produces_vacancy
                                          ? (fired->is_L ? l_bind : k_bind) : 0.0,
                                      ic_electrons,
                                      vertex, seg_buffer, uniform, rng);
                    if (fired->produces_vacancy)
                        summed += emit_vacancy_xray_deposit(dc.daughter_Z, fired->is_L,
                                      fired->l_subshell, vertex, seg_buffer, uniform, rng,
                                      ic_electrons);
                }
            }
        } else {
            if (selected_post >= 0) {
                const WeakOutcome& o = dc.weak_outcome_law.outcomes[selected_post];
                if (o.kind == WeakOutcomeKind::ElectronCapture) {
                    const double us = uniform(rng);
                    if (us < o.ec_K)
                        summed += emit_vacancy_xray_deposit(dc.daughter_Z, false,
                            -1, vertex, seg_buffer, uniform, rng, ic_electrons);
                    else if (us < o.ec_K + o.ec_L)
                        summed += emit_vacancy_xray_deposit(dc.daughter_Z, true,
                            -1, vertex, seg_buffer, uniform, rng, ic_electrons);
                }
            }
            for (const VacancyGroup& vg : dc.vacancy_groups) {
                if (vg.kind == VacancyGroupKind::InternalConversion) {
                    if (vg.gamma_member == static_cast<int>(primary_index)) continue;
                    if (vg.gamma_member >= 0 &&
                        vg.gamma_member < static_cast<int>(primary_partner_fired.size()) &&
                        primary_partner_fired[vg.gamma_member])
                        continue;  // the complementary gamma emitted
                }
                const double uv = uniform(rng);
                double av = vg.p_none;
                if (uv < av) continue;
                bool is_L = false;
                int lshell = -1;
                av += vg.p_K;
                if (uv < av) {
                    is_L = false;
                } else {
                    av += vg.p_L1;
                    if (uv < av) { is_L = true; lshell = 0; }
                    else {
                        av += vg.p_L2;
                        if (uv < av) { is_L = true; lshell = 1; }
                        else {
                            av += vg.p_L3;
                            if (uv < av) { is_L = true; lshell = 2; }
                            else {
                                av += vg.p_outer;
                                if (uv < av)
                                    summed += emit_ic_electron(
                                        vg.transition_energy_keV, 0.0,
                                        ic_electrons, vertex, seg_buffer,
                                        uniform, rng);
                                else {
                                    av += vg.p_unmodeled;
                                    // Explicit above-pair E0 remainder: advance
                                    // probability accounting, emit nothing.
                                    if (uv < av) { /* intentionally unmodelled */ }
                                }
                                continue;
                            }
                        }
                    }
                }
                summed += emit_ic_electron(vg.transition_energy_keV,
                              is_L ? l_bind : k_bind, ic_electrons,
                              vertex, seg_buffer, uniform, rng);
                summed += emit_vacancy_xray_deposit(dc.daughter_Z, is_L,
                              lshell, vertex, seg_buffer, uniform, rng,
                              ic_electrons);
            }
            for (const KVacancySource& kv : dc.k_vacancies) {
                // An IC vacancy fires only when its gamma did NOT emit; an EC
                // vacancy (gamma_member < 0) is unconditional. Condition on the
                // primary being emitted (its own IC vacancy cannot fire).
                double p = kv.prob;
                if (kv.gamma_member >= 0) {
                    if (kv.gamma_member == static_cast<int>(primary_index)) continue;
                    p *= 1.0 - pmate[kv.gamma_member];
                }
                if (p > 0.0 && uniform(rng) < p) {
                    // IC conversion electron (only for an IC vacancy, gamma_member>=0).
                    if (kv.gamma_member >= 0 && kv.gamma_member < static_cast<int>(dc.members.size()))
                        summed += emit_ic_electron(dc.members[kv.gamma_member].energy_keV,
                                      kv.is_L ? l_bind : k_bind, ic_electrons,
                                      vertex, seg_buffer, uniform, rng);
                    summed += emit_vacancy_xray_deposit(dc.daughter_Z, kv.is_L,
                                  kv.l_subshell, vertex, seg_buffer, uniform, rng,
                                  ic_electrons);
                }
            }
        }

        // With-summing score X = I_sum + Σ_k w_k·P_k. The sum-fed pair channels
        // are separate decays (primary not emitted) whose a+b deposit can land in
        // the window; photon a is cone-biased like the primary (same cone), b is
        // W(theta)-correlated / isotropic, plus the pair's conditioned partners
        // and vacancies -- so each channel carries its own summing-out.
        double x = std::abs(summed - peak_keV) < tol_keV ? 1.0 : 0.0;
        for (const CascadeSumPairChannel& ch : channels) {
            const DecayCascade& dcb = cascades[ch.branch];
            int selected_pair_post = -1;
            if (!ch.weak_outcome_posterior.empty()) {
                const double uw = uniform(rng);
                double aw = 0.0;
                for (std::size_t k = 0; k < ch.weak_outcome_posterior.size(); ++k) {
                    aw += ch.weak_outcome_posterior[k];
                    if (uw <= aw) { selected_pair_post = static_cast<int>(k); break; }
                }
                if (selected_pair_post < 0)
                    selected_pair_post = static_cast<int>(
                        ch.weak_outcome_posterior.size()) - 1;
            }
            const Eigen::Vector3d da = coned
                ? sample_cone_direction(vertex, cos_theta_max, rng)
                : sample_isotropic_dir(uniform, rng);
            double dsum = transport_cascade_member(vertex, da, ch.e_a_keV, seg_buffer, rng);
            double a2 = 0.0, a4 = 0.0;
            Eigen::Vector3d db;
            if (ch.member_a >= 0 && ch.member_b >= 0 &&
                cascade_corr_link(dcb, ch.member_a, ch.member_b, a2, a4)) {
                const double ct = sample_cos_theta_W(a2, a4, uniform, rng);
                db = direction_at_angle(da, ct, kTwoPi * uniform(rng));
            } else {
                db = sample_isotropic_dir(uniform, rng);
            }
            dsum += transport_cascade_member(vertex, db, ch.e_b_keV, seg_buffer, rng);
            std::vector<char> pair_partner_fired(dcb.members.size(), 0);
            if (ch.member_a >= 0 && ch.member_a < static_cast<int>(pair_partner_fired.size()))
                pair_partner_fired[ch.member_a] = 1;
            if (ch.member_b >= 0 && ch.member_b < static_cast<int>(pair_partner_fired.size()))
                pair_partner_fired[ch.member_b] = 1;
            for (std::size_t m = 0; m < ch.partner_prob.size(); ++m) {
                double partner_p = ch.partner_prob[m];
                if (selected_pair_post >= 0 &&
                    selected_pair_post < static_cast<int>(ch.weak_outcome_partner_prob.size()) &&
                    m < ch.weak_outcome_partner_prob[
                        static_cast<std::size_t>(selected_pair_post)].size())
                    partner_p = ch.weak_outcome_partner_prob[
                        static_cast<std::size_t>(selected_pair_post)][m];
                if (!(partner_p > 0.0)) continue;
                bool fires = false;
                if (dcb.members[m].type == CascadeParticleType::Annih511 &&
                    selected_pair_post >= 0)
                    fires = dcb.weak_outcome_law.outcomes[selected_pair_post].kind ==
                        WeakOutcomeKind::Positron;
                else
                    fires = uniform(rng) < partner_p;
                if (!fires) continue;
                pair_partner_fired[m] = 1;
                if (dcb.members[m].type == CascadeParticleType::Annih511) {
                    const Eigen::Vector3d d = sample_isotropic_dir(uniform, rng);
                    dsum += transport_cascade_member(vertex,  d, 510.998950, seg_buffer, rng);
                    dsum += transport_cascade_member(vertex, -d, 510.998950, seg_buffer, rng);
                } else {
                    dsum += transport_cascade_member(vertex, sample_isotropic_dir(uniform, rng),
                                                     dcb.members[m].energy_keV, seg_buffer, rng);
                }
            }
            const FluorescenceData* flb =
                ic_electrons ? CrossSectionData::instance().fluorescence(dcb.daughter_Z) : nullptr;
            const LFluorescenceData* flLb =
                ic_electrons ? CrossSectionData::instance().l_fluorescence(dcb.daughter_Z) : nullptr;
            const double kb = flb ? static_cast<double>(flb->k_edge_keV) : 0.0;
            const double lb = flLb ? static_cast<double>(flLb->l3_edge_keV) : 0.0;
            if (selected_pair_post >= 0) {
                const WeakOutcome& o =
                    dcb.weak_outcome_law.outcomes[selected_pair_post];
                if (o.kind == WeakOutcomeKind::ElectronCapture) {
                    const double us = uniform(rng);
                    if (us < o.ec_K)
                        dsum += emit_vacancy_xray_deposit(dcb.daughter_Z, false,
                            -1, vertex, seg_buffer, uniform, rng, ic_electrons);
                    else if (us < o.ec_K + o.ec_L)
                        dsum += emit_vacancy_xray_deposit(dcb.daughter_Z, true,
                            -1, vertex, seg_buffer, uniform, rng, ic_electrons);
                }
            }
            std::vector<int> sampled_groups;
            for (const CascadeVacancyDraw& v : ch.vacancies) {
                if (selected_pair_post >= 0 && v.exclusive_group == 0) continue;
                const CascadeVacancyDraw* fired = nullptr;
                if (v.exclusive_group < 0) {
                    if (v.prob > 0.0 && uniform(rng) < v.prob) fired = &v;
                } else if (std::find(sampled_groups.begin(), sampled_groups.end(),
                                     v.exclusive_group) == sampled_groups.end()) {
                    sampled_groups.push_back(v.exclusive_group);
                    if (v.gamma_member >= 0 &&
                        v.gamma_member < static_cast<int>(pair_partner_fired.size()) &&
                        pair_partner_fired[v.gamma_member])
                        continue;
                    double gate_scale = 1.0;
                    if (v.gamma_member >= 0 &&
                        v.gamma_member < static_cast<int>(ch.partner_prob.size())) {
                        double pg = ch.partner_prob[v.gamma_member];
                        if (selected_pair_post >= 0 &&
                            selected_pair_post < static_cast<int>(
                                ch.weak_outcome_partner_prob.size()) &&
                            v.gamma_member < static_cast<int>(
                                ch.weak_outcome_partner_prob[
                                    static_cast<std::size_t>(selected_pair_post)].size()))
                            pg = ch.weak_outcome_partner_prob[
                                static_cast<std::size_t>(selected_pair_post)]
                                [static_cast<std::size_t>(v.gamma_member)];
                        const double room = 1.0 - pg;
                        if (!(room > 0.0)) continue;
                        gate_scale = 1.0 / room;
                    }
                    const double uv = uniform(rng);
                    double av = 0.0;
                    for (const CascadeVacancyDraw& q : ch.vacancies)
                        if (q.exclusive_group == v.exclusive_group) {
                            double pq = q.prob;
                            if (selected_pair_post >= 0 &&
                                selected_pair_post < static_cast<int>(
                                    q.weak_outcome_prob.size()))
                                pq = q.weak_outcome_prob[
                                    static_cast<std::size_t>(selected_pair_post)];
                            av += pq * gate_scale;
                            if (!fired && uv < av) fired = &q;
                        }
                }
                if (fired) {
                    if (fired->trans_keV > 0.0)
                        dsum += emit_ic_electron(fired->trans_keV,
                                    fired->produces_vacancy
                                        ? (fired->is_L ? lb : kb) : 0.0,
                                    ic_electrons, vertex, seg_buffer, uniform, rng);
                    if (fired->produces_vacancy)
                        dsum += emit_vacancy_xray_deposit(dcb.daughter_Z, fired->is_L,
                                    fired->l_subshell, vertex, seg_buffer, uniform, rng,
                                    ic_electrons);
                }
            }
            if (std::abs(dsum - peak_keV) < tol_keV) x += ch.weight;
        }

        const double nosum = std::abs(dep_primary - peak_keV) < tol_keV ? 1.0 : 0.0;
        if (nosum > 0.0) ++tally.n_nosum;
        tally.sum_x  += x;
        tally.sum_xx += x * x;
        tally.sum_xd += x * nosum;

        // RNG-neutral progress flush (reads local counters + clock only).
        if (progress_flush && (e & kCheckMask) == kCheckMask) {
            const auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration<double>(now - last_flush).count() >= kFlushSec) {
                last_flush = now;
                progress_flush(tally);
            }
        }
    }
    return tally;
}

EfficiencyCalculator::CascadeFullTally EfficiencyCalculator::cascade_full_thread(
    const std::vector<DecayCascade>& cascades,
    const std::vector<double>& branch_cdf, double branch_weight_sum,
    const std::vector<CascadePeakTarget>& targets,
    const std::vector<float>& spectrum_edges,
    uint64_t num_events, uint64_t seed, bool ic_electrons,
    const std::function<void(const CascadeFullTally&)>& progress_flush) const
{
    CascadeFullTally tally;
    tally.branch_weight_sum = branch_weight_sum;
    tally.n_emit.assign(targets.size(), 0);
    tally.n_nosum.assign(targets.size(), 0);
    tally.n_sum.assign(targets.size(), 0);
    if (spectrum_edges.size() >= 2)
        tally.spectrum.assign(spectrum_edges.size() - 1, 0.0);

    // Progress flush cadence (RNG-neutral); see cascade_peak_thread.
    constexpr uint64_t kCheckMask = 0xFF;
    constexpr double kFlushSec = 0.25;
    auto last_flush = std::chrono::steady_clock::now();

    const uint64_t seed_mixed = seed * 2654435761ULL;
    const std::array<uint_least32_t, 6> seed_data{{
      static_cast<uint_least32_t>(seed), static_cast<uint_least32_t>(seed >> 32),
      static_cast<uint_least32_t>(seed_mixed), static_cast<uint_least32_t>(seed_mixed >> 32),
      static_cast<uint_least32_t>(seed ^ (seed >> 16)),
      static_cast<uint_least32_t>((seed >> 48) | (seed << 16))
    }};
    std::seed_seq seq( seed_data.begin(), seed_data.end() );
    std::mt19937_64 rng(seq);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    // Per-branch coherent fallback forest (constant; branches are few).
    std::vector<std::vector<CascadeFallbackNode>> fallback_forests(cascades.size());
    std::vector<LevelDag> level_dags;
    level_dags.reserve(cascades.size());
    for (std::size_t b = 0; b < cascades.size(); ++b) {
        level_dags.emplace_back(cascades[b].level_scheme);
        fallback_forests[b] = build_cascade_fallback_forest(cascades[b]);
    }

    std::vector<char> emitted;
    std::vector<double> dep;
    std::vector<Eigen::Vector3d> mdir;  // per-member emission directions (this event)
    std::vector<PathSegment> seg_buffer;

    for (uint64_t e = 0; e < num_events; ++e) {
        ++tally.n;
        // One decay = one location: a shared vertex per history.
        const Eigen::Vector3d vertex = sample_source_position(rng);

        // Sample one decay branch ~ branch_weight, then its emitted member set.
        const double r = uniform(rng);
        std::size_t b = static_cast<std::size_t>(
            std::upper_bound(branch_cdf.begin(), branch_cdf.end(), r) - branch_cdf.begin());
        if (b >= cascades.size()) b = cascades.size() - 1;
        const DecayCascade& dc = cascades[b];
        const std::size_t nm = dc.members.size();
        dep.assign(nm, 0.0);
        mdir.assign(nm, Eigen::Vector3d::Zero());
        emitted.assign(nm, 0);
        double D = 0.0;

        // Select the weak decay mode exactly once.  The adapter's law is
        // branch-conditional and categorical: an EC vacancy and a positron can
        // therefore never be produced by the same branch firing.  Keep the
        // individual positron endpoint instead of the old intensity-weighted
        // aggregate so source ranging follows the selected beta spectrum.
        const WeakOutcome* selected_weak = nullptr;
        if (dc.weak_outcome_law.usable()) {
            double rw = uniform(rng), aw = 0.0;
            for (const WeakOutcome& o : dc.weak_outcome_law.outcomes) {
                aw += o.selected_mass;
                if (rw <= aw) { selected_weak = &o; break; }
            }
            if (!selected_weak && !dc.weak_outcome_law.outcomes.empty())
                selected_weak = &dc.weak_outcome_law.outcomes.back();
        }
        bool fire_positron = selected_weak &&
            selected_weak->kind == WeakOutcomeKind::Positron;
        double selected_endpoint = fire_positron
            ? selected_weak->positron_endpoint_keV : 0.0;
        int selected_ec_shell = 0;  // 0 none/outer, 1 K, 2 unresolved L
        if (selected_weak && selected_weak->kind == WeakOutcomeKind::ElectronCapture) {
            const double us = uniform(rng);
            if (us < selected_weak->ec_K) selected_ec_shell = 1;
            else if (us < selected_weak->ec_K + selected_weak->ec_L)
                selected_ec_shell = 2;
        }

        // Emit the L x-ray for a vacancy in subshell `s` (0/1/2 = L1/L2/L3); a
        // negative `s` (electron-capture / Kalpha-induced vacancy, no subshell
        // resolution) samples the subshell by L-shell occupancy (L1:L2:L3 = 2:2:4).
        // At each subshell the vacancy RADIATES (prob omega_s), transfers down a
        // subshell by a radiationless Coster-Kronig transition (L1 -> L3/L2, L2 ->
        // L3), or normal-Augers (no L x-ray) -- these compete, omega + CK + Auger
        // = 1 (the Krause-yield normalization, matching G4's EADL deexcitation).
        // Radiating BEFORE testing CK is essential: e.g. ~16% of Np L1 vacancies
        // emit an L1 x-ray (Lbeta3,4/Lgamma) and must not be pre-empted by the
        // (large) L1->L3 transfer. One vacancy in => one radiating vacancy out.
        // Returns the energy deposited in the crystal (0 if Auger / no data). The
        // L line table is pre-filtered to >= 10 keV.
        auto emit_l_xray = [&](const LFluorescenceData* flL, int s) -> double {
            if (!flL) return 0.0;
            if (s < 0) { const double r = uniform(rng);
                         s = (r < 0.25) ? 0 : (r < 0.5 ? 1 : 2); }
            const float* f = flL->coster_kronig;  // {f12, f13, f23}
            for (int guard = 0; guard < 4; ++guard) {
                const LSubshellFluor& sub = flL->sub[s];
                const double r = uniform(rng);
                if (r < sub.fluorescence_yield) {  // radiate from sub[s]
                    double xi = uniform(rng), cum = 0.0;
                    double ex = sub.num_lines ? sub.line_energy_keV[0] : 0.0;
                    for (int li = 0; li < sub.num_lines; ++li) {
                        cum += sub.line_probability[li];
                        if (xi <= cum) { ex = sub.line_energy_keV[li]; break; }
                    }
                    if (ex <= 0.0) return 0.0;
                    return transport_cascade_member(
                        vertex, sample_isotropic_dir(uniform, rng), ex, seg_buffer, rng);
                }
                // Coster-Kronig transfer (radiationless), else normal Auger (done).
                const double w = sub.fluorescence_yield;
                if (s == 0) {                                  // L1
                    if (r < w + f[1]) { s = 2; continue; }     // f13: L1 -> L3
                    if (r < w + f[1] + f[0]) { s = 1; continue; } // f12: L1 -> L2
                }
                else if (s == 1) {                             // L2
                    if (r < w + f[2]) { s = 2; continue; }     // f23: L2 -> L3
                }
                return 0.0;  // normal Auger / L3 non-radiative: no L x-ray
            }
            return 0.0;
        };

        // Emit the annihilation photons for a fired Annih511 member, ranging the
        // positron (KE sampled from the beta+ spectrum up to `endpoint`) to its
        // annihilation site instead of point-annihilating at the vertex:
        //   (1) in-flight annihilation -> two non-511 photons (off the 511 peak);
        //   (2) a positron that ESCAPES the source geometry leaves no clean source
        //       511 pair (it annihilates downstream / in the crystal as a
        //       high-energy event -- not re-emitted here, consistent with the
        //       photon-only cascade), so its source 511 pair is suppressed.
        // A contained positron annihilates at rest near the vertex (the sub-mm
        // endpoint offset is negligible at the source-detector distance).
        // endpoint <= 0 or no source material => legacy back-to-back 511 at vertex.
        // Returns the energy deposited in the crystal (keV).
        auto emit_annihilation = [&](double endpoint) -> double {
            constexpr double mc2 = 510.998950;
            auto at_rest = [&]() {
                const Eigen::Vector3d d = sample_isotropic_dir(uniform, rng);
                return transport_cascade_member(vertex,  d, mc2, seg_buffer, rng)
                     + transport_cascade_member(vertex, -d, mc2, seg_buffer, rng);
            };
            const Material* mat = source_geometry_.source_material();
            if (endpoint <= 0.0 || mat == nullptr ||
                !source_geometry_.has_source_effects())
                return at_rest();  // legacy / no material to range in
            const double T = sample_beta_plus_KE(endpoint, uniform, rng);
            // (1) in-flight annihilation: two photons summing to 2*mc2 + T, split
            // by a sampled fraction, emitted isotropically -> continuum, not 511.
            if (uniform(rng) < inflight_annih_prob(*mat, T)) {
                const double Etot = 2.0 * mc2 + T;
                const double x = 0.5 + (uniform(rng) - 0.5) * (T / Etot);
                const double E1 = x * Etot;
                return transport_cascade_member(
                           vertex, sample_isotropic_dir(uniform, rng), E1, seg_buffer, rng)
                     + transport_cascade_member(
                           vertex, sample_isotropic_dir(uniform, rng), Etot - E1, seg_buffer, rng);
            }
            // (2) range the positron; escape => no clean source 511 pair.
            // Containment fast-path (mirrors SourceGeometry::process_electron):
            // a positron whose CSDA range in `mat` is below the distance to the
            // nearest source boundary cannot escape -> annihilate at rest and
            // skip the walk entirely (the common case for a soft beta+ spectrum
            // in a finite source; a net speedup).  The skin-escape gate is now
            // regime-aware at the exit boundary, so a soft positron that DOES
            // reach the surface is kept as an escape automatically -- no
            // per-caller skip flag is needed (was apply_skin_escape_gate=false,
            // which the gate's exit-energy window now subsumes).
            if (mat->density() > 0.0) {
                const double range_cm =
                    ElectronCsda::instance().range_gcm2_material(*mat, T)
                    / mat->density();
                const double d_stay = source_geometry_.min_distance_to_boundary(
                    vertex, /*include_shields=*/true);
                if (range_cm < d_stay) return at_rest();  // contained: no walk
            }
            const ElectronSourceWalkResult w =
                ElectronCsda::instance().walk_in_source_geometry(
                    source_geometry_, *mat, vertex,
                    sample_isotropic_dir(uniform, rng), T, rng);
            if (w.escaped) return 0.0;
            return at_rest();
        };

        // Shorthand for the IC/Auger electron deposit (member emit_ic_electron;
        // see its declaration). RNG-neutral when the flag is off, so the legacy
        // photon-only path stays bit-identical.
        auto ic_e = [&](double trans_keV, double binding_keV) -> double {
            return emit_ic_electron(trans_keV, binding_keV, ic_electrons,
                                    vertex, seg_buffer, uniform, rng);
        };

      if (dc.level_scheme.valid) {
        // ===== level-path realization: walk the daughter level scheme =====
        const LevelScheme& ls = dc.level_scheme;
        int n_kvac = 0;
        int n_lvac[3] = {0, 0, 0};  // resolved IC vacancies per subshell L1/L2/L3
        int n_lvac_unres = 0;       // EC / Kalpha-induced (subshell by occupancy)
        // K/L binding energies of the daughter atom, for the IC-electron KE
        // (E_IC = E_transition - E_binding(shell)). Only l3_edge is tabulated;
        // it stands in for L1/L2 (a few keV high, negligible vs the peak window).
        const FluorescenceData* fl_ic =
            ic_electrons ? CrossSectionData::instance().fluorescence(ls.daughter_Z) : nullptr;
        const LFluorescenceData* flL_ic =
            ic_electrons ? CrossSectionData::instance().l_fluorescence(ls.daughter_Z) : nullptr;
        const double k_bind = fl_ic ? static_cast<double>(fl_ic->k_edge_keV) : 0.0;
        const double l_bind = flL_ic ? static_cast<double>(flL_ic->l3_edge_keV) : 0.0;
        // Pick the directly-fed level.  Matched EC/beta+ records identify their
        // level exactly (including level zero and nonzero terminal isomers);
        // metadata-free outcomes retain the audited legacy feeding draw.
        double tot_feed = 0.0;
        for (const auto& lv : ls.levels) tot_feed += lv.feeding;
        int L = 0;
        // Not every decay of the branch enters the level scheme: some feed the
        // daughter's ground (or isomeric) state directly and emit no cascade
        // gamma. Walking unconditionally over-produces the scheme's emissions by
        // 1/entry_probability (3.5x for Th-234 -> Pa-234m).
        //
        // Drawing over max(1, tot_feed) folds that into the level draw itself
        // rather than adding a separate Bernoulli: when tot_feed < 1 the draw
        // falls past the last level with probability 1 - tot_feed and L stays 0
        // (the ground state; levels[0].feeding is identically 0 because the
        // ground state has no out-transitions), and conditional on landing, the
        // level is still picked in proportion to feeding. When tot_feed >= 1 this
        // is exactly the original draw. Same number of RNG draws either way, so
        // every branch that always enters stays bit-identical -- and does so
        // because the distribution is unchanged, not because of a floating-point
        // comparison against 1.0.
        if (selected_weak && selected_weak->fed_level >= 0 &&
            selected_weak->fed_level < static_cast<int>(ls.levels.size())) {
            L = selected_weak->fed_level;
        } else if (selected_weak && dc.weak_outcome_law.usable()) {
            // Unknown-level outcomes may use only the graph feeding left after
            // subtracting every exact weak outcome. Reusing the full aggregate
            // feeding here would double-count known EC/beta+ destinations and
            // make an unknown record borrow their cascades. The residual weights
            // are conditional on this unknown outcome and may sum below one;
            // the remainder is a terminal/direct landing (L stays zero).
            const std::vector<double> residual =
                weak_unknown_start_weights(dc, level_dags[b]);
            const double rf = uniform(rng);
            double acc = 0.0;
            for (std::size_t li = 0; li < residual.size(); ++li) {
                acc += residual[li];
                if (rf <= acc) { L = static_cast<int>(li); break; }
            }
        } else {
            double rf = uniform(rng) * std::max(1.0, tot_feed), acc = 0.0;
            for (std::size_t li = 0; li < ls.levels.size(); ++li) {
                acc += ls.levels[li].feeding;
                if (rf <= acc) { L = static_cast<int>(li); break; }
            }
        }
        // Decay-mode outcome for this realization, drawn ONCE.
        //
        // Electron capture and beta+ are competing modes of the same nuclear
        // transition, and a capture removes exactly one electron from one shell.
        // So the positron, a K capture and an L capture are mutually exclusive,
        // and sampling them with independent draws lets a single decay be both --
        // for Rb-82 (95% beta+, 5% EC) that put a spurious ~0.039 K-vacancy per
        // decay in coincidence with the 511 pair. Partitioning one uniform keeps
        // every marginal exactly as before while forbidding the overlap.
        //
        // The capture draw deliberately sits OUTSIDE the level-entry gate: a
        // capture makes its vacancy whatever level it feeds, and L is 0 for a
        // decay that did not enter the scheme -- exactly the level a
        // ground-state-fed capture populates.
        if (selected_weak) {
            if (selected_ec_shell == 1) ++n_kvac;
            else if (selected_ec_shell == 2) ++n_lvac_unres;
        } else {
            double p_pos = 0.0;
            for (std::size_t m = 0; m < nm; ++m)
                if (dc.members[m].type == CascadeParticleType::Annih511)
                    p_pos += dc.members[m].intensity;
            const bool lvl_ok = (L >= 0 && L < static_cast<int>(ls.levels.size()));
            const double ecK = lvl_ok ? ls.levels[L].feed_ecK : 0.0;
            const double ecL = lvl_ok ? ls.levels[L].feed_ecL : 0.0;
            // Partitioning one uniform requires the three probabilities to fit in
            // it. They need not: p_pos is per decay while feed_ec* is conditional
            // on landing at L, and the data gives no way to make them commensurate
            // (this compatibility path is reached only for host/legacy data with
            // no WeakOutcomeLaw). Where they do fit, partition and the exclusivity
            // is exact with every marginal preserved; where they do not, keep the
            // independent draws rather than silently truncate the last category.
            if (p_pos + ecK + ecL <= 1.0) {
                if (p_pos > 0.0 || ecK > 0.0 || ecL > 0.0) {
                    const double u = uniform(rng);
                    if (u < p_pos) fire_positron = true;
                    else if (u < p_pos + ecK) ++n_kvac;
                    else if (u < p_pos + ecK + ecL) ++n_lvac_unres;
                }
            } else {
                if (p_pos > 0.0 && uniform(rng) < p_pos) fire_positron = true;
                if (ecK > 0.0 && uniform(rng) < ecK) ++n_kvac;
                if (ecL > 0.0 && uniform(rng) < ecL) ++n_lvac_unres;
            }
        }
        // Walk down to the ground state, one transition per level.
        int prev_member = -1, guard = 0;
        while (L > 0 && L < static_cast<int>(ls.levels.size()) && guard++ < 2000) {
            const CascadeLevel& lvl = ls.levels[L];
            if (lvl.out.empty()) break;
            double tw = 0.0; for (const auto& t : lvl.out) tw += t.weight;
            if (tw <= 0.0) break;
            double rt = uniform(rng) * tw, acct = 0.0;
            const LevelOutTransition* sel = &lvl.out.front();
            for (const auto& t : lvl.out) { acct += t.weight; if (rt <= acct) { sel = &t; break; } }
            const double u = uniform(rng);
            // gamma_member is -1 for an E0 transition (no photon member). p_gamma
            // is 0 there so this branch is unreachable, but cascade_corr_link below
            // takes std::size_t and -1 would index at SIZE_MAX, so make the
            // precondition explicit rather than resting on an exact 0.0 compare.
            if (u < sel->p_gamma && sel->gamma_member >= 0
                && sel->gamma_member < static_cast<int>(nm)) {
                // Emit the gamma (W(theta)-correlated to the previous gamma).
                const int m = sel->gamma_member;
                Eigen::Vector3d d;
                double a2 = 0.0, a4 = 0.0;
                if (prev_member >= 0 && cascade_corr_link(dc, prev_member, m, a2, a4))
                    d = direction_at_angle(mdir[prev_member],
                            sample_cos_theta_W(a2, a4, uniform, rng), kTwoPi * uniform(rng));
                else
                    d = sample_isotropic_dir(uniform, rng);
                if (m >= 0 && m < static_cast<int>(nm)) {
                    mdir[m] = d; emitted[m] = 1;
                    const double dm = transport_cascade_member(
                        vertex, d, dc.members[m].energy_keV, seg_buffer, rng);
                    dep[m] = dm; D += dm; prev_member = m;
                }
            } else if (u < sel->p_gamma + sel->p_icK) {
                ++n_kvac;  // internal conversion in the K shell -> K vacancy
                D += ic_e(sel->gamma_keV, k_bind);  // conversion e-
            } else if (u < sel->p_gamma + sel->p_icK + sel->p_icL1) {
                ++n_lvac[0];  // IC in L1 -> L1 vacancy
                D += ic_e(sel->gamma_keV, l_bind);
            } else if (u < sel->p_gamma + sel->p_icK + sel->p_icL1 + sel->p_icL2) {
                ++n_lvac[1];  // IC in L2 -> L2 vacancy
                D += ic_e(sel->gamma_keV, l_bind);
            } else if (u < sel->p_gamma + sel->p_icK + sel->p_icL1 + sel->p_icL2
                                                     + sel->p_icL3) {
                ++n_lvac[2];  // IC in L3 -> L3 vacancy
                D += ic_e(sel->gamma_keV, l_bind);
            } else if (u < sel->p_gamma + sel->p_icK + sel->p_icL1
                              + sel->p_icL2 + sel->p_icL3
                              + sel->p_icUnresolved) {
                // Below-pair memberless E0 with no resolved shell. Approximate
                // the neutral-atom energy release by one full-transition-energy
                // electron, but do not invent an atomic vacancy/x-ray.
                D += ic_e(sel->gamma_keV, 0.0);
            } else if (u < sel->p_gamma + sel->p_icK + sel->p_icL1
                              + sel->p_icL2 + sel->p_icL3
                              + sel->p_icUnresolved + sel->p_unmodeled) {
                // Above-pair E0 remainder: advance the nuclear level path, but
                // emit/deposit nothing until internal-pair formation is modelled.
            }  // else: ordinary outer-shell (M+) conversion electron not tracked -- a
               // near-full-energy electron, but depositing it (validated vs
               // GEANT4) mis-shifts the low-E converting transition's OWN peak
               // and over-corrects at 2 cm; the K/L electrons already give <1%
               // agreement at 2 cm, so M+ is left out (residual only at bare-
               // source contact/2pi, where real source material absorbs it).
            L = sel->to_level;
        }

        // Transitions which could not be placed on the accepted level graph are
        // still real branch products.  Each has one absolute categorical draw;
        // this preserves gamma/conversion exclusivity and prevents the old valid-
        // scheme path from silently dropping unmatched photons.  With no topology
        // available, distinct residual transitions are independent of each other
        // and of the matched walk.
        for (const ResidualTransition& r : dc.residual_transitions) {
            const double u = uniform(rng);
            double a = r.p_none;
            if (u < a) continue;
            a += r.p_gamma;
            if (u < a) {
                const int m = r.gamma_member;
                if (m < 0 || m >= static_cast<int>(nm) || emitted[m])
                    continue;
                Eigen::Vector3d d;
                double a2 = 0.0, a4 = 0.0;
                std::size_t ref = nm;
                for (std::size_t k = 0; k < nm; ++k)
                    if (emitted[k] && cascade_corr_link(dc, k,
                            static_cast<std::size_t>(m), a2, a4)) {
                        ref = k;
                        break;
                    }
                if (ref < nm)
                    d = direction_at_angle(mdir[ref],
                        sample_cos_theta_W(a2, a4, uniform, rng),
                        kTwoPi * uniform(rng));
                else
                    d = sample_isotropic_dir(uniform, rng);
                emitted[m] = 1;
                mdir[m] = d;
                const double dm = transport_cascade_member(
                    vertex, d, dc.members[m].energy_keV, seg_buffer, rng);
                dep[m] = dm;
                D += dm;
                continue;
            }
            a += r.p_icK;
            if (u < a) { ++n_kvac; D += ic_e(r.transition_energy_keV, k_bind); continue; }
            a += r.p_icL1;
            if (u < a) { ++n_lvac[0]; D += ic_e(r.transition_energy_keV, l_bind); continue; }
            a += r.p_icL2;
            if (u < a) { ++n_lvac[1]; D += ic_e(r.transition_energy_keV, l_bind); continue; }
            a += r.p_icL3;
            if (u < a) { ++n_lvac[2]; D += ic_e(r.transition_energy_keV, l_bind); continue; }
            a += r.p_icOuter;
            if (u < a) { D += ic_e(r.transition_energy_keV, 0.0); continue; }
            a += r.p_icUnresolved;
            if (u < a) { D += ic_e(r.transition_energy_keV, 0.0); continue; }
            a += r.p_unmodeled;
            if (u < a) continue;
        }
        // Relax the K and L vacancies into x-rays (one line each; Auger
        // otherwise). A K x-ray from the L shell (Kalpha) leaves an L vacancy.
        if ((n_kvac > 0 || n_lvac[0] + n_lvac[1] + n_lvac[2] + n_lvac_unres > 0) &&
            ls.daughter_Z > 0) {
            const FluorescenceData* fl =
                CrossSectionData::instance().fluorescence(ls.daughter_Z);
            const LFluorescenceData* flL =
                CrossSectionData::instance().l_fluorescence(ls.daughter_Z);
            for (int v = 0; v < n_kvac; ++v) {
                const double ex = sample_vacancy_xray(fl, uniform, rng);
                if (ex <= 0.0) {
                    // Auger: the K binding energy is carried off by Auger
                    // electrons (previously dropped). Deposit it as a local-ish
                    // electron so E_IC + E_binding == E_transition when contained.
                    D += ic_e(k_bind, 0.0);
                    continue;
                }
                D += transport_cascade_member(
                    vertex, sample_isotropic_dir(uniform, rng), ex, seg_buffer, rng);
                if (fl && ex < 0.93 * fl->k_edge_keV) ++n_lvac_unres;  // Kalpha -> L vacancy
            }
            for (int s = 0; s < 3; ++s)
                for (int v = 0; v < n_lvac[s]; ++v) D += emit_l_xray(flL, s);
            for (int v = 0; v < n_lvac_unres; ++v) D += emit_l_xray(flL, -1);
        }
        // Kept non-K x-ray members (e.g. L x-rays the K-only vacancy model does
        // not reproduce): emitted independently (transition-group) at their
        // intensity, coincident with this decay's gammas.
        for (std::size_t m = 0; m < nm; ++m) {
            if (dc.members[m].type != CascadeParticleType::Xray) continue;
            if (uniform(rng) >= dc.members[m].intensity) continue;
            emitted[m] = 1;
            const Eigen::Vector3d d = sample_isotropic_dir(uniform, rng);
            const double dm = transport_cascade_member(
                vertex, d, dc.members[m].energy_keV, seg_buffer, rng);
            dep[m] = dm; D += dm;
        }
        // Annihilation 511 members (beta+ branch): fire by intensity (the
        // positron fraction), coincident with the level-path gammas (they sum
        // into D). The level-path walk only handles the daughter de-excitation
        // gammas, so without this the 511s would be dropped for any beta+ emitter
        // whose daughter builds a valid level scheme (Na-22, F-18, ...).
        // Fired (or not) by the single decay-mode draw above, so that a positron
        // and an electron-capture vacancy can never both occur in one decay.
        for (std::size_t m = 0; m < nm; ++m) {
            if (dc.members[m].type != CascadeParticleType::Annih511) continue;
            if (!fire_positron) continue;
            emitted[m] = 1;
            const double dm = emit_annihilation(
                selected_endpoint > 0.0 ? selected_endpoint
                                        : dc.members[m].positron_endpoint_keV);
            dep[m] = dm; D += dm;
        }
      } else {
        cascade_sample_realization(
            uniform, rng, fallback_forests[b],
            dc.weak_outcome_law.usable() ? (fire_positron ? 1 : 0) : -1,
            emitted);
        // Transport in the same parent-before-child forest order used for the
        // occurrence law. This guarantees a selected parent direction exists
        // even when its member index is numerically larger than the child's,
        // and prevents a discarded legacy link from supplying W(theta).
        for (const CascadeFallbackNode& node : fallback_forests[b]) {
            const std::size_t m = node.member;
            if (m >= nm) continue;
            if (!emitted[m]) continue;
            Eigen::Vector3d d;
            double a2 = 0.0, a4 = 0.0;
            std::size_t ref = nm;
            if (node.parent_member >= 0) {
                const std::size_t p = static_cast<std::size_t>(node.parent_member);
                // emit_annihilation samples its own back-to-back direction, so
                // mdir[p] is not a usable angular reference for a later child.
                if (p < nm && emitted[p]
                    && dc.members[p].type != CascadeParticleType::Annih511
                    && cascade_corr_link(dc, p, m, a2, a4))
                    ref = p;
            }
            if (ref < nm) {
                const double ct = sample_cos_theta_W(a2, a4, uniform, rng);
                d = direction_at_angle(mdir[ref], ct, kTwoPi * uniform(rng));
            } else {
                d = sample_isotropic_dir(uniform, rng);
            }
            mdir[m] = d;
            double dm = 0.0;
            if (dc.members[m].type == CascadeParticleType::Annih511) {
                // Positron-ranged annihilation (in-flight + escape); isotropic
                // pair, so the W(theta)-sampled `d` above is not used here.
                dm = emit_annihilation(selected_endpoint > 0.0
                    ? selected_endpoint : dc.members[m].positron_endpoint_keV);
            } else {
                dm = transport_cascade_member(vertex, d, dc.members[m].energy_keV, seg_buffer, rng);
            }
            dep[m] = dm;
            D += dm;
        }

        // Vacancy-level x-ray model: emit one x-ray (or Auger) per fired K/L
        // vacancy from the daughter element's fluorescence data, summed into D
        // alongside the gammas. One line per vacancy enforces line exclusivity;
        // multiple vacancies (EC + IC) give x-ray-x-ray summing.
        if ((!dc.k_vacancies.empty() || !dc.vacancy_groups.empty() ||
             selected_ec_shell != 0) && dc.daughter_Z > 0) {
            const FluorescenceData* fl =
                CrossSectionData::instance().fluorescence(dc.daughter_Z);
            const LFluorescenceData* flL =
                CrossSectionData::instance().l_fluorescence(dc.daughter_Z);
            const double k_bind = fl ? static_cast<double>(fl->k_edge_keV) : 0.0;
            const double l_bind = flL ? static_cast<double>(flL->l3_edge_keV) : 0.0;
            auto relax_selected = [&](bool is_L, int l_subshell,
                                      double trans_keV) {
                if (trans_keV > 0.0)
                    D += ic_e(trans_keV, is_L ? l_bind : k_bind);
                if (is_L) {
                    D += emit_l_xray(flL, l_subshell);
                    return;
                }
                const double ex = sample_vacancy_xray(fl, uniform, rng);
                if (ex <= 0.0) {
                    D += ic_e(k_bind, 0.0);
                    return;
                }
                D += transport_cascade_member(
                    vertex, sample_isotropic_dir(uniform, rng), ex, seg_buffer, rng);
                if (fl && ex < 0.93 * fl->k_edge_keV)
                    D += emit_l_xray(flL, -1);
            };

            // EC is part of the global weak categorical law, never an
            // independent vacancy source.
            if (selected_ec_shell == 1) relax_selected(false, -1, 0.0);
            else if (selected_ec_shell == 2) relax_selected(true, -1, 0.0);

            // One draw per conversion transition: none/K/L1/L2/L3/outer are
            // mutually exclusive.  For ordinary IC the draw is gated by the
            // complementary gamma outcome; a memberless E0 is unconditional.
            for (const VacancyGroup& vg : dc.vacancy_groups) {
                if (vg.kind == VacancyGroupKind::InternalConversion &&
                    vg.gamma_member >= 0 &&
                    vg.gamma_member < static_cast<int>(nm) &&
                    emitted[vg.gamma_member])
                    continue;
                const double uv = uniform(rng);
                double av = vg.p_none;
                if (uv < av) continue;
                av += vg.p_K;
                if (uv < av) { relax_selected(false, -1, vg.transition_energy_keV); continue; }
                av += vg.p_L1;
                if (uv < av) { relax_selected(true, 0, vg.transition_energy_keV); continue; }
                av += vg.p_L2;
                if (uv < av) { relax_selected(true, 1, vg.transition_energy_keV); continue; }
                av += vg.p_L3;
                if (uv < av) { relax_selected(true, 2, vg.transition_energy_keV); continue; }
                av += vg.p_outer;
                if (uv < av) {
                    D += ic_e(vg.transition_energy_keV, 0.0);
                    continue;
                }
                av += vg.p_unmodeled;
                if (uv < av)
                    continue;  // explicit unmodelled above-pair E0 outcome
            }
            for (const KVacancySource& kv : dc.k_vacancies) {
                // An IC vacancy can occur only when its gamma did NOT emit
                // (the transition either emits the gamma or converts).
                if (kv.gamma_member >= 0 &&
                    kv.gamma_member < static_cast<int>(nm) &&
                    emitted[kv.gamma_member])
                    continue;
                if (uniform(rng) >= kv.prob) continue;  // no vacancy
                // Internal-conversion electron (only for an IC vacancy; an EC
                // vacancy has gamma_member == -1 and produces no fast electron).
                if (kv.gamma_member >= 0 && kv.gamma_member < static_cast<int>(nm))
                    D += ic_e(dc.members[kv.gamma_member].energy_keV,
                              kv.is_L ? l_bind : k_bind);
                if (kv.is_L) {
                    D += emit_l_xray(flL, kv.l_subshell);
                    continue;
                }
                const double ex = sample_vacancy_xray(fl, uniform, rng);
                if (ex <= 0.0) {  // K-Auger: deposit the K binding locally.
                    D += ic_e(k_bind, 0.0);
                    continue;
                }
                D += transport_cascade_member(
                    vertex, sample_isotropic_dir(uniform, rng), ex, seg_buffer, rng);
                // K x-ray from L (Kalpha) leaves an L vacancy -> L x-ray.
                if (fl && ex < 0.93 * fl->k_edge_keV)
                    D += emit_l_xray(flL, -1);
            }
        }
      }  // end else (member-realization + independent vacancy model)

        if (!tally.spectrum.empty() && D > 0.0) {
            auto it = std::upper_bound(spectrum_edges.begin(), spectrum_edges.end(),
                                       static_cast<float>(D));
            if (it != spectrum_edges.begin()) {
                const std::size_t bin =
                    static_cast<std::size_t>(it - spectrum_edges.begin()) - 1;
                if (bin < tally.spectrum.size()) tally.spectrum[bin] += 1.0;
            }
        }

        for (std::size_t pi = 0; pi < targets.size(); ++pi) {
            const CascadePeakTarget& tg = targets[pi];
            if (!tg.found) continue;
            if (std::abs(D - tg.energy_keV) < tg.tol_keV) ++tally.n_sum[pi];
            // A fitted window can contain the same line from several mutually
            // exclusive branches, or multiple unresolved gamma identities in
            // one branch. The denominator is the total emitted-line yield, not
            // the yield of only the dominant identity. Allow >1 count/history
            // when two matching photons truly co-emit: the unsummed yield then
            // contains two photons while the pulse-height numerator remains one
            // event.
            for (const CascadePeakTarget::PrimaryIdentity& id
                 : tg.primary_matches) {
                if (id.branch != b || id.member >= nm || !emitted[id.member])
                    continue;
                ++tally.n_emit[pi];
                if (std::abs(dep[id.member] - tg.energy_keV) < tg.tol_keV)
                    ++tally.n_nosum[pi];
            }
        }

        // RNG-neutral progress flush (reads local counters + clock only).
        if (progress_flush && (e & kCheckMask) == kCheckMask) {
            const auto now = std::chrono::steady_clock::now();
            if (std::chrono::duration<double>(now - last_flush).count() >= kFlushSec) {
                last_flush = now;
                progress_flush(tally);
            }
        }
    }
    return tally;
}

std::vector<EfficiencyCalculator::CascadePeakTarget>
EfficiencyCalculator::cascade_locate_targets(const CascadeConfig& config) const
{
    std::vector<CascadePeakTarget> targets;
    targets.reserve(config.peaks.size());
    for (const PeakWindow& peak : config.peaks) {
        CascadePeakTarget tg;
        tg.energy_keV = peak.energy_keV;
        tg.tol_keV = peak.tolerance_keV;
        double best = -1.0;
        for (std::size_t b = 0; b < config.cascades.size(); ++b)
            for (std::size_t mi = 0; mi < config.cascades[b].members.size(); ++mi) {
                const CascadeMember& mem = config.cascades[b].members[mi];
                if (mem.type != CascadeParticleType::Gamma) continue;
                if (std::abs(mem.energy_keV - peak.energy_keV) <= peak.tolerance_keV) {
                    const DecayCascade& dc = config.cascades[b];
                    double emission = std::clamp(mem.intensity, 0.0, 1.0);
                    const LevelDag dag(dc.level_scheme);
                    if (dag.valid) {
                        const int t = dag.transition_of(static_cast<int>(mi));
                        if (t >= 0)
                            emission = weak_pass_probability(dc, dag, t) *
                                       dag.ts[t].p_gamma;
                        else if (const int r = residual_transition_of(
                                     dc, static_cast<int>(mi)); r >= 0)
                            emission = std::clamp(dc.residual_transitions[
                                static_cast<std::size_t>(r)].p_gamma, 0.0, 1.0);
                    }
                    const double score = std::max(0.0, dc.branch_weight) * emission;
                    if (!(score > 0.0)) continue;
                    tg.primary_matches.push_back({b, mi});
                    if (score > best) {
                        best = score; tg.branch = b; tg.member = mi; tg.found = true;
                    }
                }
            }
        if (tg.found) {
            for (const CascadePeakTarget::PrimaryIdentity& id
                 : tg.primary_matches)
                if (!config.cascades[id.branch].level_scheme.valid)
                    tg.summing_model_complete = false;
            for (std::size_t b = 0; b < config.cascades.size(); ++b) {
                const int excluded = b == tg.branch
                    ? static_cast<int>(tg.member) : -1;
                if (cascade_invalid_branch_can_feed_peak(
                        config.cascades[b], peak.energy_keV,
                        peak.tolerance_keV, excluded))
                    tg.summing_model_complete = false;
            }
        }
        targets.push_back(tg);
    }
    return targets;
}

PeakCascadeResult EfficiencyCalculator::full_peak_result(
    uint64_t n_emit, uint64_t n_nosum, uint64_t n_sum,
    double energy_keV, bool found)
{
    PeakCascadeResult pr;
    pr.energy_keV = energy_keV;
    pr.found = found;
    const double ne = static_cast<double>(n_emit);
    if (found && ne > 0.0) {
        const double es = n_sum   / ne;  // ε with summing
        const double en = n_nosum / ne;  // ε no summing (single-gamma FEP)
        pr.eff_with_summing = es;
        pr.eff_no_summing   = en;
        pr.eff_with_summing_unc = es * std::sqrt(
            (n_sum > 0 ? 1.0 / n_sum : 0.0) + 1.0 / ne);
        pr.eff_no_summing_unc = (n_nosum > 0)
            ? en * std::sqrt(1.0 / n_nosum + 1.0 / ne) : 0.0;
        if (en > 0.0) {
            pr.summing_factor = es / en;
            const double rel = (n_sum > 0 ? 1.0 / n_sum : 0.0)
                             + (n_nosum > 0 ? 1.0 / n_nosum : 0.0);
            pr.summing_factor_unc = pr.summing_factor * std::sqrt(rel);
        }
    }
    return pr;
}

PeakCascadeResult EfficiencyCalculator::conditional_peak_result(
    const CascadePeakTally& agg, double energy_keV, bool found, double prim_w)
{
    PeakCascadeResult pr;
    pr.energy_keV = energy_keV;
    pr.found = found;
    const double N = static_cast<double>(agg.n);
    if (found && N > 0.0) {
        const double p_no = agg.n_nosum / N;
        const double xbar = agg.sum_x / N;                 // E[X], with-summing score
        const double var_x = std::max(0.0, agg.sum_xx / N - xbar * xbar);
        pr.eff_no_summing       = prim_w * p_no;
        pr.eff_with_summing     = prim_w * xbar;
        pr.eff_no_summing_unc   = prim_w * std::sqrt(std::max(0.0, p_no * (1.0 - p_no) / N));
        pr.eff_with_summing_unc = prim_w * std::sqrt(var_x / N);
        if (p_no > 0.0 && xbar > 0.0) {
            const double k = xbar / p_no;
            pr.summing_factor = k;
            // Correlated ratio variance (numerator/denominator share events).
            // Reduces to the old integer-count formula when X in {0,1}.
            const double cov_xd = agg.sum_xd / N - xbar * p_no;
            const double term = var_x / (N * xbar * xbar)
                              + (1.0 - p_no) / (N * p_no)
                              - 2.0 * cov_xd / (N * xbar * p_no);
            pr.summing_factor_unc = std::sqrt(std::max(0.0, k * k * term));
        } else {
            pr.summing_factor = (p_no > 0.0) ? 0.0 : 1.0;
            pr.summing_factor_unc = 0.0;
        }
    }
    return pr;
}

CascadeResult EfficiencyCalculator::compute_cascade(const CascadeConfig& config) const
{
    validate_cascade_member_references(
        config.cascades, "EfficiencyCalculator::compute_cascade");
    CascadeResult result;
    result.num_events_per_peak = config.num_events;

    unsigned nthreads = config.num_threads ? config.num_threads
                                           : std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;

    // Map each requested peak to its dominant (branch, gamma member); shared
    // by both methods.
    const std::vector<CascadePeakTarget> targets = cascade_locate_targets(config);

    // Progress bookkeeping. Only exercised when a callback is set; the flush
    // lambdas below never touch the RNG, so results are bit-identical with or
    // without a callback.
    const auto t_start = std::chrono::steady_clock::now();
    std::mutex prog_mutex;
    auto last_fire = t_start;  // first intermediate fire ~1 s in

    // ---- Full-realization estimator (all summing channels + summed spectrum) ----
    if (config.method == CascadeMethod::FullRealization) {
        // Branch selection CDF (∝ branch_weight) and R = Σ branch_weight.
        double R = 0.0;
        for (const DecayCascade& dc : config.cascades)
            R += std::max(0.0, dc.branch_weight);
        if (!(R > 0.0))
            throw std::invalid_argument(
                "EfficiencyCalculator::compute_cascade: FullRealization "
                "requires at least one positive branch weight");
        if (config.num_events > 0 && config.spectrum_bin_edges.size() >= 2)
            for (const DecayCascade& dc : config.cascades)
                if (dc.branch_weight > 0.0 && !dc.level_scheme.valid)
                    result.summed_spectrum_model_complete = false;
        std::vector<double> cdf(config.cascades.size(), 1.0);
        double acc = 0.0;
        for (std::size_t b = 0; b < config.cascades.size(); ++b) {
            acc += std::max(0.0, config.cascades[b].branch_weight);
            cdf[b] = acc / R;
        }

        // Per-thread progress slots (counters only — never the spectrum, so a
        // flush stays O(#peaks) even with fine spectrum bins).
        std::vector<CascadeFullTally> slots;
        if (config.progress_callback) {
            slots.resize(nthreads);
            for (CascadeFullTally& s : slots) {
                s.n_emit.assign(targets.size(), 0);
                s.n_nosum.assign(targets.size(), 0);
                s.n_sum.assign(targets.size(), 0);
            }
        }

        const uint64_t per = config.num_events / nthreads;
        std::vector<std::future<CascadeFullTally>> futures;
        futures.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            const uint64_t ne = per + (t == nthreads - 1
                ? config.num_events - per * nthreads : 0);
            const uint64_t seed = 0x9E3779B97F4A7C15ULL
                ^ (static_cast<uint64_t>(t) * 0xD1B54A32D192ED03ULL);
            const std::vector<DecayCascade>& casc = config.cascades;
            const std::vector<float>& edges = config.spectrum_bin_edges;
            std::function<void(const CascadeFullTally&)> flush;
            if (config.progress_callback) {
                flush = [&, t](const CascadeFullTally& local) {
                    std::lock_guard<std::mutex> lock(prog_mutex);
                    slots[t].n = local.n;              // counters only, no spectrum
                    slots[t].n_emit  = local.n_emit;
                    slots[t].n_nosum = local.n_nosum;
                    slots[t].n_sum   = local.n_sum;
                    const auto now = std::chrono::steady_clock::now();
                    if (std::chrono::duration<double>(now - last_fire).count() < 1.0)
                        return;
                    last_fire = now;
                    CascadeProgress p;
                    for (const CascadeFullTally& s : slots) p.num_events += s.n;
                    p.elapsed = now - t_start;
                    p.frac_complete = std::min(1.0,
                        static_cast<double>(p.num_events) /
                        static_cast<double>(std::max<uint64_t>(1, config.num_events)));
                    p.peaks.reserve(targets.size());
                    for (std::size_t pi = 0; pi < targets.size(); ++pi) {
                        uint64_t n_emit = 0, n_nosum = 0, n_sum = 0;
                        for (const CascadeFullTally& s : slots) {
                            n_emit += s.n_emit[pi]; n_nosum += s.n_nosum[pi];
                            n_sum += s.n_sum[pi];
                        }
                        p.peaks.push_back(full_peak_result(n_emit, n_nosum, n_sum,
                            targets[pi].energy_keV, targets[pi].found));
                        p.peaks.back().summing_model_complete =
                            targets[pi].summing_model_complete;
                    }
                    config.progress_callback(p);
                };
            }
            const bool ic_e = config.enable_ic_electrons;
            futures.push_back(std::async(std::launch::async,
                [this, &casc, cdf, R, &targets, &edges, ne, seed, ic_e, flush] {
                    return cascade_full_thread(casc, cdf, R, targets, edges, ne, seed,
                                               ic_e, flush);
                }));
        }

        CascadeFullTally agg;
        agg.n_emit.assign(targets.size(), 0);
        agg.n_nosum.assign(targets.size(), 0);
        agg.n_sum.assign(targets.size(), 0);
        if (config.spectrum_bin_edges.size() >= 2)
            agg.spectrum.assign(config.spectrum_bin_edges.size() - 1, 0.0);
        for (auto& f : futures) {
            const CascadeFullTally t = f.get();
            agg.n += t.n;
            for (std::size_t i = 0; i < targets.size(); ++i) {
                agg.n_emit[i]  += t.n_emit[i];
                agg.n_nosum[i] += t.n_nosum[i];
                agg.n_sum[i]   += t.n_sum[i];
            }
            for (std::size_t bn = 0; bn < agg.spectrum.size() && bn < t.spectrum.size(); ++bn)
                agg.spectrum[bn] += t.spectrum[bn];
        }

        for (std::size_t pi = 0; pi < targets.size(); ++pi)
            result.peaks.push_back(full_peak_result(
                agg.n_emit[pi], agg.n_nosum[pi], agg.n_sum[pi],
                targets[pi].energy_keV, targets[pi].found));
        for (std::size_t pi = 0; pi < targets.size(); ++pi)
            result.peaks[pi].summing_model_complete =
                targets[pi].summing_model_complete;

        // Summed spectrum, normalised to counts per parent decay (R · count / N).
        if (!agg.spectrum.empty() && agg.n > 0) {
            const double scale = R / static_cast<double>(agg.n);
            result.summed_spectrum.resize(agg.spectrum.size());
            result.summed_spectrum_uncertainty.resize(agg.spectrum.size());
            for (std::size_t bn = 0; bn < agg.spectrum.size(); ++bn) {
                result.summed_spectrum[bn] =
                    static_cast<float>(scale * agg.spectrum[bn]);
                result.summed_spectrum_uncertainty[bn] =
                    static_cast<float>(scale * std::sqrt(agg.spectrum[bn]));
            }
        }

        // Completion fire: exactly once, peaks identical to the returned result.
        if (config.progress_callback) {
            CascadeProgress p;
            p.num_events = agg.n;
            p.elapsed = std::chrono::steady_clock::now() - t_start;
            p.frac_complete = 1.0;
            p.peaks = result.peaks;
            config.progress_callback(p);
        }
        return result;
    }

    // Cone for the primary gamma. Only valid for a fixed point vertex with no
    // source scatter; extended/shielded sources emit isotropically (a scattered
    // photon can reach the detector from outside any cone, so cone biasing would
    // be biased there).
    double cone_half = kPi;
    if (use_cone_sampling_ && source_type_ == SourceType::Point
        && !source_geometry_.has_source_effects())
        cone_half = compute_cone_half_angle(source_position_);
    const bool coned = (cone_half < kPi - 0.01);
    const double cos_max = std::cos(cone_half);
    const double prim_w = coned ? 0.5 * (1.0 - cos_max) : 1.0;  // Omega/4pi

    // Conditional progress bookkeeping: a persistent running snapshot across
    // the sequential per-peak fan-outs (completed peaks keep their final
    // values; the in-flight peak carries its running estimate).
    std::vector<PeakCascadeResult> snap;
    uint64_t total_events = 0;  // num_events × found peaks
    uint64_t events_done = 0;   // histories from completed peaks
    std::vector<CascadePeakTally> slots;
    if (config.progress_callback) {
        snap.resize(config.peaks.size());
        for (std::size_t i = 0; i < snap.size(); ++i) {
            snap[i].energy_keV = config.peaks[i].energy_keV;
            snap[i].found = targets[i].found;
            snap[i].summing_model_complete =
                targets[i].summing_model_complete;
        }
        for (const CascadePeakTarget& tg : targets)
            if (tg.found) total_events += config.num_events;
        slots.resize(nthreads);
    }

    std::size_t peak_index = 0;
    for (const PeakWindow& peak : config.peaks) {
        const CascadePeakTarget& tg = targets[peak_index];
        if (!tg.found) {
            PeakCascadeResult pr;
            pr.energy_keV = peak.energy_keV;
            pr.found = false;
            result.peaks.push_back(pr);
            ++peak_index;
            continue;
        }

        // A single conditioned stream cannot represent an unresolved window
        // with several positive-yield primary identities without either
        // repeating sum-fed channels or omitting branch-specific summing-out.
        // Use the exact full-history estimator for this uncommon case. The
        // subcall has no spectrum or callback, so recursion terminates and the
        // outer callback retains its documented one-completion behavior.
        if (tg.primary_matches.size() > 1) {
            CascadeConfig full_cfg = config;
            full_cfg.peaks = {peak};
            full_cfg.method = CascadeMethod::FullRealization;
            full_cfg.spectrum_bin_edges.clear();
            full_cfg.progress_callback = {};
            const CascadeResult full = compute_cascade(full_cfg);
            PeakCascadeResult pr = full.peaks.front();
            if (config.progress_callback) {
                std::lock_guard<std::mutex> lock(prog_mutex);
                snap[peak_index] = pr;
                events_done += config.num_events;
                for (CascadePeakTally& s : slots) s = CascadePeakTally{};
            }
            result.peaks.push_back(pr);
            ++peak_index;
            continue;
        }
        const std::size_t best_branch = tg.branch;
        const std::size_t best_mi = tg.member;

        // Analytic sum-peak-fed pair channels (pairs a+b that deposit into this
        // window without the peak's own gamma) and the primary's coincident
        // vacancy x-rays. Immutable, computed once, shared across threads.
        const std::vector<CascadeSumPairChannel> channels =
            cascade_sum_pair_channels(config.cascades, best_branch, best_mi,
                                      peak.energy_keV, peak.tolerance_keV);
        const std::vector<CascadeVacancyDraw> prim_vacs =
            cascade_level_vacancies(config.cascades[best_branch], best_mi);

        const uint64_t per = config.num_events / nthreads;
        std::vector<std::future<CascadePeakTally>> futures;
        futures.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            const uint64_t ne = per + (t == nthreads - 1
                ? config.num_events - per * nthreads : 0);
            const uint64_t seed = 0x9E3779B97F4A7C15ULL
                ^ (static_cast<uint64_t>(peak_index + 1) << 32)
                ^ (static_cast<uint64_t>(t) * 0xD1B54A32D192ED03ULL);
            const std::size_t bp = best_branch, mip = best_mi;
            const double pkv = peak.energy_keV, tolv = peak.tolerance_keV;
            std::function<void(const CascadePeakTally&)> flush;
            if (config.progress_callback) {
                flush = [&, t, peak_index, pkv](const CascadePeakTally& local) {
                    std::lock_guard<std::mutex> lock(prog_mutex);
                    slots[t] = local;
                    const auto now = std::chrono::steady_clock::now();
                    if (std::chrono::duration<double>(now - last_fire).count() < 1.0)
                        return;
                    last_fire = now;
                    CascadePeakTally cur;
                    for (const CascadePeakTally& s : slots) {
                        cur.n += s.n; cur.n_nosum += s.n_nosum;
                        cur.sum_x += s.sum_x; cur.sum_xx += s.sum_xx;
                        cur.sum_xd += s.sum_xd;
                    }
                    snap[peak_index] = conditional_peak_result(cur, pkv, true, prim_w);
                    snap[peak_index].summing_model_complete =
                        targets[peak_index].summing_model_complete;
                    CascadeProgress p;
                    p.num_events = events_done + cur.n;
                    p.elapsed = now - t_start;
                    p.frac_complete = total_events
                        ? std::min(1.0, static_cast<double>(p.num_events) /
                                        static_cast<double>(total_events))
                        : 1.0;
                    p.peaks = snap;
                    config.progress_callback(p);
                };
            }
            const bool ic_e = config.enable_ic_electrons;
            futures.push_back(std::async(std::launch::async,
                [this, &config, bp, mip, pkv, tolv, cos_max, coned,
                 &prim_vacs, &channels, ne, seed, ic_e, flush] {
                    return cascade_peak_thread(config.cascades, bp, mip, pkv, tolv,
                                               cos_max, coned, prim_vacs, channels,
                                               ne, seed, ic_e, flush);
                }));
        }

        CascadePeakTally agg;
        for (auto& f : futures) {
            const CascadePeakTally t = f.get();
            agg.n += t.n; agg.n_nosum += t.n_nosum;
            agg.sum_x += t.sum_x; agg.sum_xx += t.sum_xx; agg.sum_xd += t.sum_xd;
        }

        PeakCascadeResult pr =
            conditional_peak_result(agg, peak.energy_keV, true, prim_w);
        pr.summing_model_complete = tg.summing_model_complete;
        if (config.progress_callback) {
            // Peak boundary: record the final values and reset the per-thread
            // slots for the next fan-out (no forced fire — throttle governs).
            std::lock_guard<std::mutex> lock(prog_mutex);
            snap[peak_index] = pr;
            events_done += agg.n;
            for (CascadePeakTally& s : slots) s = CascadePeakTally{};
        }
        result.peaks.push_back(pr);
        ++peak_index;
    }

    // Completion fire: exactly once, peaks identical to the returned result.
    if (config.progress_callback) {
        CascadeProgress p;
        p.num_events = events_done;
        p.elapsed = std::chrono::steady_clock::now() - t_start;
        p.frac_complete = 1.0;
        p.peaks = result.peaks;
        config.progress_callback(p);
    }
    return result;
}

// --- Main compute method (SimulationConfig) ---

double EfficiencyCalculator::central_ray_optical_depth(double energy_keV) const
{
    if (source_type_ == SourceType::Marinelli) {
        return std::numeric_limits<double>::infinity();
    }

    Eigen::Vector3d from_pos = source_position_;
    if (source_type_ == SourceType::Cylindrical) from_pos = cyl_src_.center;
    if (source_type_ == SourceType::Rectangular) from_pos = rect_src_.center;

    if (from_pos.norm() < 1e-9) {
        return std::numeric_limits<double>::infinity();
    }
    Eigen::Vector3d dir = -from_pos.normalized();

    auto segments = geometry_.trace_ray(from_pos, dir);
    if (segments.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const double energy_MeV = energy_keV * 1e-3;
    double tau = 0.0;
    for (const auto& seg : segments) {
        if (!seg.material) continue;
        double len = seg.t_end - std::max(seg.t_start, 0.0);
        if (len <= 0.0) continue;
        MacroscopicXS xs = seg.material->macroscopic_xs(energy_MeV);
        if (!transport_config_.enable_rayleigh) xs.mu_rs = 0.0;
        if (!transport_config_.enable_pair_production || energy_keV < 1022.0) {
            xs.mu_pp = 0.0;
        }
        tau += xs.mu_total() * len;
    }
    return tau;
}

BiasingConfig EfficiencyCalculator::compute_effective_biasing(
    const SimulationConfig& config) const
{
    if (biasing_manual_) {
        return biasing_;
    }

    BiasingConfig b;
    const auto& term = config.termination;

    // FEP precision matters unless the run targets ONLY total efficiency.
    const bool fep_matters = (term.target_fep_rel_precision > 0.0) ||
                             (term.target_total_rel_precision <= 0.0);

    // Forced first interaction: wins when the analog pass-through
    // probability through the detector is high (thin/low-tau geometries);
    // in thick crystals it improves total-efficiency FOM but degrades FEP
    // FOM (every event pays full transport cost).
    const double tau = central_ray_optical_depth(config.energy_keV);
    const double tau_threshold = fep_matters ? 0.7 : 4.0;
    b.force_detector_interaction = (tau < tau_threshold);

    // Two-stream stratified estimator (+ Compton-angle biasing) for
    // full-mode source-effect runs. Supersedes the mixture angular biasing
    // it replaced (mixture traded total-efficiency FOM x0.04-0.6 for its
    // FEP gains; the two-stream split improves both: measured joint
    // min-FOM x1.6-17 vs analog across the benchmark + stress sets).
    // Excluded for Marinelli (cone degenerate for most positions, re-entry
    // physics) and for close geometries where the detector subtends a
    // large solid angle from the source center (omega/4pi >= 0.15;
    // measured neutral-to-negative there, e.g. a water puck 1 cm from the
    // face, while analog is already efficient).
    if (source_geometry_.has_source_effects() && !fep_only_mode_ &&
        source_type_ != SourceType::Marinelli) {
        Eigen::Vector3d ref_pos = source_position_;
        if (source_type_ == SourceType::Cylindrical) ref_pos = cyl_src_.center;
        if (source_type_ == SourceType::Rectangular) ref_pos = rect_src_.center;
        const double cone_half = compute_cone_half_angle(ref_pos);
        const double omega_frac = 0.5 * (1.0 - std::cos(cone_half));
        if (cone_half < kPi - 0.01 && omega_frac < 0.15) {
            const bool total_or_spectrum =
                (term.target_total_rel_precision > 0.0) ||
                !config.energy_bin_edges.empty() ||
                (term.target_fep_rel_precision <= 0.0);
            b.two_stream = true;
            // f = 0.25 balances FEP and total (total FOM stays >= analog);
            // f = 0.5 when only FEP precision is targeted. gamma = 0.3 is
            // the measured safe optimum for the scattered-in total term
            // (gamma >= 0.5 degrades large unshielded sources); it only
            // helps totals, so it is off for FEP-only targets.
            b.two_stream_direct_fraction = total_or_spectrum ? 0.25 : 0.5;
            b.compton_cone_fraction = total_or_spectrum ? 0.3 : 0.0;
        }
    }
    return b;
}

EfficiencyResult EfficiencyCalculator::compute(const SimulationConfig& config)
{
    auto t_start = std::chrono::steady_clock::now();

    // Resolve event-biasing settings for this run (manual or auto-enabled).
    active_biasing_ = compute_effective_biasing(config);

    unsigned num_threads = config.num_threads;
    if (num_threads == 0) {
        num_threads = std::max(1u, std::thread::hardware_concurrency());
    }

    const uint64_t batch_size = std::max(config.batch_size, uint64_t(100));
    const auto& term = config.termination;

    // For very small max_events, use single thread
    if (term.max_events > 0 && term.max_events < 1000) {
        num_threads = 1;
    }

    // Set up shared state
    SimulationState state;
    state.termination = term;
    state.progress_callback = config.progress_callback;
    state.start_time = t_start;
    state.last_callback_time = t_start;  // first throttled fire ~1 s in
    if (config.energy_bin_edges.size() >= 2) {
        state.merged.pulse_height_counts.resize(config.energy_bin_edges.size() - 1, 0.0);
        state.merged.pulse_height_counts_sq.resize(config.energy_bin_edges.size() - 1, 0.0);
    }

    // Pass detector bounding geometry to source for electron filtering and
    // for Compton-angle biasing (the bias cone subtends these bounds).
    if ((source_electron_transport_ ||
         active_biasing_.compton_cone_fraction > 0.0) &&
        source_geometry_.has_source_effects()) {
        source_geometry_.set_detector_bounds(
            geometry_.outer_bounding_radius(),
            geometry_.outer_z_extent().first,
            geometry_.outer_z_extent().second);
    }

    // Configure (or reset) Compton-angle biasing for this run.
    {
        SourceGeometry::ComptonBiasConfig cb;
        cb.cone_fraction =
            std::clamp(active_biasing_.compton_cone_fraction, 0.0, 0.9);
        source_geometry_.set_compton_bias(cb);
    }

    // Seed generator. A nonzero config.seed makes the run reproducible
    // (deterministic per-thread seeds); 0 falls back to std::random_device.
    std::mt19937_64 seed_rng;
    if (config.seed != 0) {
        seed_rng.seed(config.seed);
    } else {
        std::random_device rd;
        seed_rng.seed(rd());
    }

    // Thread function: run batches until stop
    auto thread_fn = [&](uint64_t initial_seed) {
        const uint64_t is_mixed = initial_seed * 2654435761ULL;
        const std::array<uint_least32_t, 6> is_data{{
          static_cast<uint_least32_t>(initial_seed), static_cast<uint_least32_t>(initial_seed >> 32),
          static_cast<uint_least32_t>(is_mixed), static_cast<uint_least32_t>(is_mixed >> 32),
          static_cast<uint_least32_t>(initial_seed ^ (initial_seed >> 16)),
          static_cast<uint_least32_t>((initial_seed >> 48) | (initial_seed << 16))
        }};
        std::seed_seq seq( is_data.begin(), is_data.end() );
        std::mt19937_64 rng(seq);

        double cpu_last = thread_cpu_seconds();

        // Geometric batch growth: the caller's batch_size sizes the FIRST
        // wave (threads x batch ~= min_events, so cheap runs stop right at
        // the precision check) and later waves double toward the classic
        // 10k so long runs are not taxed by per-batch overhead.
        uint64_t grow_batch = batch_size;
        const uint64_t batch_cap = std::max<uint64_t>(batch_size, 10000);

        while (!state.stop_flag.load(std::memory_order_acquire)) {
            // Determine batch size (don't overshoot max_events)
            uint64_t this_batch = grow_batch;
            grow_batch = std::min(grow_batch * 2, batch_cap);
            if (term.max_events > 0) {
                uint64_t done;
                { std::lock_guard<std::mutex> lock(state.mutex); done = state.merged.num_events; }
                if (done >= term.max_events) break;
                this_batch = std::min(this_batch, term.max_events - done);
            }
            if (this_batch == 0) break;

            // Run batch using existing simulate_thread with a fresh per-batch seed
            uint64_t batch_seed = rng();
            auto batch_tally = simulate_thread(
                config.energy_keV, this_batch, batch_seed, config.energy_bin_edges);

            // This thread's CPU cost since the previous merge (covers the batch
            // plus loop overhead; nothing is lost between batches).
            const double cpu_now = thread_cpu_seconds();
            const double batch_cpu_s = cpu_now - cpu_last;
            cpu_last = cpu_now;

            // Merge and check termination
            state.merge_and_check(batch_tally, batch_cpu_s);
        }
    };

    // Launch worker threads
    std::vector<std::thread> threads;
    for (unsigned i = 0; i < num_threads; ++i) {
        threads.emplace_back(thread_fn, seed_rng());
    }

    for (auto& t : threads) {
        t.join();
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    // Build result from merged tally
    EfficiencyResult result;
    result.num_events_simulated = state.merged.num_events;
    result.wall_time_seconds = elapsed;
    result.cpu_time_seconds = state.cpu_seconds;
    result.stop_reason = state.stop_reason;
    result.pp_secondaries_processed = state.merged.dbg_pp_secondary;
    result.pp_secondary_deposit_keV = state.merged.dbg_dep_pp;
    result.pp_only_any_weight = state.merged.dbg_pp_only_any_w;
    result.src_diag = state.merged.dbg_src;
    result.forced_absorption_fraction =
        state.merged.num_events > 0
            ? static_cast<double>(state.merged.num_forced_absorption) /
                  static_cast<double>(state.merged.num_events)
            : 0.0;

    {
        const TallyEstimates est = compute_estimates(state.merged);
        result.full_energy_peak_efficiency = est.eps_fep;
        result.fep_uncertainty = est.sig_fep;
        result.total_efficiency = est.eps_tot;
        result.total_uncertainty = est.sig_tot;
    }

    // Completion fire: exactly once, payload identical to the returned result.
    if (config.progress_callback) {
        Progress p;
        p.num_events = result.num_events_simulated;
        p.elapsed = t_end - t_start;
        p.frac_complete = 1.0;
        p.fep_efficiency = result.full_energy_peak_efficiency;
        p.fep_uncertainty = result.fep_uncertainty;
        p.total_efficiency = result.total_efficiency;
        p.total_uncertainty = result.total_uncertainty;
        config.progress_callback(p);
    }

    if (!state.merged.pulse_height_counts.empty() && state.merged.sum_weights > 0.0) {
        result.pulse_height_distribution.resize(state.merged.pulse_height_counts.size());
        result.pulse_height_uncertainty.resize(state.merged.pulse_height_counts.size());
        for (size_t i = 0; i < state.merged.pulse_height_counts.size(); ++i) {
            result.pulse_height_distribution[i] =
                static_cast<float>(state.merged.pulse_height_counts[i] / state.merged.sum_weights);
            // Per-bin 1-sigma: sqrt(Sum w^2) / sum_weights (sum_weights == N).
            // Analog (w==1) reduces to sqrt(count)/N == Poisson; weighted
            // sampling differs by the in-bin weight dispersion.
            result.pulse_height_uncertainty[i] =
                static_cast<float>(std::sqrt(state.merged.pulse_height_counts_sq[i]) /
                                   state.merged.sum_weights);
        }
    }

    // Print Marinelli re-entry diagnostics if any data was collected
    if (source_type_ == SourceType::Marinelli &&
        state.merged.dbg_initial_hit > 0) {
        auto& m = state.merged;
        double N = static_cast<double>(m.num_events);
        std::fprintf(stderr,
            "\n=== Marinelli Re-entry Diagnostics ===\n"
            "Events: %llu\n"
            "Initial crystal hits:  %llu (%.4f%%)\n"
            "Miss-bounce recoveries: %llu (%.4f%%)\n"
            "Post-crystal re-entries: %llu (%.4f per initial hit)\n"
            "Escaped secondary hits:  %llu (%.4f per initial hit)\n"
            "PP annihilation gammas:  %llu\n"
            "\nEnergy deposits (keV sum / event):\n"
            "  Initial:     %.2f\n"
            "  Re-entry:    %.2f\n"
            "  Secondary:   %.2f\n"
            "  PP 511:      %.2f\n"
            "  Total:       %.2f\n"
            "\nRe-entry pipeline funnel:\n"
            "  Escaped crystal:       %llu\n"
            "  Can re-enter water:    %llu (%.1f%% of escaped)\n"
            "  Wall passed:           %llu\n"
            "  Wall scattered:        %llu\n"
            "  Wall absorbed:         %llu\n"
            "  Water survived:        %llu (%.1f%% of wall-passed)\n"
            "  Trace hit crystal:     %llu (%.1f%% of water-survived)\n"
            "  Trace missed crystal:  %llu\n"
            "===================================\n\n",
            (unsigned long long)m.num_events,
            (unsigned long long)m.dbg_initial_hit, 100.0 * m.dbg_initial_hit / N,
            (unsigned long long)m.dbg_miss_bounce_hit, 100.0 * m.dbg_miss_bounce_hit / N,
            (unsigned long long)m.dbg_reentry_hit,
            m.dbg_initial_hit > 0 ? (double)m.dbg_reentry_hit / m.dbg_initial_hit : 0.0,
            (unsigned long long)m.dbg_secondary_hit,
            m.dbg_initial_hit > 0 ? (double)m.dbg_secondary_hit / m.dbg_initial_hit : 0.0,
            (unsigned long long)m.dbg_pp_secondary,
            m.dbg_dep_initial / N,
            m.dbg_dep_reentry / N,
            m.dbg_dep_secondary / N,
            m.dbg_dep_pp / N,
            (m.dbg_dep_initial + m.dbg_dep_reentry + m.dbg_dep_secondary + m.dbg_dep_pp) / N,
            (unsigned long long)m.dbg_n_escaped,
            (unsigned long long)m.dbg_n_can_reenter,
            m.dbg_n_escaped > 0 ? 100.0 * m.dbg_n_can_reenter / m.dbg_n_escaped : 0.0,
            (unsigned long long)m.dbg_n_wall_pass,
            (unsigned long long)m.dbg_n_wall_scatter,
            (unsigned long long)m.dbg_n_wall_absorb,
            (unsigned long long)m.dbg_n_water_survived,
            m.dbg_n_wall_pass > 0 ? 100.0 * m.dbg_n_water_survived / m.dbg_n_wall_pass : 0.0,
            (unsigned long long)m.dbg_n_trace_hit,
            m.dbg_n_water_survived > 0 ? 100.0 * m.dbg_n_trace_hit / m.dbg_n_water_survived : 0.0,
            (unsigned long long)m.dbg_n_trace_miss);
    }

    return result;
}

// --- Old compute() wrapper ---

EfficiencyResult EfficiencyCalculator::compute(
    double energy_keV,
    uint64_t num_events,
    unsigned num_threads,
    const std::vector<float>& energy_bin_edges)
{
    SimulationConfig config;
    config.energy_keV = energy_keV;
    config.termination.max_events = num_events;
    config.num_threads = num_threads;
    config.energy_bin_edges = energy_bin_edges;
    config.batch_size = std::max(uint64_t(10000),
                                 num_events / 10); // check termination ~10x per run
    return compute(config);
}

std::vector<EfficiencyResult> EfficiencyCalculator::compute_batch(
    const std::vector<double>& energies_keV,
    uint64_t num_events_per_energy,
    unsigned num_threads,
    const std::vector<float>& energy_bin_edges)
{
    std::vector<EfficiencyResult> results;
    results.reserve(energies_keV.size());
    for (double e : energies_keV) {
        results.push_back(compute(e, num_events_per_energy, num_threads, energy_bin_edges));
    }
    return results;
}

// --- Cascade correction ---

namespace {
/// Find the cache entry whose energy is closest to E within `tol` keV, or
/// nullptr if none qualifies. The cache is a std::map sorted by energy, so the
/// only candidates are the two keys bracketing E (lower_bound and its
/// predecessor). This tolerance match replaces the previous bit-exact `find`
/// so a rounded/round-tripped energy key does not silently miss (review C1).
const EfficiencyResult* find_cached_efficiency(
        const std::map<double, EfficiencyResult>& cache, double E, double tol)
{
    if (cache.empty())
        return nullptr;
    const EfficiencyResult* best = nullptr;
    double best_d = tol;
    auto hi = cache.lower_bound(E);            // first key >= E
    if (hi != cache.end()) {
        const double d = std::abs(hi->first - E);
        if (d <= best_d) { best_d = d; best = &hi->second; }
    }
    if (hi != cache.begin()) {
        auto lo = std::prev(hi);               // last key < E
        const double d = std::abs(lo->first - E);
        if (d <= best_d) { best_d = d; best = &lo->second; }
    }
    return best;
}
} // namespace

CascadeCorrectionResult EfficiencyCalculator::cascade_correction(
    double primary_fep,
    double primary_total,
    const std::vector<CoincidentGamma>& coincident,
    const std::map<double, EfficiencyResult>& efficiency_cache,
    double primary_fep_uncertainty,
    double match_tolerance_keV)
{
    return cascade_correction(primary_fep, primary_total, coincident,
                              efficiency_cache, {}, primary_fep_uncertainty,
                              match_tolerance_keV);
}

CascadeCorrectionResult EfficiencyCalculator::cascade_correction(
    double primary_fep,
    [[maybe_unused]] double primary_total,
    const std::vector<CoincidentGamma>& coincident,
    const std::map<double, EfficiencyResult>& efficiency_cache,
    const std::vector<SummingInPair>& summing_in_pairs,
    double primary_fep_uncertainty,
    double match_tolerance_keV)
{
    CascadeCorrectionResult res{};
    const double tol = std::max(match_tolerance_keV, 0.0);

    // NOTE: this is the flat pairwise-product form (spec Eq. 10). The C_out
    // product below is only exact when the supplied `coincidence_fraction`s are
    // marginal AND the coincident gammas are mutually independent; for gammas
    // that are mutually EXCLUSIVE de-excitations of a shared level it over-sums
    // (multiplies survival factors that cannot co-occur). The DecayCascade-driven
    // path compute_cascade_analytic() (cascade/AnalyticCascade.h) uses the exact
    // level-scheme survival DP instead (Eq. 10'), which zeroes exclusive partners
    // and adds vacancy-x-ray / sum-fed channels; when InterSpec feeds this
    // function DP-derived coincidence fractions the flat product becomes correct.

    // --- Summing-out: C_out = Π_j (1 - f_j ε_total(E_j)) ---
    // Retain (f, t=ε_total, σ_t) per matched gamma for the uncertainty below.
    double c_out = 1.0;
    struct OutTerm { double f, t, sig_t; };
    std::vector<OutTerm> out_terms;
    out_terms.reserve(coincident.size());
    for (const auto& cg : coincident) {
        const EfficiencyResult* e =
            find_cached_efficiency(efficiency_cache, cg.energy_keV, tol);
        if (e) {
            const double t = e->total_efficiency;
            c_out *= (1.0 - cg.coincidence_fraction * t);
            out_terms.push_back({cg.coincidence_fraction, t, e->total_uncertainty});
        } else {
            res.unmatched_energies.push_back(cg.energy_keV);
        }
    }

    // --- Summing-in: S = Σ f_ab ε_FEP(a) ε_FEP(b);  C_in = S / ε_FEP(primary) ---
    double sum_in = 0.0;         // absolute peak gain Σ f ε_a ε_b
    double var_sum_in = 0.0;     // its variance
    for (const auto& pair : summing_in_pairs) {
        const EfficiencyResult* ea =
            find_cached_efficiency(efficiency_cache, pair.energy_a_keV, tol);
        const EfficiencyResult* eb =
            find_cached_efficiency(efficiency_cache, pair.energy_b_keV, tol);
        if (ea && eb) {
            const double a = ea->full_energy_peak_efficiency;
            const double b = eb->full_energy_peak_efficiency;
            const double f = pair.joint_fraction;
            sum_in += f * a * b;
            const double sa = ea->fep_uncertainty;
            const double sb = eb->fep_uncertainty;
            var_sum_in += f * f * (b * b * sa * sa + a * a * sb * sb);
        } else {
            if (!ea) res.unmatched_energies.push_back(pair.energy_a_keV);
            if (!eb) res.unmatched_energies.push_back(pair.energy_b_keV);
        }
    }

    const bool have_primary = (primary_fep > 1e-15);
    const double c_in = have_primary ? (sum_in / primary_fep) : 0.0;
    const double c_net = c_out + c_in;

    res.summing_out_factor = c_out;
    res.summing_in_term    = c_in;
    res.net_correction     = c_net;
    res.corrected_fep      = primary_fep * c_net;  // identical to the prior formula

    // --- First-order uncertainty on corrected_fep (independent estimates) ---
    // corrected_fep = ε_FEP·C_out + Σ f ε_a ε_b  (the ε_FEP cancels in C_in), so
    //   d/dε_FEP = C_out ; d/dt_j = -ε_FEP·f_j·C_out/(1 - f_j t_j) ; ε_a,ε_b → var_sum_in.
    double var = (c_out * primary_fep_uncertainty) *
                 (c_out * primary_fep_uncertainty);
    for (const auto& ot : out_terms) {
        const double denom = 1.0 - ot.f * ot.t;
        if (std::abs(denom) > 1e-12) {
            const double d_dt = -primary_fep * ot.f * c_out / denom;
            var += d_dt * d_dt * ot.sig_t * ot.sig_t;
        }
    }
    if (have_primary)
        var += var_sum_in;   // the summing-in gain only enters corrected_fep here
    res.corrected_fep_uncertainty = std::sqrt(std::max(0.0, var));

    return res;
}

// --- GEANT4 Export ---

void EfficiencyCalculator::export_geant4_gdml(const std::string& filename,
                                               bool vacuum_world) const {
    const SourceGeometry* sg = source_geometry_.has_source_effects() ? &source_geometry_ : nullptr;
    write_gdml(geometry_, filename, sg, vacuum_world);
}

void EfficiencyCalculator::export_geant4_macro(const std::string& filename,
                                               double energy_keV,
                                               uint64_t num_events) const {
    if (source_type_ == SourceType::Marinelli) {
        write_geant4_macro_marinelli(source_geometry_, energy_keV, num_events, filename);
    } else if (source_type_ == SourceType::Cylindrical) {
        write_geant4_macro_cylindrical(
            source_geometry_.cyl_center(), source_geometry_.cyl_radius(),
            source_geometry_.cyl_half_length(),
            energy_keV, num_events, filename);
    } else if (source_type_ == SourceType::Rectangular) {
        write_geant4_macro_rectangular(
            source_geometry_.rect_center(), source_geometry_.rect_half_dims(),
            energy_keV, num_events, filename);
    } else if (source_type_ == SourceType::Spherical) {
        write_geant4_macro_spherical(
            source_geometry_.sphere_center(), source_geometry_.sphere_radius(),
            energy_keV, num_events, filename);
    } else {
        write_geant4_macro(source_position_, energy_keV, num_events, filename);
    }
}

} // namespace ceelo
