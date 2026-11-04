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

/// @file env_registered_e0_gate.cpp
/// @brief Held-out registered-E0 gate for the environment-scatter fit campaign.
///
/// Registered-E0 validation is a mandatory publication artifact. The Python
/// evidence consumer (`env_scatter_fit registered-e0-report`) runs this program
/// with --json-report and independently recomputes every acceptance metric
/// from the raw node maps it writes.
///
/// A bake's accuracy budgets apply AT its nodes, while interpolation error is
/// delivered to every consumer query that does not land exactly on a node.
/// Method: bake the candidate ladder AND its geometric midpoints, then compare
/// the midpoint truth against interpolation from the ladder nodes only.
///
///  1. **A matrix of interpolation schemes on the SAME baked nodes** -- x in
///     {axis, log axis} x y in {s_hat, ln s_hat} x {linear, PCHIP}. All eight
///     are pure post-processing of one bake, so choosing the representation
///     costs no extra MC. This matters because the curve is physically
///     two-regime on the distance axis (near field: s_hat ~ d^p as the fixed
///     ground geometry is normalized by a 1/d^2 direct term; far field: air
///     buildup, near-linear in optical depth), and the right variables can
///     linearize a regime that a raw-value spline has to bend through.
///     Note PCHIP is invariant under affine rescaling of x, so interpolating
///     in tau_air = mu_air * d is IDENTICAL to interpolating in d -- only the
///     log-vs-linear choice is a real degree of freedom.
///
///  2. **Per-node precision targeting.** Fixed histories leave far nodes at
///     ~1.3-1.5% while near nodes sit at 0.4%, so a 2-5% interpolation miss
///     out there is not resolvable (z ~ 1.4) and the gate cannot tell
///     "interpolation is fine" from "we cannot see it". Each node therefore
///     runs in batches until its relative sigma meets --precision.
///
/// Also reported: interpolation error per COARSE ENERGY BAND (mu-summed), the
/// shape metric -- a scheme that nails the integral while distorting the
/// spectrum is not acceptable for a kernel whose product is a spectrum.
///
///   env_registered_e0_gate [--axis d|h|e0] [--e0 662] [--precision 0.0025]
///                [--mode air|field --material soil|concrete]
///                [--absolute-precision 0.001]
///                [--absolute-min-frac 0.01]
///                [--native-energy-bins 384 --native-mu-bins 64]
///                [--batch 1000000] [--max-histories 40000000]
///                [--threads N] [--scoring-radius 4]
///                [--json-report FILE --config-sha256 HEX]

#include "geometry/EnvironmentGeometry.h"
#include "io/Pchip.h"
#include "materials/Material.h"
#include "transport/NextEventEstimator.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace ceelo;

namespace {

// Coarse energy bands for the shape metric (keV, lower edges + top).
const double kBandEdges[] = {30.0, 60.0, 120.0, 250.0, 500.0, 1000.0, 3000.0};
constexpr int kNBands = 6;

// Publication gate grid.  Energy bands are predefined ranges of registered
// channel INDICES, so every E0 contributes the same 12 coordinates even when a
// cell has a zero Monte-Carlo estimate.  The eight equal-width mu bands are a
// fixed, detector-frame angular partition with edges -1,-0.75,...,+1.
constexpr int kJointEnergyBands = 12;
constexpr int kJointMuBands = 8;
constexpr int kJointCells = kJointEnergyBands * kJointMuBands;

int joint_energy_band(int energy_bin, int n_energy)
{
    return std::min(kJointEnergyBands - 1,
                    energy_bin * kJointEnergyBands / n_energy);
}

int joint_mu_band(int mu_bin, int n_mu)
{
    return std::min(kJointMuBands - 1, mu_bin * kJointMuBands / n_mu);
}

int joint_index(int energy_band, int mu_band)
{
    return energy_band * kJointMuBands + mu_band;
}

bool is_lower_sha256(const std::string& value)
{
    if (value.size() != 64) return false;
    for (char c : value)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}

int band_of(double e_low_keV)
{
    for (int b = kNBands - 1; b >= 0; --b)
        if (e_low_keV >= kBandEdges[b]) return b;
    return 0;
}

/// How the shape bands are defined when comparing nodes of DIFFERENT E0.
///
/// Absolute-energy bands are meaningless across E0 (a 59.5 keV source fills
/// only [E_cut, 59.5] keV while a 2614 keV source fills everything), so the
/// bands must live in a coordinate that maps each node's own spectrum onto a
/// common [0,1]. Two candidates:
///
///   Rel       u = ln(E/E_cut) / ln(E0/E_cut) -- log-uniform over the node's
///             whole support. Simple, but it does NOT align the spectrum's
///             dominant feature: the Compton backscatter edge sits at
///             E_bs/E0 = 1/(1 + 2 E0/511), which sweeps from 0.81 at 59.5 keV
///             to 0.28 at 662 keV to 0.089 at 2614 keV. Under Rel the edge
///             therefore slides across the bands as E0 changes -- precisely
///             the misalignment that makes an E0 interpolation smear features.
///
///   Landmark  a two-block coordinate that PINS the edge: [E_cut, E_bs] maps
///             to [0, 0.5] and [E_bs, E0] maps to [0.5, 1], log-uniform inside
///             each block. E_bs lands at u = 0.5 and E0 at u = 1 for EVERY
///             node, so the same u means the same physical place in the
///             spectrum. This is the coordinate the landmark energy grid
///             already builds, which is a second and independent reason to
///             adopt it.
///
/// (511 keV is deliberately NOT a third block here: the annihilation continuum
/// only exists above the 1022 keV pair threshold, so a 511 block would change
/// the coordinate's STRUCTURE partway along the E0 axis -- exactly the kind of
/// discontinuity this scaling exists to remove. The backscatter edge, by
/// contrast, exists for every E0.)
enum class EScale { Absolute, Rel, Landmark };

double backscatter_edge(double e0_keV)
{
    return e0_keV / (1.0 + 2.0 * e0_keV / 510.998950);
}

/// Map a photon energy onto [0,1] within its node's own spectrum.
double scaled_u(double e_keV, double e0_keV, double e_cut, EScale mode)
{
    if (mode == EScale::Rel) {
        if (!(e0_keV > e_cut) || !(e_keV > 0.0)) return 0.0;
        return std::log(std::max(e_keV, e_cut) / e_cut) /
               std::log(e0_keV / e_cut);
    }
    // Landmark: [e_cut, E_bs] -> [0, 0.5], [E_bs, E0] -> [0.5, 1].
    const double ebs = backscatter_edge(e0_keV);
    if (!(ebs > e_cut) || !(e0_keV > ebs)) return 0.0;
    if (e_keV <= ebs)
        return 0.5 * std::log(std::max(e_keV, e_cut) / e_cut) /
               std::log(ebs / e_cut);
    return 0.5 + 0.5 * std::log(e_keV / ebs) / std::log(e0_keV / ebs);
}

/// Band index for a scaled coordinate: kNBands equal bins on [0,1]. With
/// kNBands = 6 there is an edge exactly at u = 0.5, so under Landmark the
/// backscatter edge is itself a band boundary.
int band_of_u(double u)
{
    int b = static_cast<int>(u * kNBands);
    if (b < 0) b = 0;
    if (b >= kNBands) b = kNBands - 1;
    return b;
}

struct NodeResult {
    double coord = 0.0;
    double s_hat = 0.0;
    double s_hat_sigma = 0.0;
    std::vector<double> band;        ///< [kNBands] s_hat, mu-summed
    std::vector<double> band_sigma;  ///< BATCH-MEANS sigma (see below)
    /// The same band sigma computed the naive way (per-bin sigmas added in
    /// quadrature), kept only to MEASURE how badly quadrature understates.
    std::vector<double> band_sigma_quad;
    /// Fixed 12 registered-energy bands x 8 detector-frame mu bands for the
    /// total continuum. Exact atoms remain separate, as they do in the archive.
    std::vector<double> joint;
    std::vector<double> joint_sigma;
    int rounds = 0;
    uint64_t histories = 0;
    bool total_precision_met = false;
    bool absolute_precision_met = false;
    bool relative_band_precision_met = false;
    double max_significant_absolute_se_fraction = 0.0;

    /// Relative sigma of the worst band carrying at least `min_frac` of the
    /// flux. Bands below that are ignored: driving a 0.1%-of-flux band to a
    /// tight sigma costs enormously and changes no delivered spectrum.
    double worst_band_rel_sigma(double min_frac) const
    {
        double worst = 0.0;
        for (std::size_t b = 0; b < band.size(); ++b) {
            if (!(s_hat > 0.0) || band[b] < min_frac * s_hat) continue;
            if (band[b] > 0.0)
                worst = std::max(worst, band_sigma[b] / band[b]);
        }
        return worst;
    }
};

struct JointMetric {
    double coord = 0.0;
    double integrated_relative_error = 0.0;
    double joint_tv = 0.0;
    double energy_marginal_tv = 0.0;
    double angular_marginal_tv = 0.0;
    double fluence_coverage_fraction = 0.0;
};

/// Bake one node, running batches until the relative sigma of the total meets
/// `precision` (or `max_histories` is exhausted).
///
/// BAND UNCERTAINTIES COME FROM BATCH MEANS, NOT QUADRATURE. Summing per-bin
/// sigmas in quadrature across the bins of a band assumes the bins are
/// independent within a history. They are not: one primary scores at many
/// vertices into many bins at once, and a history whose walk happens to pass
/// near the detector scores large in ALL of them, so the covariance is
/// POSITIVE and quadrature UNDER-states the band's uncertainty -- the unsafe
/// direction, by about 1.5x for the scattered total in this regime. Each
/// round here is an
/// independent set of histories with its own seed, so the spread of the
/// round-level band values is an exact uncertainty estimate that needs no
/// assumption about correlation at all. Both are reported so the size of the
/// quadrature error is measured rather than assumed.
NodeResult run_node(double e0, double d, double h_src, double h_det,
                    double scoring_radius, double precision, uint64_t batch,
                    uint64_t max_histories, int threads, uint64_t salt,
                    double band_precision = 0.0, double band_min_frac = 0.01,
                    double absolute_precision = 0.001,
                    double absolute_min_frac = 0.01,
                    int native_energy_bins = 384, int native_mu_bins = 64,
                    const std::string& geometry_mode = "field",
                    const std::string& material_name = "soil",
                    EScale escale = EScale::Absolute, double e_cut = 30.0)
{
    const Material air = make_Air();
    const Material slab = material_name == "concrete"
        ? make_Concrete() : make_Soil();
    const Eigen::Vector3d src(0.0, 0.0, -d);
    EnvironmentGeometry env = geometry_mode == "air"
        ? EnvironmentGeometry(air)
        : EnvironmentGeometry::field(air, slab, src, h_src, h_det);
    if (geometry_mode == "field" && env.num_slabs() == 0) {
        // field() silently degenerates to open air when |h_src - h_det| >= d
        // -- a DIFFERENT geometry family; interpolating across that boundary
        // is meaningless (a ladder whose end node crosses it shows a spurious
        // -21% "error" at the last h midpoint). Fail loudly instead.
        std::fprintf(stderr,
                     "ERROR: degenerate field geometry (d=%g h_src=%g "
                     "h_det=%g): |h_src-h_det| >= d leaves no slab\n",
                     d, h_src, h_det);
        std::exit(1);
    }

    EnvScatterConfig cfg;
    cfg.scoring_radius_cm = scoring_radius;
    // Lowering the floor is provably additive (Compton strictly loses energy,
    // Rayleigh preserves it, so nothing re-enters from below) -- the tally
    // floor and the walk floor must move together or flux is transported and
    // then silently dropped.
    cfg.axes.e_min_keV = e_cut;
    cfg.e_floor_keV = e_cut;
    // REGISTERED per-node axis spanning THIS node's [e_cut, E0], with the
    // backscatter edge always at the block boundary halfway through. Bin k
    // then means the same
    // fractional place in the spectrum at every E0, so a band is just a range
    // of bin INDICES and no coordinate transform is needed at comparison time.
    // A fixed absolute grid would instead give a 59 keV source only ~15
    // channels against ~47 for a 2614 keV one, which reads as a spurious
    // support mismatch across E0.
    if (escale != EScale::Absolute)
        cfg.axes = make_registered_landmark_axes(
            e0, native_energy_bins, native_mu_bins, e_cut);
    const NextEventEstimator est(env, cfg, src);

    const FluenceAxes& ax_p = cfg.axes;
    // Per-band accumulation has to be recomputed each round to test the
    // band-precision stopping rule, so hoist the reducer.
    auto reduce_bands = [&](const EnvScatterResult& res, double direct,
                            std::vector<double>& band,
                            std::vector<double>& bsig) {
        band.assign(kNBands, 0.0);
        bsig.assign(kNBands, 0.0);
        if (!(direct > 0.0)) return;
        for (int ie = 0; ie < ax_p.n_energy; ++ie) {
            int b;
            if (escale == EScale::Absolute) {
                b = band_of(ax_p.energy_low(ie));
            } else if (escale == EScale::Landmark) {
                // Registered axis: the band IS a bin-index range.
                b = std::min(kNBands - 1, ie * kNBands / ax_p.n_energy);
            } else {
                const double em = std::sqrt(ax_p.energy_low(ie) *
                                            ax_p.energy_high(ie));
                b = band_of_u(scaled_u(em, e0, e_cut, escale));
            }
            for (int im = 0; im < ax_p.n_mu; ++im) {
                band[b] += res.tally.mean(ie, im) / direct;
                const double s = res.tally.sigma(ie, im) / direct;
                bsig[b] += s * s;
            }
        }
        for (int b = 0; b < kNBands; ++b) bsig[b] = std::sqrt(bsig[b]);
    };

    // Reduce the total continuum onto a fixed registered 12x8 diagnostic grid.
    // Empty cells are deliberately retained as zero entries. Exact atoms stay
    // separate and are not smeared into continuum cells.
    auto reduce_joint = [&](const EnvScatterResult& res, double direct,
                            std::vector<double>& joint) {
        joint.assign(kJointCells, 0.0);
        if (!(direct > 0.0)) return;
        const FluenceAxes& ax = res.tally.axes();
        for (int ie = 0; ie < ax.n_energy; ++ie) {
            const int eb = joint_energy_band(ie, ax.n_energy);
            for (int im = 0; im < ax.n_mu; ++im) {
                const int mb = joint_mu_band(im, ax.n_mu);
                joint[joint_index(eb, mb)] +=
                    res.tally.mean(ie, im) / direct;
            }
        }
    };

    // Minimum rounds before the batch-means variance is meaningful (a spread
    // estimated from a handful of batches is itself very noisy).
    constexpr int kMinRounds = 12;

    EnvScatterResult total(cfg.axes, e0);
    std::vector<std::vector<double>> round_band;   // [round][band]
    std::vector<std::vector<double>> round_joint;  // [round][12*8]
    std::vector<double> round_total;               // [round]
    uint64_t done = 0;
    uint64_t round = 0;
    double direct = 0.0;

    // Batch-means sigma of a per-round series: sigma of the MEAN over rounds.
    auto batch_sigma = [](const std::vector<double>& v) {
        const int n = static_cast<int>(v.size());
        if (n < 2) return 0.0;
        double m = 0.0;
        for (double x : v) m += x;
        m /= n;
        double s2 = 0.0;
        for (double x : v) s2 += (x - m) * (x - m);
        return std::sqrt(s2 / (n - 1) / n);
    };
    auto series_mean = [](const std::vector<double>& v) {
        double m = 0.0;
        for (double x : v) m += x;
        return v.empty() ? 0.0 : m / v.size();
    };

    while (true) {
        const uint64_t per = std::max<uint64_t>(1, batch / threads);
        std::vector<EnvScatterResult> parts(
            threads, EnvScatterResult(cfg.axes, e0));
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; ++t) {
            pool.emplace_back([&, t]() {
                est.run(e0, per,
                        0x9E3779B97F4A7C15ULL * (t + 1) + salt * 7919u +
                            round * 104729u,
                        parts[t]);
            });
        }
        for (auto& th : pool) th.join();

        // This ROUND on its own -- an independent replicate.
        EnvScatterResult rd(cfg.axes, e0);
        for (const auto& p : parts) {
            rd.tally.merge(p.tally);
            rd.direct_transmission = p.direct_transmission;
        }
        direct = rd.direct_transmission / (4.0 * M_PI * d * d);
        std::vector<double> bb, bs;
        reduce_bands(rd, direct, bb, bs);
        round_band.push_back(bb);
        std::vector<double> jj;
        reduce_joint(rd, direct, jj);
        round_joint.push_back(std::move(jj));
        round_total.push_back(direct > 0.0 ? rd.tally.total_mean() / direct : 0.0);

        total.tally.merge(rd.tally);
        total.direct_transmission = rd.direct_transmission;
        done += per * threads;
        ++round;

        if (static_cast<int>(round) >= kMinRounds) {
            const double m = series_mean(round_total);
            const double s = batch_sigma(round_total);
            const bool total_ok = (m > 0.0 && s / m <= precision);
            bool bands_ok = true;
            if (band_precision > 0.0) {
                for (int b = 0; b < kNBands; ++b) {
                    std::vector<double> series;
                    for (const auto& rb : round_band) series.push_back(rb[b]);
                    const double bm = series_mean(series);
                    if (!(m > 0.0) || bm < band_min_frac * m) continue;
                    if (bm > 0.0 && batch_sigma(series) / bm > band_precision)
                        bands_ok = false;
                }
            }
            bool absolute_ok = true;
            if (absolute_precision > 0.0 && m > 0.0) {
                auto check_series = [&](const std::vector<double>& series) {
                    const double value = series_mean(series);
                    if (value >= absolute_min_frac * m &&
                        batch_sigma(series) > absolute_precision * m)
                        absolute_ok = false;
                };
                for (int cell = 0; cell < kJointCells; ++cell) {
                    std::vector<double> series;
                    for (const auto& rj : round_joint)
                        series.push_back(rj[cell]);
                    check_series(series);
                }
                for (int eb = 0; eb < kJointEnergyBands; ++eb) {
                    std::vector<double> series(round_joint.size(), 0.0);
                    for (std::size_t ir = 0; ir < round_joint.size(); ++ir)
                        for (int mb = 0; mb < kJointMuBands; ++mb)
                            series[ir] += round_joint[ir][joint_index(eb, mb)];
                    check_series(series);
                }
                for (int mb = 0; mb < kJointMuBands; ++mb) {
                    std::vector<double> series(round_joint.size(), 0.0);
                    for (std::size_t ir = 0; ir < round_joint.size(); ++ir)
                        for (int eb = 0; eb < kJointEnergyBands; ++eb)
                            series[ir] += round_joint[ir][joint_index(eb, mb)];
                    check_series(series);
                }
            }
            if (total_ok && bands_ok && absolute_ok) break;
        }
        if (done >= max_histories) break;
    }
    total.direct_fluence = direct;

    NodeResult r;
    r.histories = done;
    r.rounds = static_cast<int>(round);
    // Values from the full merged tally (best estimate); uncertainties from the
    // spread of the independent rounds (correlation-proof).
    r.s_hat = total.tally.total_mean() / total.direct_fluence;
    r.s_hat_sigma = batch_sigma(round_total);
    const FluenceAxes& ax = cfg.axes;
    r.band.assign(kNBands, 0.0);
    r.band_sigma.assign(kNBands, 0.0);
    r.band_sigma_quad.assign(kNBands, 0.0);
    for (int ie = 0; ie < ax.n_energy; ++ie) {
        // MUST match reduce_bands' assignment exactly. Any other banding here
        // would silently compare two different quantities: the batch-means
        // sigma of one banding against the values and quadrature sigma of
        // another.
        int b;
        if (escale == EScale::Absolute) {
            b = band_of(ax.energy_low(ie));
        } else if (escale == EScale::Landmark) {
            b = std::min(kNBands - 1, ie * kNBands / ax.n_energy);
        } else {
            const double em = std::sqrt(ax.energy_low(ie) * ax.energy_high(ie));
            b = band_of_u(scaled_u(em, e0, e_cut, escale));
        }
        for (int im = 0; im < ax.n_mu; ++im) {
            r.band[b] += total.tally.mean(ie, im) / total.direct_fluence;
            const double s = total.tally.sigma(ie, im) / total.direct_fluence;
            r.band_sigma_quad[b] += s * s;
        }
    }
    for (int b = 0; b < kNBands; ++b) {
        r.band_sigma_quad[b] = std::sqrt(r.band_sigma_quad[b]);
        std::vector<double> series;
        for (const auto& rb : round_band) series.push_back(rb[b]);
        r.band_sigma[b] = batch_sigma(series);
    }
    reduce_joint(total, total.direct_fluence, r.joint);
    r.joint_sigma.assign(kJointCells, 0.0);
    for (int cell = 0; cell < kJointCells; ++cell) {
        std::vector<double> series;
        for (const auto& rj : round_joint) series.push_back(rj[cell]);
        r.joint_sigma[cell] = batch_sigma(series);
    }

    const double total_mean = series_mean(round_total);
    r.total_precision_met = r.rounds >= kMinRounds && total_mean > 0.0 &&
        batch_sigma(round_total) <= precision * total_mean;
    r.relative_band_precision_met = r.rounds >= kMinRounds;
    if (band_precision > 0.0 && total_mean > 0.0) {
        for (int b = 0; b < kNBands; ++b) {
            std::vector<double> series;
            for (const auto& rb : round_band) series.push_back(rb[b]);
            const double value = series_mean(series);
            if (value >= band_min_frac * total_mean && value > 0.0 &&
                batch_sigma(series) > band_precision * value)
                r.relative_band_precision_met = false;
        }
    }
    r.absolute_precision_met = r.rounds >= kMinRounds;
    r.max_significant_absolute_se_fraction = 0.0;
    if (absolute_precision > 0.0 && total_mean > 0.0) {
        auto record_series = [&](const std::vector<double>& series) {
            const double value = series_mean(series);
            if (value < absolute_min_frac * total_mean) return;
            const double fraction = batch_sigma(series) / total_mean;
            r.max_significant_absolute_se_fraction = std::max(
                r.max_significant_absolute_se_fraction, fraction);
            if (fraction > absolute_precision) r.absolute_precision_met = false;
        };
        for (int cell = 0; cell < kJointCells; ++cell) {
            std::vector<double> series;
            for (const auto& rj : round_joint) series.push_back(rj[cell]);
            record_series(series);
        }
        for (int eb = 0; eb < kJointEnergyBands; ++eb) {
            std::vector<double> series(round_joint.size(), 0.0);
            for (std::size_t ir = 0; ir < round_joint.size(); ++ir)
                for (int mb = 0; mb < kJointMuBands; ++mb)
                    series[ir] += round_joint[ir][joint_index(eb, mb)];
            record_series(series);
        }
        for (int mb = 0; mb < kJointMuBands; ++mb) {
            std::vector<double> series(round_joint.size(), 0.0);
            for (std::size_t ir = 0; ir < round_joint.size(); ++ir)
                for (int eb = 0; eb < kJointEnergyBands; ++eb)
                    series[ir] += round_joint[ir][joint_index(eb, mb)];
            record_series(series);
        }
    }
    return r;
}

// --- interpolation scheme matrix -------------------------------------------

enum class XT { Lin, Log };
enum class YT { Lin, Log };
enum class MT { Linear, Pchip };

struct Scheme {
    const char* name;
    XT xt;
    YT yt;
    MT mt;
};

const Scheme kSchemes[] = {
    {"lin(s) vs lin(x)",   XT::Lin, YT::Lin, MT::Linear},
    {"lin(s) vs log(x)",   XT::Log, YT::Lin, MT::Linear},
    {"log(s) vs lin(x)",   XT::Lin, YT::Log, MT::Linear},
    {"log(s) vs log(x)",   XT::Log, YT::Log, MT::Linear},
    {"PCHIP s  vs lin(x)", XT::Lin, YT::Lin, MT::Pchip},
    {"PCHIP s  vs log(x)", XT::Log, YT::Lin, MT::Pchip},
    {"PCHIP ln s vs lin(x)", XT::Lin, YT::Log, MT::Pchip},
    {"PCHIP ln s vs log(x)", XT::Log, YT::Log, MT::Pchip},
};
constexpr int kNSchemes = 8;

/// Interpolate y(x) from the node arrays under one scheme. `skip` is the node
/// index pair bracketing the query (used only by the linear variants).
double interp_scheme(const Scheme& s, const std::vector<double>& xs,
                     const std::vector<double>& ys, double x)
{
    std::vector<double> tx(xs.size()), ty(ys.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
        tx[i] = (s.xt == XT::Log) ? std::log(xs[i]) : xs[i];
        ty[i] = (s.yt == YT::Log) ? std::log(std::max(ys[i], 1e-300)) : ys[i];
    }
    const double qx = (s.xt == XT::Log) ? std::log(x) : x;

    double v = 0.0;
    if (s.mt == MT::Pchip) {
        Pchip p(tx, ty);
        v = p(qx);
    } else {
        // piecewise linear in the transformed variables
        std::size_t i = 0;
        while (i + 2 < tx.size() && tx[i + 1] < qx) ++i;
        const double f = (qx - tx[i]) / (tx[i + 1] - tx[i]);
        v = ty[i] + f * (ty[i + 1] - ty[i]);
    }
    return (s.yt == YT::Log) ? std::exp(v) : v;
}

/// 1-sigma uncertainty of an INTERPOLATED value, propagated from the grid
/// nodes' own MC noise.
///
/// This is the correction that makes a shape verdict trustworthy. Comparing
/// (interpolated - truth) against the TRUTH's sigma alone silently assumes the
/// interpolant is exact, but it is built from noisy nodes; for a band holding
/// a small share of the flux the node noise is the LARGER term, so that naive z
/// overstates significance. Propagated by resampling because PCHIP is not
/// linear in the node
/// values (its slope estimates use a weighted harmonic mean plus monotonicity
/// clamps), so an analytic Jacobian would be an approximation of its own.
///
/// Nodes are INDEPENDENT MC runs (distinct seed salts), so no cross terms.
double interp_sigma(const Scheme& s, const std::vector<double>& xs,
                    const std::vector<double>& ys,
                    const std::vector<double>& ysig, double x,
                    std::mt19937_64& rng, int draws = 512)
{
    std::vector<double> pert(ys.size());
    double sum = 0.0, sum2 = 0.0;
    for (int k = 0; k < draws; ++k) {
        for (std::size_t i = 0; i < ys.size(); ++i) {
            std::normal_distribution<double> g(ys[i], ysig[i]);
            // Keep positive: the log transforms are undefined otherwise, and a
            // fluence cannot be negative.
            pert[i] = std::max(g(rng), 1e-300);
        }
        const double v = interp_scheme(s, xs, pert, x);
        sum += v;
        sum2 += v * v;
    }
    const double m = sum / draws;
    return std::sqrt(std::max(0.0, sum2 / draws - m * m));
}

} // namespace

int main(int argc, char** argv)
{
    std::string axis = "d";
    std::string json_report, config_sha256;
    std::string geometry_mode = "field", material_name = "soil";
    double e0 = 662.0, scoring_radius = 4.0, precision = 0.0025;
    double band_precision = 0.0, band_min_frac = 0.01;
    double absolute_precision = 0.001, absolute_min_frac = 0.01;
    bool null_test = false;
    EScale escale = EScale::Absolute;
    double e_cut = 30.0;
    uint64_t batch = 1000000, max_histories = 40000000;
    int native_energy_bins = 384, native_mu_bins = 64;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() { return (i + 1 < argc) ? argv[++i] : ""; };
        if (a == "--axis") axis = next();
        else if (a == "--mode") geometry_mode = next();
        else if (a == "--material") material_name = next();
        else if (a == "--e0") e0 = std::atof(next());
        else if (a == "--precision") precision = std::atof(next());
        else if (a == "--band-precision") band_precision = std::atof(next());
        else if (a == "--band-min-frac") band_min_frac = std::atof(next());
        else if (a == "--absolute-precision")
            absolute_precision = std::atof(next());
        else if (a == "--absolute-min-frac")
            absolute_min_frac = std::atof(next());
        else if (a == "--json-report") json_report = next();
        else if (a == "--config-sha256") config_sha256 = next();
        else if (a == "--native-energy-bins")
            native_energy_bins = std::atoi(next());
        else if (a == "--native-mu-bins") native_mu_bins = std::atoi(next());
        else if (a == "--null") null_test = true;
        else if (a == "--escale") {
            const std::string v = next();
            escale = (v == "landmark") ? EScale::Landmark
                   : (v == "rel")      ? EScale::Rel
                                       : EScale::Absolute;
        }
        else if (a == "--e-cut") e_cut = std::atof(next());
        else if (a == "--batch") batch = std::strtoull(next(), nullptr, 10);
        else if (a == "--max-histories") max_histories = std::strtoull(next(), nullptr, 10);
        else if (a == "--threads") threads = std::atoi(next());
        else if (a == "--scoring-radius") scoring_radius = std::atof(next());
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 1; }
    }
    if (threads <= 0 || batch == 0 || max_histories == 0 ||
        !(precision > 0.0) || !(absolute_precision > 0.0) ||
        absolute_min_frac < 0.0 || absolute_min_frac > 1.0 ||
        band_precision < 0.0 || band_min_frac < 0.0 || band_min_frac > 1.0) {
        std::fprintf(stderr, "invalid precision, history, or thread option\n");
        return 1;
    }
    if (native_energy_bins <= 0 || native_mu_bins <= 0 ||
        native_energy_bins % kJointEnergyBands != 0 ||
        native_mu_bins % kJointMuBands != 0) {
        std::fprintf(stderr, "native bins must be positive multiples of the "
                             "12x8 publication grid\n");
        return 1;
    }
    if ((geometry_mode != "air" && geometry_mode != "field") ||
        (material_name != "soil" && material_name != "concrete")) {
        std::fprintf(stderr, "--mode must be air or field and --material must "
                             "be soil or concrete\n");
        return 1;
    }
    if (!json_report.empty() && !is_lower_sha256(config_sha256)) {
        std::fprintf(stderr, "--json-report requires --config-sha256 as 64 "
                             "lowercase hexadecimal characters\n");
        return 1;
    }
    if (!json_report.empty() &&
        (axis != "e0" || escale != EScale::Landmark)) {
        std::fprintf(stderr, "--json-report is publication evidence and "
                             "requires --axis e0 --escale landmark\n");
        return 1;
    }

    std::vector<double> grid;
    if (axis == "d") {
        for (int i = 0; i < 8; ++i) grid.push_back(50.0 * std::pow(1.93, i));
    } else if (axis == "h") {
        // Stop before the |h_src - h_det| >= d validity boundary at d = 500.
        for (int i = 0; i < 5; ++i) grid.push_back(30.0 * std::pow(1.93, i));
    } else if (axis == "e0") {
        grid = {59.5, 122.0, 239.0, 478.0, 955.0, 1836.0, 2614.0};
    } else {
        std::fprintf(stderr, "--axis must be d, h, or e0\n");
        return 1;
    }
    std::vector<double> mid;
    for (std::size_t i = 0; i + 1 < grid.size(); ++i)
        mid.push_back(std::sqrt(grid[i] * grid[i + 1]));


    const double d_fix = 500.0, h_fix = 100.0;

    // NULL TEST: bake the SAME point twice with different seeds and compare.
    // The true difference is exactly zero, so every z here is pure measurement
    // error -- if the uncertainties are right, |z| should look like N(0,1) and
    // NOT systematically exceed ~2. This validates the batch-means sigma end to
    // end, which matters because the quantity being validated (a band's
    // uncertainty under intra-history correlation) is precisely what the naive
    // quadrature estimate gets wrong.
    if (null_test) {
        const double c = grid[grid.size() / 2];
        std::printf("NULL TEST at %s = %.4g (same geometry, two seeds)\n",
                    axis.c_str(), c);
        NodeResult a, b;
        if (axis == "d") {
            a = run_node(e0, c, h_fix, h_fix, scoring_radius, precision, batch,
                         max_histories, threads, 101, band_precision, band_min_frac,
                         absolute_precision, absolute_min_frac,
                         native_energy_bins, native_mu_bins,
                         geometry_mode, material_name,
                         escale, e_cut);
            b = run_node(e0, c, h_fix, h_fix, scoring_radius, precision, batch,
                         max_histories, threads, 202, band_precision, band_min_frac,
                         absolute_precision, absolute_min_frac,
                         native_energy_bins, native_mu_bins,
                         geometry_mode, material_name,
                         escale, e_cut);
        } else if (axis == "h") {
            a = run_node(e0, d_fix, c, h_fix, scoring_radius, precision, batch,
                         max_histories, threads, 101, band_precision, band_min_frac,
                         absolute_precision, absolute_min_frac,
                         native_energy_bins, native_mu_bins,
                         geometry_mode, material_name,
                         escale, e_cut);
            b = run_node(e0, d_fix, c, h_fix, scoring_radius, precision, batch,
                         max_histories, threads, 202, band_precision, band_min_frac,
                         absolute_precision, absolute_min_frac,
                         native_energy_bins, native_mu_bins,
                         geometry_mode, material_name,
                         escale, e_cut);
        } else {
            a = run_node(c, d_fix, h_fix, h_fix, scoring_radius, precision, batch,
                         max_histories, threads, 101, band_precision, band_min_frac,
                         absolute_precision, absolute_min_frac,
                         native_energy_bins, native_mu_bins,
                         geometry_mode, material_name,
                         escale, e_cut);
            b = run_node(c, d_fix, h_fix, h_fix, scoring_radius, precision, batch,
                         max_histories, threads, 202, band_precision, band_min_frac,
                         absolute_precision, absolute_min_frac,
                         native_energy_bins, native_mu_bins,
                         geometry_mode, material_name,
                         escale, e_cut);
        }
        const double dz = (a.s_hat - b.s_hat) /
            std::sqrt(a.s_hat_sigma * a.s_hat_sigma +
                      b.s_hat_sigma * b.s_hat_sigma);
        std::printf("  total    : %.5f +/- %.5f  vs  %.5f +/- %.5f   "
                    "diff %+.2f%%  z %+.2f\n",
                    a.s_hat, a.s_hat_sigma, b.s_hat, b.s_hat_sigma,
                    100.0 * (a.s_hat - b.s_hat) / b.s_hat, dz);
        for (int bb = 0; bb < kNBands; ++bb) {
            const double sd = std::sqrt(a.band_sigma[bb] * a.band_sigma[bb] +
                                        b.band_sigma[bb] * b.band_sigma[bb]);
            const double sq = std::sqrt(a.band_sigma_quad[bb] * a.band_sigma_quad[bb] +
                                        b.band_sigma_quad[bb] * b.band_sigma_quad[bb]);
            const double diff = a.band[bb] - b.band[bb];
            std::printf("  %4.0f-%-5.0f: diff %+7.3f%% total   z(batch-means) %+5.2f"
                        "   z(quadrature) %+5.2f   [%.1f%% of flux]\n",
                        kBandEdges[bb], kBandEdges[bb + 1],
                        100.0 * diff / std::max(b.s_hat, 1e-300),
                        (sd > 0.0) ? diff / sd : 0.0,
                        (sq > 0.0) ? diff / sq : 0.0,
                        100.0 * a.band[bb] / a.s_hat);
        }
        std::printf("\nIf the batch-means z values look like N(0,1) while the "
                    "quadrature ones are systematically larger,\nthe quadrature "
                    "estimate is the broken one -- which is the claim being "
                    "tested.\n");
        return 0;
    }
    auto bake = [&](double c, uint64_t salt) {
        if (axis == "d")
            return run_node(e0, c, h_fix, h_fix, scoring_radius, precision,
                            batch, max_histories, threads, salt,
                            band_precision, band_min_frac, absolute_precision,
                            absolute_min_frac, native_energy_bins,
                            native_mu_bins, geometry_mode, material_name,
                            escale, e_cut);
        if (axis == "h")
            return run_node(e0, d_fix, c, h_fix, scoring_radius, precision,
                            batch, max_histories, threads, salt,
                            band_precision, band_min_frac, absolute_precision,
                            absolute_min_frac, native_energy_bins,
                            native_mu_bins, geometry_mode, material_name,
                            escale, e_cut);
        return run_node(c, d_fix, h_fix, h_fix, scoring_radius, precision,
                        batch, max_histories, threads, salt,
                        band_precision, band_min_frac, absolute_precision,
                        absolute_min_frac, native_energy_bins, native_mu_bins,
                        geometry_mode, material_name, escale, e_cut);
    };

    std::printf("held-out gate: mode=%s material=%s axis=%s e0=%.1f "
                "escale=%s e_cut=%.0f keV "
                "precision=%.3f%%; relative-band=%.3f%%; "
                "absolute-SE<=%.3f%% total (joint cells/marginals >=%.1f%%) "
                "R=%.0f cm\n",
                geometry_mode.c_str(),
                geometry_mode == "air" ? "air" : material_name.c_str(),
                axis.c_str(), e0,
                escale == EScale::Absolute ? "absolute" :
                escale == EScale::Rel ? "rel" : "landmark",
                e_cut, 100.0 * precision, 100.0 * band_precision,
                100.0 * absolute_precision, 100.0 * absolute_min_frac,
                scoring_radius);
    std::printf("native axis: %d registered energy x %d mu bins; "
                "publication reduction: %dx%d\n",
                native_energy_bins, native_mu_bins,
                kJointEnergyBands, kJointMuBands);
    std::printf("grid:");
    for (double g : grid) std::printf(" %.4g", g);
    std::printf("\n\n");

    std::vector<NodeResult> gr, mr;
    uint64_t salt = 1;
    for (double g : grid) {
        gr.push_back(bake(g, salt++));
        gr.back().coord = g;
        std::printf("  node %8.4g: s_hat = %.5f +/- %.5f (%.2f%%)  "
                    "[%.1fM; max absSE %.2f%%T; %s]\n",
                    g, gr.back().s_hat, gr.back().s_hat_sigma,
                    100.0 * gr.back().s_hat_sigma / gr.back().s_hat,
                    gr.back().histories / 1e6,
                    100.0 * gr.back().max_significant_absolute_se_fraction,
                    gr.back().total_precision_met &&
                            gr.back().absolute_precision_met &&
                            gr.back().relative_band_precision_met
                        ? "targets met" : "capped before targets");
        std::fflush(stdout);
    }
    for (double m : mid) {
        mr.push_back(bake(m, salt++));
        mr.back().coord = m;
        std::printf("  MID  %8.4g: s_hat = %.5f +/- %.5f (%.2f%%)  "
                    "[%.1fM; max absSE %.2f%%T; %s]\n",
                    m, mr.back().s_hat, mr.back().s_hat_sigma,
                    100.0 * mr.back().s_hat_sigma / mr.back().s_hat,
                    mr.back().histories / 1e6,
                    100.0 * mr.back().max_significant_absolute_se_fraction,
                    mr.back().total_precision_met &&
                            mr.back().absolute_precision_met &&
                            mr.back().relative_band_precision_met
                        ? "targets met" : "capped before targets");
        std::fflush(stdout);
    }

    std::vector<double> xs, ys;
    for (const auto& n : gr) { xs.push_back(n.coord); ys.push_back(n.s_hat); }

    // --- scheme matrix on the total ---
    std::printf("\nINTERPOLATION SCHEME MATRIX (total s_hat, held-out midpoints)\n");
    std::printf("%-22s %8s %8s   per-midpoint %%err\n", "scheme", "worst", "rms");
    int best = -1;
    double best_worst = 1e300;
    for (int k = 0; k < kNSchemes; ++k) {
        double worst = 0.0, sum2 = 0.0;
        std::printf("%-22s", kSchemes[k].name);
        std::vector<double> errs;
        for (std::size_t i = 0; i < mr.size(); ++i) {
            const double t = mr[i].s_hat;
            const double p = interp_scheme(kSchemes[k], xs, ys, mr[i].coord);
            const double e = 100.0 * (p - t) / t;
            errs.push_back(e);
            worst = std::max(worst, std::abs(e));
            sum2 += e * e;
        }
        const double rms = std::sqrt(sum2 / errs.size());
        std::printf(" %7.2f%% %7.2f%%  ", worst, rms);
        for (double e : errs) std::printf("%+6.2f ", e);
        std::printf("\n");
        if (worst < best_worst) { best_worst = worst; best = k; }
    }
    std::printf("\nbest scheme: %s (worst %.2f%%)\n", kSchemes[best].name,
                best_worst);

    // --- z-scores for the best scheme (is the miss resolved?) ---
    std::printf("\n%10s %12s %10s %9s %8s\n", "midpoint", "truth", "sigma",
                "err(best)", "z");
    for (std::size_t i = 0; i < mr.size(); ++i) {
        const double t = mr[i].s_hat, s = mr[i].s_hat_sigma;
        const double p = interp_scheme(kSchemes[best], xs, ys, mr[i].coord);
        std::printf("%10.4g %12.5f %10.5f %8.2f%% %+8.1f\n", mr[i].coord, t, s,
                    100.0 * (p - t) / t, (p - t) / s);
    }

    // Every one of the predefined 12x8 cells participates, including cells
    // whose observed mass is exactly zero at one or more training nodes. Shape
    // interpolation is in linear mass space; log(0) must never become a hidden
    // coverage filter. The x transform and linear/PCHIP family still follow the
    // total's selected scheme.
    Scheme shape_scheme = kSchemes[best];
    shape_scheme.yt = YT::Lin;
    std::vector<JointMetric> joint_metrics;
    double max_integrated = 0.0, max_joint = 0.0;
    double max_energy = 0.0, max_angular = 0.0;
    std::printf("\nREGISTERED JOINT CONTINUUM SHAPE (12 energy-index x 8 mu bands)\n");
    std::printf("%10s %11s %10s %10s %10s %10s\n", "midpoint",
                "|total err|", "joint TV", "energy TV", "mu TV", "cont/total");
    for (std::size_t i = 0; i < mr.size(); ++i) {
        std::vector<double> predicted(kJointCells, 0.0);
        for (int cell = 0; cell < kJointCells; ++cell) {
            std::vector<double> cell_values;
            cell_values.reserve(gr.size());
            for (const auto& node : gr)
                cell_values.push_back(std::max(0.0, node.joint[cell]));
            // PCHIP and linear interpolation are shape preserving between
            // nonnegative nodes. Clamp only roundoff-sized negatives so the
            // normalized metric remains a probability distance.
            predicted[cell] = std::max(
                0.0, interp_scheme(shape_scheme, xs, cell_values, mr[i].coord));
        }

        double predicted_total = 0.0, truth_total = 0.0;
        for (int cell = 0; cell < kJointCells; ++cell) {
            predicted_total += predicted[cell];
            truth_total += mr[i].joint[cell];
        }
        JointMetric metric;
        metric.coord = mr[i].coord;
        const double predicted_integral = interp_scheme(
            kSchemes[best], xs, ys, mr[i].coord);
        metric.integrated_relative_error = std::abs(
            predicted_integral - mr[i].s_hat) /
            std::max(mr[i].s_hat, 1e-300);
        metric.fluence_coverage_fraction = truth_total /
            std::max(mr[i].s_hat, 1e-300);

        std::vector<double> predicted_energy(kJointEnergyBands, 0.0);
        std::vector<double> truth_energy(kJointEnergyBands, 0.0);
        std::vector<double> predicted_mu(kJointMuBands, 0.0);
        std::vector<double> truth_mu(kJointMuBands, 0.0);
        if (predicted_total > 0.0 && truth_total > 0.0) {
            for (int eb = 0; eb < kJointEnergyBands; ++eb) {
                for (int mb = 0; mb < kJointMuBands; ++mb) {
                    const int cell = joint_index(eb, mb);
                    metric.joint_tv += std::abs(
                        predicted[cell] / predicted_total -
                        mr[i].joint[cell] / truth_total);
                    predicted_energy[eb] += predicted[cell];
                    truth_energy[eb] += mr[i].joint[cell];
                    predicted_mu[mb] += predicted[cell];
                    truth_mu[mb] += mr[i].joint[cell];
                }
            }
            metric.joint_tv *= 0.5;
            for (int eb = 0; eb < kJointEnergyBands; ++eb)
                metric.energy_marginal_tv += std::abs(
                    predicted_energy[eb] / predicted_total -
                    truth_energy[eb] / truth_total);
            for (int mb = 0; mb < kJointMuBands; ++mb)
                metric.angular_marginal_tv += std::abs(
                    predicted_mu[mb] / predicted_total -
                    truth_mu[mb] / truth_total);
            metric.energy_marginal_tv *= 0.5;
            metric.angular_marginal_tv *= 0.5;
        } else {
            metric.joint_tv = metric.energy_marginal_tv =
                metric.angular_marginal_tv = 1.0;
        }
        joint_metrics.push_back(metric);
        max_integrated = std::max(
            max_integrated, metric.integrated_relative_error);
        max_joint = std::max(max_joint, metric.joint_tv);
        max_energy = std::max(max_energy, metric.energy_marginal_tv);
        max_angular = std::max(max_angular, metric.angular_marginal_tv);
        std::printf("%10.4g %10.3f%% %9.3f%% %9.3f%% %9.3f%% %9.3f%%\n",
                    metric.coord, 100.0 * metric.integrated_relative_error,
                    100.0 * metric.joint_tv,
                    100.0 * metric.energy_marginal_tv,
                    100.0 * metric.angular_marginal_tv,
                    100.0 * metric.fluence_coverage_fraction);
    }
    std::printf("maximum    %10.3f%% %9.3f%% %9.3f%% %9.3f%%\n",
                100.0 * max_integrated, 100.0 * max_joint,
                100.0 * max_energy, 100.0 * max_angular);

    // --- SHAPE ---------------------------------------------------------
    // Two ways a consumer can reconstruct a band, compared head to head:
    //   DIRECT   interpolate the band's own value across nodes.
    //   FRACTION interpolate the band's SHARE of the total, and multiply by
    //            the separately interpolated total. This is the "interpolate
    //            the shape, not the magnitude" idea: if the spectrum's SHAPE
    //            varies more smoothly with geometry than its magnitude does,
    //            the fraction is the better-behaved quantity and the accurate
    //            total (measured above at <=1%) carries the scale.
    // Both are reported with FULL error propagation (node noise + truth noise),
    // so a large % error in a sparse band is not mistaken for a real defect.
    // How badly does adding per-bin sigmas in quadrature understate a band's
    // uncertainty? Quadrature assumes bins are independent within a history;
    // they are positively correlated, so it reads low. Measuring the factor
    // here keeps the shape verdict honest AND documents the size of the
    // understatement (about 1.5x for the scattered total).
    std::printf("\nBAND SIGMA: batch-means / quadrature (>1 means quadrature "
                "UNDER-states)\n%-14s", "band keV");
    for (const auto& n : gr) std::printf("%8.4g ", n.coord);
    std::printf("\n");
    for (int b = 0; b < kNBands; ++b) {
        bool any = false;
        for (const auto& n : gr) if (n.band[b] > 0.0) any = true;
        if (!any) continue;
        if (escale == EScale::Absolute)
            std::printf("%6.0f-%-7.0f", kBandEdges[b], kBandEdges[b + 1]);
        else
            std::printf("u %.2f-%-8.2f", double(b) / kNBands,
                        double(b + 1) / kNBands);
        for (const auto& n : gr) {
            const double q = n.band_sigma_quad[b];
            std::printf("%8.2f ", (q > 0.0) ? n.band_sigma[b] / q : 0.0);
        }
        std::printf("\n");
    }

    std::mt19937_64 prng(20260805);
    // Spectral bands may have a zero Monte-Carlo estimate at one node even
    // though their registered support exists.  Log-y interpolation cannot
    // represent that observation, and requiring every node to be positive
    // would cut the E0 comparison to a small fraction of the spectrum.  Shape
    // is therefore interpolated in linear mass/fraction space while retaining
    // the selected x transform and linear/PCHIP family.
    std::printf("\nSHAPE, scheme = %s, z includes NODE noise + truth noise\n",
                kSchemes[best].name);
    std::printf("  (band flux share in brackets; |z| < ~2 means the miss is "
                "not resolved)\n\n");

    std::vector<double> tot_y, tot_sig;
    for (const auto& n : gr) { tot_y.push_back(n.s_hat); tot_sig.push_back(n.s_hat_sigma); }

    for (int b = 0; b < kNBands; ++b) {
        std::vector<double> by, bsig, fy, fsig;
        bool ok = true;
        for (const auto& n : gr) {
            if (!(n.s_hat > 0.0)) { ok = false; break; }
            by.push_back(std::max(0.0, n.band[b]));
            bsig.push_back(n.band_sigma[b]);
            const double f = std::max(0.0, n.band[b]) / n.s_hat;
            fy.push_back(f);
            // sigma of a RATIO whose numerator is part of its denominator:
            // with T = b + r and r independent of b,
            //   Var(b/T) = (r^2 Var(b) + b^2 Var(r)) / T^4,  Var(r)=Var(T)-Var(b).
            // Treating b and T as independent would overstate this and unfairly
            // penalize the fraction method.
            const double vb = n.band_sigma[b] * n.band_sigma[b];
            const double vT = n.s_hat_sigma * n.s_hat_sigma;
            const double r = n.s_hat - n.band[b];
            const double vr = std::max(0.0, vT - vb);
            const double T2 = n.s_hat * n.s_hat;
            fsig.push_back(std::sqrt(std::max(0.0,
                (r * r * vb + n.band[b] * n.band[b] * vr) / (T2 * T2))));
        }
        if (!ok) continue;

        const double share = 100.0 * mr[0].band[b] / std::max(mr[0].s_hat, 1e-300);
        if (escale == EScale::Absolute)
            std::printf("  %4.0f-%-5.0f keV [%.1f%% of flux]\n",
                        kBandEdges[b], kBandEdges[b + 1], share);
        else
            std::printf("  u %.2f-%.2f [%.1f%% of flux]\n",
                        double(b) / kNBands, double(b + 1) / kNBands, share);

        for (int method = 0; method < 2; ++method) {
            std::printf("    %-8s", method == 0 ? "direct" : "fraction");
            for (std::size_t i = 0; i < mr.size(); ++i) {
                const double t = mr[i].band[b];
                const double ts = mr[i].band_sigma[b];
                double p = 0.0, ps = 0.0;
                if (method == 0) {
                    p = interp_scheme(shape_scheme, xs, by, mr[i].coord);
                    ps = interp_sigma(shape_scheme, xs, by, bsig,
                                      mr[i].coord, prng);
                } else {
                    const double f = interp_scheme(shape_scheme, xs, fy,
                                                   mr[i].coord);
                    const double fs = interp_sigma(shape_scheme, xs, fy, fsig,
                                                   mr[i].coord, prng);
                    const double T = interp_scheme(kSchemes[best], xs, tot_y,
                                                   mr[i].coord);
                    const double Ts = interp_sigma(kSchemes[best], xs, tot_y,
                                                   tot_sig, mr[i].coord, prng);
                    p = f * T;
                    ps = std::abs(p) * std::sqrt(
                        (f > 0.0 ? (fs / f) * (fs / f) : 0.0) +
                        (T > 0.0 ? (Ts / T) * (Ts / T) : 0.0));
                }
                const double denom = std::sqrt(ps * ps + ts * ts);
                const double z = (denom > 0.0) ? (p - t) / denom : 0.0;
                std::printf(" %+7.3f%%T(z%+5.1f)",
                            100.0 * (p - t) /
                                std::max(mr[i].s_hat, 1e-300), z);
            }
            std::printf("\n");
        }
    }

    // --- single-number shape distortion: fraction of the spectrum misplaced --
    // 0.5 * sum_b |p_b/P - t_b/T| is the share of the spectrum in the wrong
    // band after removing any overall normalization error -- the same metric
    // hpge_binning_study uses, so the two are directly comparable. The floor is
    // what MC noise alone contributes (E|noise| = sigma sqrt(2/pi)).
    std::printf("\nSPECTRUM DISTORTION (%% of spectrum in the wrong band, "
                "normalization removed)\n");
    std::printf("%-12s", "midpoint");
    for (const auto& m : mr) std::printf("%9.4g ", m.coord);
    std::printf("\n");
    for (int method = 0; method < 2; ++method) {
        std::printf("%-12s", method == 0 ? "direct" : "fraction");
        for (std::size_t i = 0; i < mr.size(); ++i) {
            std::vector<double> pb(kNBands, 0.0);
            std::vector<bool> used(kNBands, false);
            double P = 0.0, T = 0.0;
            for (int b = 0; b < kNBands; ++b) {
                std::vector<double> by, fy;
                bool ok = true;
                for (const auto& n : gr) {
                    if (!(n.s_hat > 0.0)) { ok = false; break; }
                    by.push_back(std::max(0.0, n.band[b]));
                    fy.push_back(std::max(0.0, n.band[b]) / n.s_hat);
                }
                // Registered bands have support at every E0.  A zero observed
                // mass is a valid noisy estimate, not a reason to remove that
                // band from the shape normalization.
                if (!ok) continue;
                pb[b] = (method == 0)
                    ? interp_scheme(shape_scheme, xs, by, mr[i].coord)
                    : interp_scheme(shape_scheme, xs, fy, mr[i].coord) *
                      interp_scheme(kSchemes[best], xs, tot_y, mr[i].coord);
                used[b] = true;
                P += pb[b];
                T += mr[i].band[b];
            }
            double l1 = 0.0;
            for (int b = 0; b < kNBands; ++b) {
                if (!used[b] || P <= 0.0 || T <= 0.0) continue;
                l1 += std::abs(pb[b] / P - mr[i].band[b] / T);
            }
            // Flag when the comparison covers only part of the spectrum, which
            // is the normal case on the E0 axis (nodes have different upper
            // limits) and makes the number a partial-spectrum statement.
            double all_T = 0.0;
            for (int b = 0; b < kNBands; ++b) all_T += mr[i].band[b];
            const double cover = (all_T > 0.0) ? T / all_T : 0.0;
            if (cover > 0.999) std::printf("%8.2f%%  ", 50.0 * l1);
            else std::printf("%6.2f%%[%.0f%%] ", 50.0 * l1, 100.0 * cover);
        }
        std::printf("\n");
    }
    std::printf("%-12s", "noise floor");
    for (std::size_t i = 0; i < mr.size(); ++i) {
        double T = 0.0, floor = 0.0;
        for (int b = 0; b < kNBands; ++b) T += mr[i].band[b];
        for (int b = 0; b < kNBands; ++b)
            floor += std::sqrt(2.0 / M_PI) * mr[i].band_sigma[b];
        std::printf("%8.2f%%  ", (T > 0.0) ? 50.0 * floor / T : 0.0);
    }
    std::printf("\n");

    if (!json_report.empty()) {
        bool all_stopping = true;
        for (const auto& node : gr)
            all_stopping = all_stopping && node.total_precision_met &&
                node.absolute_precision_met && node.relative_band_precision_met;
        for (const auto& node : mr)
            all_stopping = all_stopping && node.total_precision_met &&
                node.absolute_precision_met && node.relative_band_precision_met;
        const bool canonical_noise_targets = precision <= 0.0025 &&
            absolute_precision <= 0.001 && absolute_min_frac <= 0.01;
        const bool canonical_axes = native_energy_bins == 384 &&
            native_mu_bins == 64 && std::abs(e_cut - 10.0) <= 1e-12;
        const bool passed = axis == "e0" && escale == EScale::Landmark &&
            canonical_axes && canonical_noise_targets && all_stopping &&
            max_integrated <= 0.01 && max_joint <= 0.01 &&
            max_energy <= 0.005 && max_angular <= 0.005;

        std::ostringstream out;
        out << std::setprecision(17)
            << "{\n"
            << "  \"schema\": \"CEEFIT-REGISTERED-E0-GATE-v1\",\n"
            << "  \"config_sha256\": \"" << config_sha256 << "\",\n"
            << "  \"passed\": " << (passed ? "true" : "false") << ",\n"
            << "  \"predefined_index_bands\": true,\n"
            << "  \"statistically_empty_bands_included\": true,\n"
            << "  \"shape_coverage_fraction\": 1.0,\n"
            << "  \"max_integrated_relative_error\": " << max_integrated << ",\n"
            << "  \"max_joint_tv\": " << max_joint << ",\n"
            << "  \"max_energy_marginal_tv\": " << max_energy << ",\n"
            << "  \"max_angular_marginal_tv\": " << max_angular << ",\n"
            << "  \"all_node_stopping_targets_met\": "
            << (all_stopping ? "true" : "false") << ",\n"
            << "  \"best_scheme\": \"" << kSchemes[best].name << "\",\n"
            << "  \"observable\": \"total_continuum; exact lines remain atoms\",\n"
            << "  \"gate_configuration\": {\n"
            << "    \"geometry_tier\": \"" << geometry_mode << "\",\n"
            << "    \"material\": \""
            << (geometry_mode == "air" ? "Air" :
                material_name == "concrete" ? "Concrete" : "Soil")
            << "\",\n"
            << "    \"distance_cm\": " << d_fix << ",\n"
            << "    \"source_height_cm\": " << h_fix << ",\n"
            << "    \"detector_height_cm\": " << h_fix << ",\n"
            << "    \"scoring_radius_cm\": " << scoring_radius << ",\n"
            << "    \"energy_cut_keV\": " << e_cut << "\n"
            << "  },\n"
            << "  \"grid\": {\n"
            << "    \"native_registered_energy_bins\": "
            << native_energy_bins << ",\n"
            << "    \"native_mu_bins\": " << native_mu_bins << ",\n"
            << "    \"registered_energy_bands\": " << kJointEnergyBands << ",\n"
            << "    \"registered_energy_edges\": [";
        for (int edge = 0; edge <= kJointEnergyBands; ++edge) {
            if (edge) out << ", ";
            out << double(edge) / kJointEnergyBands;
        }
        out << "],\n"
            << "    \"mu_bands\": " << kJointMuBands << ",\n"
            << "    \"mu_edges\": [";
        for (int edge = 0; edge <= kJointMuBands; ++edge) {
            if (edge) out << ", ";
            out << -1.0 + 2.0 * edge / kJointMuBands;
        }
        out << "]\n"
            << "  },\n"
            << "  \"stopping_targets\": {\n"
            << "    \"total_relative_se\": " << precision << ",\n"
            << "    \"absolute_se_fraction_of_total\": "
            << absolute_precision << ",\n"
            << "    \"significant_fraction_of_total\": "
            << absolute_min_frac << ",\n"
            << "    \"batch_histories\": " << batch << ",\n"
            << "    \"max_histories_per_node\": " << max_histories << "\n"
            << "  },\n"
            << "  \"acceptance_thresholds\": {\n"
            << "    \"max_integrated_relative_error\": 0.01,\n"
            << "    \"max_joint_tv\": 0.01,\n"
            << "    \"max_energy_marginal_tv\": 0.005,\n"
            << "    \"max_angular_marginal_tv\": 0.005\n"
            << "  },\n"
            << "  \"midpoints\": [\n";
        for (std::size_t i = 0; i < joint_metrics.size(); ++i) {
            const auto& metric = joint_metrics[i];
            out << "    {\"energy_keV\": " << metric.coord
                << ", \"integrated_relative_error\": "
                << metric.integrated_relative_error
                << ", \"joint_tv\": " << metric.joint_tv
                << ", \"energy_marginal_tv\": "
                << metric.energy_marginal_tv
                << ", \"angular_marginal_tv\": "
                << metric.angular_marginal_tv
                << ", \"continuum_fraction_of_scattered\": "
                << metric.fluence_coverage_fraction << "}";
            out << (i + 1 == joint_metrics.size() ? "\n" : ",\n");
        }
        out << "  ],\n"
            << "  \"nodes\": [\n";
        bool first_node = true;
        auto write_node = [&](const NodeResult& node, const char* role) {
            if (!first_node) out << ",\n";
            first_node = false;
            out << "    {\"role\": \"" << role << "\", \"energy_keV\": "
                << node.coord << ", \"histories\": " << node.histories
                << ", \"total_continuum_s_hat\": " << node.s_hat
                << ", \"joint_cell_masses\": [";
            for (std::size_t cell = 0; cell < node.joint.size(); ++cell) {
                if (cell) out << ", ";
                out << node.joint[cell];
            }
            out << "]"
                << ", \"total_relative_se\": "
                << node.s_hat_sigma / std::max(node.s_hat, 1e-300)
                << ", \"max_significant_absolute_se_fraction\": "
                << node.max_significant_absolute_se_fraction
                << ", \"total_target_met\": "
                << (node.total_precision_met ? "true" : "false")
                << ", \"absolute_target_met\": "
                << (node.absolute_precision_met ? "true" : "false")
                << ", \"relative_band_target_met\": "
                << (node.relative_band_precision_met ? "true" : "false") << "}";
        };
        for (const auto& node : gr) write_node(node, "grid");
        for (const auto& node : mr) write_node(node, "heldout");
        out << "\n  ]\n}\n";

        const std::string temporary = json_report + ".tmp";
        {
            std::ofstream file(temporary, std::ios::binary | std::ios::trunc);
            if (!file || !(file << out.str()) || !file.flush()) {
                std::fprintf(stderr, "failed to write JSON report %s\n",
                             temporary.c_str());
                return 2;
            }
        }
        if (std::rename(temporary.c_str(), json_report.c_str()) != 0) {
            const int saved_errno = errno;
            std::remove(temporary.c_str());
            std::fprintf(stderr, "failed to publish JSON report %s: %s\n",
                         json_report.c_str(), std::strerror(saved_errno));
            return 2;
        }
        std::printf("JSON report: %s (%s)\n", json_report.c_str(),
                    passed ? "PASS" : "FAIL");
    }
    return 0;
}
