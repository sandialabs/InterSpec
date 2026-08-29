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

/// @file ResponseGenerator.cpp
/// @brief MC generation driver for DetectorResponse (see ResponseGenerator.h).

#include "io/ResponseGenerator.h"

#include "cross_sections/CrossSectionData.h"
#include "efficiency/EfficiencyCalculator.h"
#include "io/EfficiencyTransfer.h"  // kTotalScatterInRecapture

#include <Eigen/Dense>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <thread>

namespace ceelo {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Deterministic per-node seed (never 0); stage separates scan families.
uint64_t node_seed(uint64_t base, uint32_t stage, uint32_t node) {
    uint64_t s = base * 0x9E3779B97F4A7C15ULL;
    s ^= (uint64_t)stage << 40;
    s ^= (uint64_t)node;
    s ^= s >> 30;
    s *= 0xBF58476D1CE4E5B9ULL;
    s ^= s >> 27;
    return s | 1ULL;
}

// Halton low-discrepancy value: index i (>=0), prime base.
double halton(uint64_t i, uint32_t base) {
    double f = 1.0, r = 0.0;
    for (uint64_t n = i + 1; n > 0; n /= base) {
        f /= base;
        r += f * static_cast<double>(n % base);
    }
    return r;
}

// Per-node termination budget + batch sizing from the generation options.
// The CPU cap is the real, load-invariant budget; the wall cap is only a
// stuck-run backstop. The first merged wave is threads * batch_size events,
// so tie the batch size to min_events -- cheap low-E nodes previously
// overshot min_events ~5x through first-wave granularity.
void apply_node_budget(SimulationConfig& cfg, const GenerationOptions& opts,
                       double target_prec) {
    cfg.termination.target_fep_rel_precision = target_prec;
    cfg.termination.max_events = opts.max_events_per_node;
    cfg.termination.max_cpu_seconds = opts.max_cpu_seconds_per_node;
    const unsigned threads = opts.num_threads
        ? opts.num_threads
        : std::max(1u, std::thread::hardware_concurrency());
    if (opts.max_cpu_seconds_per_node > 0.0)
        cfg.termination.max_wall_seconds =
            10.0 * opts.max_cpu_seconds_per_node / threads;
    cfg.termination.min_events = opts.min_events_per_node;
    cfg.num_threads = opts.num_threads;
    cfg.batch_size = std::max<uint64_t>(
        2000, opts.min_events_per_node / std::max(1u, threads));
}

// Per-node target FEP precision: the graded map when set, else the scalar.
// Null map => exactly the historical scalar path (bit-identical goldens).
double resolve_node_precision(const GenerationOptions& opts, uint32_t stage,
                              double E_keV) {
    return opts.node_precision ? opts.node_precision(stage, E_keV)
                               : opts.node_fep_precision;
}

// Shared per-run state: configured calculator + progress/cancel bookkeeping.
struct Runner {
    const GeometryDescriptor& gd;
    const GenerationOptions& opts;
    std::vector<std::unique_ptr<Material>> owned;
    int nodes_done = 0;
    int nodes_total = 1;

    Runner(const GeometryDescriptor& g, const GenerationOptions& o)
        : gd(g), opts(o) {
        owned.reserve(gd.materials.size());
        for (const MaterialSpec& spec : gd.materials)
            owned.push_back(std::make_unique<Material>(spec.to_material()));
    }

    const Material* mat(int idx) const {
        if (idx < 0 || idx >= static_cast<int>(owned.size()))
            throw std::runtime_error("ResponseGenerator: bad material index");
        return owned[static_cast<size_t>(idx)].get();
    }

    void configure(EfficiencyCalculator& calc) const {
        // Same mapping as the public helper, but reusing this Runner's
        // already-instantiated materials (one instantiation per run, not per
        // node).
        calc.set_detector(gd.shape, mat(gd.crystal_material_index),
                          gd.dimensions_cm);
        if (gd.bore) calc.set_bore_hole(gd.bore->radius, gd.bore->depth);
        if (gd.dead_layer)
            calc.set_dead_layer(gd.dead_layer->front, gd.dead_layer->side,
                                gd.dead_layer->back);
        for (const LayerSpec& l : gd.layers)
            calc.add_attenuator(mat(l.material_index), l.front_thickness_cm,
                                l.side_thickness_cm, l.z_start_cm, l.z_end_cm);
        if (gd.collimator)
            calc.add_collimator(mat(gd.collimator->material_index),
                                gd.collimator->side_thickness_cm,
                                gd.collimator->z_start_cm,
                                gd.collimator->z_end_cm);
    }

    void tick(const std::string& stage) {
        ++nodes_done;
        if (opts.cancel && opts.cancel->load())
            throw GenerationCancelled();
        if (opts.progress)
            opts.progress(std::min(1.0, double(nodes_done) / nodes_total), stage);
    }

    // fep_only: FEP-only transport (skips the eps_tot tally) for nodes whose
    // total efficiency is never consumed. Currently unused pending the
    // bore-hole kill-on-exit fix (see the stage-3 call site).
    EfficiencyResult run_node(double energy_keV, const Eigen::Vector3d& src,
                              uint32_t stage, uint32_t node,
                              bool fep_only = false) {
        EfficiencyCalculator calc;
        configure(calc);
        calc.set_point_source(src);
        calc.enable_fep_only_mode(fep_only);
        SimulationConfig cfg;
        cfg.energy_keV = energy_keV;
        apply_node_budget(cfg, opts,
                          resolve_node_precision(opts, stage, energy_keV));
        cfg.seed = node_seed(opts.base_seed, stage, node);
        EfficiencyResult r = calc.compute(cfg);
        if (opts.stats_out) {
            NodeStat ns;
            ns.stage = stage;
            ns.energy_keV = energy_keV;
            const double d = src.norm();
            ns.d_cm = d;
            ns.cos_theta = (d > 0.0) ? (-src.z() / d) : 1.0;
            ns.events = r.num_events_simulated;
            ns.cpu_s = r.cpu_time_seconds;
            ns.wall_s = r.wall_time_seconds;
            ns.stop = static_cast<uint8_t>(r.stop_reason);
            ns.fep_rel_prec = (r.full_energy_peak_efficiency > 0.0)
                ? r.fep_uncertainty / r.full_energy_peak_efficiency
                : 0.0;
            opts.stats_out->add(ns);
        }
        return r;
    }
};

// One MC-scored probe point at (E, d, ct, phi) at the given precision + seed.
// Shared by every probe/validation bank (plain, adversarial, structured). Fills
// stats_out (as a stage-`stage` node) when the caller supplies one.
ProbePoint run_probe_point(Runner& run, const GenerationOptions& opts, double E,
                           double d, double ct, double phi_deg, double precision,
                           uint64_t seed, uint8_t tag, uint32_t stage) {
    const Eigen::Vector3d src = source_position(d, ct, phi_deg * kPi / 180.0);
    ProbePoint row;
    row.energy_keV = E;
    row.d_cm = d;
    row.cos_theta = ct;
    row.phi_deg = phi_deg;
    row.seed = seed;
    row.tag = tag;
    EfficiencyCalculator calc;
    run.configure(calc);
    calc.set_point_source(src);
    SimulationConfig cfg;
    cfg.energy_keV = E;
    apply_node_budget(cfg, opts, precision);
    cfg.seed = seed;
    const EfficiencyResult r = calc.compute(cfg);
    row.eps_fep = r.full_energy_peak_efficiency;
    row.fep_unc = r.fep_uncertainty;
    row.eps_tot = r.total_efficiency;
    row.tot_unc = r.total_uncertainty;
    if (opts.stats_out) {
        NodeStat ns;
        ns.stage = stage;
        ns.energy_keV = E;
        ns.d_cm = d;
        ns.cos_theta = ct;
        ns.events = r.num_events_simulated;
        ns.cpu_s = r.cpu_time_seconds;
        ns.wall_s = r.wall_time_seconds;
        ns.stop = static_cast<uint8_t>(r.stop_reason);
        ns.fep_rel_prec = (r.full_energy_peak_efficiency > 0.0)
            ? r.fep_uncertainty / r.full_energy_peak_efficiency : 0.0;
        opts.stats_out->add(ns);
    }
    return row;
}

// Crystal K-edges within the energy range (the mandatory segment breaks).
std::vector<double> crystal_edges(const GeometryDescriptor& gd, double e_min,
                                  double e_max) {
    std::vector<double> edges;
    if (gd.crystal_material_index < 0) return edges;
    const CrossSectionData& db = CrossSectionData::instance();
    const MaterialSpec& crystal =
        gd.materials[static_cast<size_t>(gd.crystal_material_index)];
    for (const MaterialComponent& c : crystal.composition) {
        if (const FluorescenceData* f = db.fluorescence(c.Z)) {
            const double e = f->k_edge_keV;
            if (e > e_min * 1.02 && e < e_max * 0.98) edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end(),
                            [](double a, double b) { return b - a < 0.01; }),
                edges.end());
    return edges;
}

size_t segment_of(double e, const std::vector<double>& edges) {
    size_t s = 0;
    for (const double x : edges)
        if (e > x) ++s;
    return s;
}

// Greedy 1-D node selection on (x = lnE, curves of ln eta), segment-aware:
// segment endpoints forced, then repeatedly add the worst-interpolated point.
// Noise-aware + tolerance-driven: `max_nodes` is a CAP, not a target - a point
// is only added when its interpolation error genuinely exceeds both `interp_tol`
// and `noise_k` times its local MC sigma. Resolving structure below the MC noise
// floor would just fit scatter, so the greedy stops early there; the cap can be
// set generously and the accuracy target decides how many nodes are actually used.
std::vector<size_t> greedy_energy_nodes(const std::vector<double>& E,
                                        const std::vector<double>& ln_eta,
                                        const std::vector<double>& sigma,
                                        const std::vector<double>& edges,
                                        size_t max_nodes,
                                        double interp_tol, double noise_k) {
    const size_t n = E.size();
    std::set<size_t> sel;
    sel.insert(0);
    sel.insert(n - 1);
    // Segment endpoints (edge flanks) forced.
    for (size_t i = 0; i + 1 < n; ++i)
        if (segment_of(E[i], edges) != segment_of(E[i + 1], edges)) {
            sel.insert(i);
            sel.insert(i + 1);
        }

    auto interp_err = [&](const std::set<size_t>& nodes, size_t q) -> double {
        // Piecewise eval with per-segment PCHIP over the selected nodes.
        const size_t seg_q = segment_of(E[q], edges);
        std::vector<double> xs, ys;
        for (const size_t i : nodes)
            if (segment_of(E[i], edges) == seg_q) {
                xs.push_back(std::log(E[i]));
                ys.push_back(ln_eta[i]);
            }
        if (xs.size() < 2) return 0.0;
        const Pchip f(xs, ys);
        return std::fabs(f(std::log(E[q])) - ln_eta[q]);
    };

    while (sel.size() < std::min(max_nodes, n)) {
        double worst = -1.0;
        size_t worst_i = SIZE_MAX;
        for (size_t q = 0; q < n; ++q) {
            if (sel.count(q)) continue;
            const double e = interp_err(sel, q);
            const double floor = std::max(interp_tol, noise_k * sigma[q]);
            if (e > floor && e > worst) { worst = e; worst_i = q; }
        }
        if (worst_i == SIZE_MAX) break;   // nothing left above tol + MC noise
        sel.insert(worst_i);
    }
    return std::vector<size_t>(sel.begin(), sel.end());
}

// Greedy shared cos-theta node selection across several ln-eta curves.
// Noise-aware + tolerance-driven like greedy_energy_nodes: `max_nodes` is a cap;
// a scan angle is added only when its worst interpolation error (over the shape-
// energy curves) exceeds both `interp_tol` and `noise_k` times that curve's local
// MC sigma. `sigma_curves` parallels `curves`; pass empty to disable the noise gate.
std::vector<size_t> greedy_ct_nodes(const std::vector<double>& cts,
                                    const std::vector<std::vector<double>>& curves,
                                    const std::vector<std::vector<double>>& sigma_curves,
                                    size_t max_nodes,
                                    double interp_tol, double noise_k) {
    const size_t n = cts.size();
    std::set<size_t> sel{0, n - 1, n / 2};
    while (sel.size() < std::min(max_nodes, n)) {
        std::vector<size_t> ns(sel.begin(), sel.end());
        std::vector<double> xs(ns.size());
        for (size_t k = 0; k < ns.size(); ++k) xs[k] = cts[ns[k]];
        double worst = -1.0;
        size_t worst_i = SIZE_MAX;
        for (size_t ci = 0; ci < curves.size(); ++ci) {
            const std::vector<double>& y = curves[ci];
            std::vector<double> yn(ns.size());
            for (size_t k = 0; k < ns.size(); ++k) yn[k] = y[ns[k]];
            const Pchip f(xs, yn);
            for (size_t q = 0; q < n; ++q) {
                if (sel.count(q)) continue;
                const double e = std::fabs(f(cts[q]) - y[q]);
                const double sig = sigma_curves.empty() ? 0.0 : sigma_curves[ci][q];
                const double floor = std::max(interp_tol, noise_k * sig);
                if (e > floor && e > worst) { worst = e; worst_i = q; }
            }
        }
        if (worst_i == SIZE_MAX) break;   // nothing left above tol + MC noise
        sel.insert(worst_i);
    }
    return std::vector<size_t>(sel.begin(), sel.end());
}

// Locally-smoothed node value: weighted linear fit of ln eta vs lnE through
// the dense scan within the node's segment (tri-cube window to the nearest
// selected neighbors) -- "fit through slightly denser data rather than
// interpolating raw noisy nodes" (campaign-measured win).
double smoothed_node_value(size_t node_idx, const std::vector<size_t>& nodes,
                           const std::vector<double>& E,
                           const std::vector<double>& ln_eta,
                           const std::vector<double>& frac_sig,
                           const std::vector<double>& edges) {
    const double x0 = std::log(E[node_idx]);
    const size_t seg0 = segment_of(E[node_idx], edges);
    // Window: distance to the nearest selected neighbor on each side (in ln E).
    double w_lo = 1e300, w_hi = 1e300;
    for (const size_t j : nodes) {
        if (j == node_idx) continue;
        const double dx = std::log(E[j]) - x0;
        if (dx < 0.0) w_lo = std::min(w_lo, -dx);
        else w_hi = std::min(w_hi, dx);
    }
    const double w = std::max(1e-6, 0.9 * std::min(w_lo, w_hi));

    double s_w = 0.0, s_wx = 0.0, s_wy = 0.0, s_wxx = 0.0, s_wxy = 0.0;
    int n_used = 0;
    for (size_t i = 0; i < E.size(); ++i) {
        if (segment_of(E[i], edges) != seg0) continue;
        const double dx = std::log(E[i]) - x0;
        const double u = std::fabs(dx) / w;
        if (u >= 1.0) continue;
        const double tri = std::pow(1.0 - u * u * u, 3);
        const double sg = std::max(frac_sig[i], 1e-4);
        const double wt = tri / (sg * sg);
        s_w += wt;
        s_wx += wt * dx;
        s_wy += wt * ln_eta[i];
        s_wxx += wt * dx * dx;
        s_wxy += wt * dx * ln_eta[i];
        ++n_used;
    }
    if (n_used < 3) return ln_eta[node_idx];
    const double det = s_w * s_wxx - s_wx * s_wx;
    if (std::fabs(det) < 1e-30) return s_wy / s_w;
    return (s_wxx * s_wy - s_wx * s_wxy) / det;  // intercept at dx = 0
}

// ln-E linear interpolation between shape-scan energies.
struct ShapeInterp {
    std::vector<double> ln_es;
    void set(const std::vector<double>& energies) {
        ln_es.resize(energies.size());
        for (size_t i = 0; i < energies.size(); ++i)
            ln_es[i] = std::log(energies[i]);
    }
    // returns (i0, i1, weight)
    void bracket(double energy_keV, size_t& i0, size_t& i1, double& w) const {
        const double x = std::log(energy_keV);
        if (ln_es.size() < 2 || x <= ln_es.front()) { i0 = i1 = 0; w = 0; return; }
        if (x >= ln_es.back()) { i0 = i1 = ln_es.size() - 1; w = 0; return; }
        i1 = static_cast<size_t>(
            std::upper_bound(ln_es.begin(), ln_es.end(), x) - ln_es.begin());
        i0 = i1 - 1;
        w = (x - ln_es[i0]) / (ln_es[i1] - ln_es[i0]);
    }
};

// One eps_tot sample at a FEP node: (E, eps_tot, frac_sigma, K_noRayleigh, K).
struct TotPoint { double E, tot, tot_sig, k_nors, k_tot; };

enum class TotDecision { Kernel, BCurve, NeedTable };

// eps_tot tier auto-choice (spec Eq. 4): eps_tot is free at every FEP node, so
// score the compact tiers held-out and ship the most compact whose 95th-pct
// error <= 1%. Sets the tier + b(E) fields for KernelExact/BCurve; returns
// NeedTable when the eta_tot table is required (the caller supplies angular
// data). b(E) is filled whenever a table is NOT chosen, so a caller without
// angular data (transfer-flat) can force BCurve as the best angle-flat total.
TotDecision decide_tot_tier(DetectorResponse& resp,
                            const std::vector<TotPoint>& tot_points,
                            const GenerationOptions& opts) {
    auto pct95 = [](std::vector<double> v) {
        if (v.empty()) return 1e300;
        std::sort(v.begin(), v.end());
        return v[static_cast<size_t>(0.95 * (v.size() - 1))];
    };
    std::vector<double> err_kernel;
    for (const TotPoint& p : tot_points)
        if (p.tot > 0.0 && p.tot_sig < 0.05)
            err_kernel.push_back(std::fabs(p.k_nors / p.tot - 1.0));

    // b(E): 8 log-spaced hat knots fit to ln(tot / K_nors) (linear LSQ).
    std::vector<double> b_knots;
    const int nb = 8;
    for (int i = 0; i < nb; ++i)
        b_knots.push_back(std::log(opts.e_min_keV) +
                          (std::log(opts.e_max_keV) - std::log(opts.e_min_keV)) *
                              double(i) / (nb - 1));
    Eigen::MatrixXd AtA = Eigen::MatrixXd::Zero(nb, nb);
    Eigen::VectorXd Atb = Eigen::VectorXd::Zero(nb);
    for (const TotPoint& p : tot_points) {
        if (p.tot <= 0.0 || p.k_nors <= 0.0 || p.tot_sig > 0.05) continue;
        const std::vector<double> B = hat_basis(std::log(p.E), b_knots);
        const double y = std::log(p.tot / p.k_nors);
        const double wt = 1.0 / std::max(p.tot_sig, 1e-3);
        for (int i = 0; i < nb; ++i) {
            if (B[i] == 0.0) continue;
            Atb(i) += wt * wt * B[i] * y;
            for (int j = 0; j < nb; ++j)
                AtA(i, j) += wt * wt * B[i] * B[j];
        }
    }
    for (int i = 0; i < nb; ++i) AtA(i, i) += 1e-9;  // regularize empty knots
    const Eigen::VectorXd bc = AtA.ldlt().solve(Atb);

    std::vector<double> err_b;
    for (const TotPoint& p : tot_points) {
        if (p.tot <= 0.0 || p.k_nors <= 0.0 || p.tot_sig > 0.05) continue;
        const std::vector<double> B = hat_basis(std::log(p.E), b_knots);
        double lnb = 0.0;
        for (int i = 0; i < nb; ++i) lnb += B[i] * bc(i);
        err_b.push_back(std::fabs(std::exp(lnb) * p.k_nors / p.tot - 1.0));
    }

    if (pct95(err_kernel) <= 0.01) {
        resp.tot_eff.tier = TotEffTier::KernelExact;
        return TotDecision::Kernel;
    }
    resp.tot_eff.b_energies_keV.clear();
    resp.tot_eff.ln_b.clear();
    for (int i = 0; i < nb; ++i) {
        resp.tot_eff.b_energies_keV.push_back(std::exp(b_knots[i]));
        resp.tot_eff.ln_b.push_back(bc(i));
    }
    if (pct95(err_b) <= 0.01) {
        resp.tot_eff.tier = TotEffTier::BCurve;
        return TotDecision::BCurve;
    }
    return TotDecision::NeedTable;
}

// ===========================================================================
// Closed-loop refinement (D-b): incremental table edits + re-certify.
//
// Assembly-reuse approach: rather than re-running the whole stage-1/2/3
// assembly, the loop performs INCREMENTAL, deterministic table edits on the
// already-built response (spec option 2 -- lower risk):
//   * insert an energy row into eta_fep (reusing the failing probe's own
//     node-precision MC, keeping the interpolated angular shape);
//   * insert a cos-theta column into eta_fep (fresh MC across the E-nodes);
//   * insert a distance plane into the near-field table (fresh MC across
//     (E, ct)).
// Each edit re-finalizes only the affected table. Determinism: same options
// -> same probe/refine sequence and seeds; bit-exact with num_threads = 1
// (the standard library caveat -- refinement pass/fail thresholds can flip on
// reproducible-in-distribution multithreaded MC values, exactly as the greedy
// node selector already does).
// ===========================================================================

// Insert a backbone energy node at E_star with on-axis level `ln_eta_onaxis`
// (from the reused probe MC) and fractional sigma `frac_sig`, keeping the
// currently-modeled angular shape. Pure (no MC). Returns false if capped or too
// close to an existing node.
bool insert_eta_energy_node(DetectorResponse& resp, double E_star,
                            double ln_eta_onaxis, double frac_sig,
                            int max_energy_nodes) {
    EtaTable& t = resp.eta_fep;
    if (static_cast<int>(t.energies_keV.size()) >= max_energy_nodes) return false;
    if (E_star <= t.energies_keV.front() || E_star >= t.energies_keV.back())
        return false;
    for (const double e : t.energies_keV)
        if (std::fabs(std::log(E_star / e)) < 0.02) return false;

    const size_t nc = t.cos_thetas.size();
    const size_t np = t.phis_deg.empty() ? 1 : t.phis_deg.size();
    const size_t stride = nc * np;
    // Angular shape at E_star from the CURRENT finalized table (relative to the
    // on-axis node = max cos-theta), added onto the freshly-measured level.
    std::vector<double> new_ln(stride), new_sig(stride);
    bool clamped = false;
    const double onaxis_ref = t.eval_ln(E_star, t.cos_thetas.back(), 0.0, clamped);
    for (size_t c = 0; c < nc; ++c) {
        const double ctc = t.cos_thetas[c];
        for (size_t p = 0; p < np; ++p) {
            const double phip = (np > 1) ? t.phis_deg[p] : 0.0;
            const double sh = t.eval_ln(E_star, ctc, phip, clamped) - onaxis_ref;
            new_ln[c * np + p] = ln_eta_onaxis + sh;
            const double nsig = t.node_frac_sigma(E_star, ctc, phip);
            new_sig[c * np + p] = std::sqrt(frac_sig * frac_sig + nsig * nsig);
        }
    }
    const size_t pos = static_cast<size_t>(
        std::lower_bound(t.energies_keV.begin(), t.energies_keV.end(), E_star) -
        t.energies_keV.begin());
    t.ln_eta.insert(t.ln_eta.begin() + pos * stride, new_ln.begin(), new_ln.end());
    t.frac_sigma.insert(t.frac_sigma.begin() + pos * stride, new_sig.begin(),
                        new_sig.end());
    t.energies_keV.insert(t.energies_keV.begin() + pos, E_star);
    t.finalize();
    return true;
}

// Insert a cos-theta column at ct_star: fresh MC at (each E-node, ct_star) at
// d_far. Returns false if capped / too close to an existing node. Accumulates
// MC CPU into `cpu_used`; node_id advances the deterministic refine seed.
bool insert_eta_ct_node(DetectorResponse& resp, Runner& run,
                        const GenerationOptions& ref_opts, double ct_star,
                        double d_far, uint32_t seed_stage, uint32_t& node_id,
                        double& cpu_used, int max_ct_nodes) {
    EtaTable& t = resp.eta_fep;
    const size_t nc = t.cos_thetas.size();
    if (static_cast<int>(nc) >= max_ct_nodes) return false;
    if (ct_star <= t.cos_thetas.front() || ct_star >= t.cos_thetas.back())
        return false;
    for (const double c : t.cos_thetas)
        if (std::fabs(c - ct_star) < 0.01) return false;

    const size_t nE = t.energies_keV.size();
    const size_t np = t.phis_deg.empty() ? 1 : t.phis_deg.size();
    std::vector<double> col_ln(nE * np), col_sig(nE * np);
    for (size_t e = 0; e < nE; ++e) {
        const double E = t.energies_keV[e];
        for (size_t p = 0; p < np; ++p) {
            const double phip = (np > 1) ? t.phis_deg[p] : 0.0;
            const Eigen::Vector3d src =
                source_position(d_far, ct_star, phip * kPi / 180.0);
            EfficiencyCalculator calc;
            run.configure(calc);
            calc.set_point_source(src);
            SimulationConfig cfg;
            cfg.energy_keV = E;
            apply_node_budget(cfg, ref_opts, ref_opts.node_fep_precision);
            cfg.seed = node_seed(ref_opts.base_seed, seed_stage, node_id++);
            const EfficiencyResult r = calc.compute(cfg);
            cpu_used += r.cpu_time_seconds;
            const ApertureQuadrature q = resp.make_quadrature(src);
            const double K = resp.kernel_K(E, q, MuChoice::Total);
            const double eps = std::max(r.full_energy_peak_efficiency, 1e-300);
            col_ln[e * np + p] = std::log(eps / std::max(K, 1e-300));
            col_sig[e * np + p] = (r.full_energy_peak_efficiency > 0.0)
                ? r.fep_uncertainty / r.full_energy_peak_efficiency : 0.5;
        }
    }
    const size_t pos = static_cast<size_t>(
        std::lower_bound(t.cos_thetas.begin(), t.cos_thetas.end(), ct_star) -
        t.cos_thetas.begin());
    const size_t new_nc = nc + 1;
    std::vector<double> nl(nE * new_nc * np), nsg(nE * new_nc * np);
    for (size_t e = 0; e < nE; ++e)
        for (size_t c = 0; c < new_nc; ++c)
            for (size_t p = 0; p < np; ++p) {
                double lv, sv;
                if (c < pos) {
                    lv = t.ln_eta[(e * nc + c) * np + p];
                    sv = t.frac_sigma[(e * nc + c) * np + p];
                } else if (c == pos) {
                    lv = col_ln[e * np + p];
                    sv = col_sig[e * np + p];
                } else {
                    lv = t.ln_eta[(e * nc + (c - 1)) * np + p];
                    sv = t.frac_sigma[(e * nc + (c - 1)) * np + p];
                }
                nl[(e * new_nc + c) * np + p] = lv;
                nsg[(e * new_nc + c) * np + p] = sv;
            }
    t.cos_thetas.insert(t.cos_thetas.begin() + pos, ct_star);
    t.ln_eta = std::move(nl);
    t.frac_sigma = std::move(nsg);
    t.finalize();
    return true;
}

// Insert a near-field distance plane at d_star: fresh MC at (each E, each ct)
// -> ln N. Returns false if no near table / too close to an existing plane.
bool insert_near_d_plane(DetectorResponse& resp, Runner& run,
                         const GenerationOptions& ref_opts, double d_star,
                         uint32_t seed_stage, uint32_t& node_id,
                         double& cpu_used) {
    NearFieldModel& nf = resp.near_field;
    if (nf.empty()) return false;
    if (d_star <= nf.dists_cm.front() || d_star >= nf.dists_cm.back())
        return false;
    for (const double d : nf.dists_cm)
        if (std::fabs(std::log(d_star / std::max(d, 1e-9))) < 0.03) return false;

    const size_t nE = nf.energies_keV.size();
    const size_t nc = nf.cos_thetas.size();
    const size_t nd = nf.dists_cm.size();
    std::vector<double> plane_ln(nE * nc), plane_sig(nE * nc);
    for (size_t ie = 0; ie < nE; ++ie) {
        const double E = nf.energies_keV[ie];
        for (size_t ic = 0; ic < nc; ++ic) {
            const double ct = nf.cos_thetas[ic];
            const Eigen::Vector3d src = source_position(d_star, ct);
            EfficiencyCalculator calc;
            run.configure(calc);
            calc.set_point_source(src);
            SimulationConfig cfg;
            cfg.energy_keV = E;
            apply_node_budget(cfg, ref_opts, ref_opts.node_fep_precision);
            cfg.seed = node_seed(ref_opts.base_seed, seed_stage, node_id++);
            const EfficiencyResult r = calc.compute(cfg);
            cpu_used += r.cpu_time_seconds;
            const ApertureQuadrature q = resp.make_quadrature(src);
            const double K = resp.kernel_K(E, q, MuChoice::Total);
            const size_t idx = ie * nc + ic;
            if (r.full_energy_peak_efficiency <= 0.0 || K <= 0.0) {
                plane_ln[idx] = 0.0;
                plane_sig[idx] = 0.5;
                continue;
            }
            bool clamped = false;
            const double ln_eta_far = resp.eta_fep.eval_ln(E, ct, 0.0, clamped);
            plane_ln[idx] =
                std::log(r.full_energy_peak_efficiency / K) - ln_eta_far;
            plane_sig[idx] = r.fep_uncertainty / r.full_energy_peak_efficiency;
        }
    }
    const size_t pos = static_cast<size_t>(
        std::lower_bound(nf.dists_cm.begin(), nf.dists_cm.end(), d_star) -
        nf.dists_cm.begin());
    const size_t new_nd = nd + 1;
    std::vector<double> nl(nE * nc * new_nd), nsg(nE * nc * new_nd);
    for (size_t ie = 0; ie < nE; ++ie)
        for (size_t ic = 0; ic < nc; ++ic)
            for (size_t id = 0; id < new_nd; ++id) {
                double lv, sv;
                if (id < pos) {
                    lv = nf.ln_n[(ie * nc + ic) * nd + id];
                    sv = nf.frac_sigma[(ie * nc + ic) * nd + id];
                } else if (id == pos) {
                    lv = plane_ln[ie * nc + ic];
                    sv = plane_sig[ie * nc + ic];
                } else {
                    lv = nf.ln_n[(ie * nc + ic) * nd + (id - 1)];
                    sv = nf.frac_sigma[(ie * nc + ic) * nd + (id - 1)];
                }
                nl[(ie * nc + ic) * new_nd + id] = lv;
                nsg[(ie * nc + ic) * new_nd + id] = sv;
            }
    nf.dists_cm.insert(nf.dists_cm.begin() + pos, d_star);
    nf.ln_n = std::move(nl);
    nf.frac_sigma = std::move(nsg);
    nf.finalize();
    return true;
}

// The closed-loop driver (dispatched from generate() when opts.closed_loop).
std::shared_ptr<DetectorResponse> generate_closed_loop(
    const GeometryDescriptor& gd, const GenerationOptions& opts_in) {
    auto scaled_progress = [&opts_in](double lo, double hi) {
        return [lo, hi, &opts_in](double f, const std::string& s) {
            if (opts_in.progress)
                opts_in.progress(lo + (hi - lo) * std::min(1.0, std::max(0.0, f)),
                                 s);
        };
    };

    // ---- 1. INITIAL GRID via the ordinary (non-closed-loop) generate() -------
    GenerationOptions init = opts_in;
    init.closed_loop = false;
    if (opts_in.initial_grid == InitialGrid::Coarse) {
        // Memo naive 5 shape-E x 6 ct start: thin the SCAN resolution so the
        // initial grid is sparse (higher initial p95), but LEAVE the node caps
        // (n_energy_nodes / n_cos_theta_nodes) generous -- the refinement loop
        // needs headroom to add nodes back. (The near-field d grid is left at
        // the generator default; coarsening it would touch the untouched
        // fixed-grid generate() path.)
        init.n_shape_energies = std::min(init.n_shape_energies, 5);
        init.n_cos_theta_scan = std::min(init.n_cos_theta_scan, 6);
        init.n_energy_scan = std::min(init.n_energy_scan, 16);
    }
    GenerationStats init_stats;
    init.stats_out = &init_stats;
    init.progress = scaled_progress(0.0, 0.60);
    std::shared_ptr<DetectorResponse> resp = ResponseGenerator::generate(gd, init);
    const double init_cpu = init_stats.total_cpu_s;
    const double refine_budget =
        opts_in.refine_cpu_budget_frac * std::max(init_cpu, 1.0);

    // Probe passes: uniform, structured precision 0.008 (EnergyGap planner drops
    // to node precision internally for reuse).
    GenerationOptions probe_opts = opts_in;
    probe_opts.closed_loop = false;
    probe_opts.node_precision = nullptr;
    probe_opts.precision_profile = PrecisionProfile::Uniform;
    probe_opts.node_fep_precision = 0.008;
    probe_opts.progress = nullptr;
    probe_opts.stats_out = nullptr;
    probe_opts.cancel = opts_in.cancel;

    // Refine MC: uniform at (tight) node precision so inserted nodes are solid.
    GenerationOptions ref_opts = opts_in;
    ref_opts.closed_loop = false;
    ref_opts.node_precision = nullptr;
    ref_opts.precision_profile = PrecisionProfile::Uniform;
    ref_opts.node_fep_precision = std::min(opts_in.node_fep_precision, 0.003);
    ref_opts.progress = nullptr;
    ref_opts.stats_out = nullptr;
    ref_opts.cancel = opts_in.cancel;
    Runner ref_run(gd, ref_opts);

    const double a = gd.transverse_half_extent();
    const double d_far = std::max(opts_in.far_distance_a * a, 10.0);
    const Eigen::Vector3d far_pos(0.0, 0.0, -d_far);

    double refine_cpu = 0.0;
    uint32_t refine_node_id = 0;
    int iters_run = 0;
    bool converged = false;
    bool aborted = false;
    const int n_per_family = 8;

    for (int k = 0; k < std::max(0, opts_in.max_refine_iters); ++k) {
        if (opts_in.cancel && opts_in.cancel->load())
            throw GenerationCancelled();
        iters_run = k + 1;
        if (opts_in.progress)
            opts_in.progress(0.60 + 0.28 * double(k) /
                                        std::max(1, opts_in.max_refine_iters),
                             "Refining");

        // (a) structured probe pass -- disjoint per-iteration seed offset.
        const int start = 20000 + k * 1000;
        const std::vector<ProbePoint> probes =
            ResponseGenerator::structured_probe_bank(
                gd, probe_opts, *resp, kAllProbeFamilies, n_per_family, start);

        // (b) noise-aware fail rule.
        struct Fail { ProbeFamily fam; double E, d, ct, rel, eps_fep, fep_unc; };
        std::vector<Fail> fails;
        int n_random = 0, n_random_fail = 0;
        for (const ProbePoint& p : probes) {
            if (p.eps_fep <= 0.0 || p.fep_unc / p.eps_fep > 0.10) continue;
            const Eigen::Vector3d src =
                source_position(p.d_cm, p.cos_theta, p.phi_deg * kPi / 180.0);
            const EffResult m = resp->eps_fep_at(p.energy_keV, src);
            const double rel = std::fabs(m.value / p.eps_fep - 1.0);
            const double mc_sig = p.fep_unc / p.eps_fep;
            const double node_sig = (m.value > 0.0) ? m.sigma / m.value : 0.0;
            const double thr = std::max(
                opts_in.cert_tol,
                2.0 * std::sqrt(mc_sig * mc_sig + node_sig * node_sig));
            const ProbeFamily fam = static_cast<ProbeFamily>(p.tag);
            if (fam == ProbeFamily::Random) {
                ++n_random;
                if (rel > thr) ++n_random_fail;
                continue;
            }
            if (rel > thr) fails.push_back({fam, p.energy_keV, p.d_cm,
                                            p.cos_theta, rel, p.eps_fep,
                                            p.fep_unc});
        }

        // (c) Random failures only implicate the MODEL FORM when the structured
        // grid is CLEAN (else they are grid-coverage gaps the structured
        // families are already refining). So the canary is evaluated only when
        // there are no structured failures.
        if (fails.empty()) {
            // >10% random fail with a clean structured grid => model-form
            // breakage: LOUD fallback to the fixed full grid. Require >=2
            // failures so a single unlucky probe in a small bank can't abort.
            if (n_random_fail >= 2 && n_random_fail * 10 > n_random) {
                std::fprintf(stderr,
                    "[CeeLo closed-loop] MODEL-FORM CANARY: %d/%d random probes "
                    "fail (>10%%) with a clean structured grid; aborting closed "
                    "loop, falling back to the fixed full-grid path "
                    "(certificate.converged=false).\n",
                    n_random_fail, n_random);
                aborted = true;
                break;
            }
            // A few random fails but structured clean => minor model-form:
            // inflate the fep floors, record, and treat the grid as converged.
            if (n_random_fail > 0) {
                resp->floors.fep_far *= 1.25;
                resp->floors.fep_near *= 1.25;
            }
            converged = true;
            break;
        }

        // (d) attribution + refine. Energy inserts are free (reuse probe MC);
        // ct / near inserts cost MC and respect the refine budget.
        bool did_refine = false;
        int energy_inserts = 0;
        const ApertureQuadrature q_far = resp->make_quadrature(far_pos);

        // Count ShapeE fails for the expensive-refinement gate.
        int shapeE_fails = 0;
        double shapeE_worst = 0.0;
        for (const Fail& f : fails)
            if (f.fam == ProbeFamily::ShapeEGap) {
                ++shapeE_fails;
                shapeE_worst = std::max(shapeE_worst, f.rel);
            }
        const bool shapeE_gate =
            (shapeE_fails >= 2) || (shapeE_worst > 2.0 * opts_in.cert_tol);

        // EnergyGap fails -> insert backbone nodes REUSING the probe's own
        // node-precision MC (the probe was on-axis @ d_far, so its eps_fep IS
        // the backbone level -- zero extra MC).
        for (const Fail& f : fails) {
            if (f.fam != ProbeFamily::EnergyGap) continue;
            if (energy_inserts >= 4) break;  // per-iteration cap
            if (f.eps_fep <= 0.0) continue;
            const double K = resp->kernel_K(f.E, q_far, MuChoice::Total);
            if (K <= 0.0) continue;
            const double ln_eta = std::log(f.eps_fep / K);
            const double fsig = f.fep_unc / f.eps_fep;
            if (insert_eta_energy_node(*resp, f.E, ln_eta, fsig,
                                       opts_in.n_energy_nodes)) {
                did_refine = true;
                ++energy_inserts;
            }
        }

        // ShapeEGap (far) fails, gated -> one on-axis energy insert at the gap.
        if (shapeE_gate && refine_cpu < refine_budget) {
            for (const Fail& f : fails) {
                if (f.fam != ProbeFamily::ShapeEGap) continue;
                if (std::fabs(f.d - d_far) > 1e-6) continue;  // far only
                if (energy_inserts >= 4) break;
                const double K = resp->kernel_K(f.E, q_far, MuChoice::Total);
                if (K <= 0.0) continue;
                EfficiencyCalculator calc;
                ref_run.configure(calc);
                calc.set_point_source(source_position(d_far, 1.0));
                SimulationConfig cfg;
                cfg.energy_keV = f.E;
                apply_node_budget(cfg, ref_opts, ref_opts.node_fep_precision);
                cfg.seed = node_seed(ref_opts.base_seed, /*stage=*/32,
                                     refine_node_id++);
                const EfficiencyResult r = calc.compute(cfg);
                refine_cpu += r.cpu_time_seconds;
                if (r.full_energy_peak_efficiency <= 0.0) continue;
                const double ln_eta =
                    std::log(r.full_energy_peak_efficiency / K);
                const double fsig =
                    r.fep_uncertainty / r.full_energy_peak_efficiency;
                if (insert_eta_energy_node(*resp, f.E, ln_eta, fsig,
                                           opts_in.n_energy_nodes)) {
                    did_refine = true;
                    ++energy_inserts;
                }
            }
        }

        // CtGap fails -> one ct insert (widest failing gap) per iteration.
        if (refine_cpu < refine_budget) {
            double worst_ct = -1.0, worst_ct_rel = -1.0;
            for (const Fail& f : fails)
                if (f.fam == ProbeFamily::CtGap && f.rel > worst_ct_rel) {
                    worst_ct_rel = f.rel;
                    worst_ct = f.ct;
                }
            if (worst_ct > 0.0 &&
                insert_eta_ct_node(*resp, ref_run, ref_opts, worst_ct, d_far,
                                   /*seed_stage=*/33, refine_node_id, refine_cpu,
                                   opts_in.n_cos_theta_nodes))
                did_refine = true;
        }

        // NearCell fails -> one near d-plane (widest failing d) per iteration.
        if (refine_cpu < refine_budget) {
            double worst_d = -1.0, worst_d_rel = -1.0;
            for (const Fail& f : fails)
                if (f.fam == ProbeFamily::NearCell && f.rel > worst_d_rel) {
                    worst_d_rel = f.rel;
                    worst_d = f.d;
                }
            if (worst_d > 0.0 &&
                insert_near_d_plane(*resp, ref_run, ref_opts, worst_d,
                                    /*seed_stage=*/34, refine_node_id,
                                    refine_cpu))
                did_refine = true;
        }

        if (!did_refine) break;              // nothing actionable
        if (refine_cpu >= refine_budget) break;  // budget exhausted
    }

    // ---- 2. Model-form canary fallback: fixed full-grid path -----------------
    if (aborted) {
        GenerationOptions full = opts_in;
        full.closed_loop = false;
        full.initial_grid = InitialGrid::Full;
        full.stats_out = nullptr;
        full.progress = scaled_progress(0.60, 0.90);
        resp = ResponseGenerator::generate(gd, full);
    }

    // ---- 3. CERTIFICATE PASS (full structured family + random @0.005) --------
    if (opts_in.progress) opts_in.progress(0.90, "Certifying");
    GenerationOptions cert_opts = opts_in;
    cert_opts.closed_loop = false;
    cert_opts.progress = nullptr;
    cert_opts.stats_out = nullptr;
    ResponseGenerator::certify(*resp, gd, cert_opts, opts_in.n_cert_probes,
                               /*seed_offset=*/7000, kAllStructuredFamilies);
    resp->certificate.iterations = iters_run;
    // Convergence is HONEST: the refinement loop must have found nothing more to
    // fix (per the noise-aware rule) AND the independent certificate must
    // actually meet the accuracy gate. When the node precision is relaxed enough
    // that the noise-aware tolerance floor sits above cert_tol, the loop stops
    // refining but the certificate p95 stays above cert_tol -- that is a
    // precision-limited response, honestly reported as NOT converged (the
    // stored fep_p95 tells the user to tighten precision, not add nodes).
    resp->certificate.converged =
        converged && !aborted &&
        (resp->certificate.fep_p95 <= opts_in.cert_tol);
    resp->certificate.cpu_seconds += refine_cpu;
    if (opts_in.progress) opts_in.progress(1.0, "Done");
    return resp;
}

}  // namespace

// ---------------------------------------------------------------------------

std::function<double(uint32_t, double)> relax_mild_precision_map(double base) {
    // D0 policy memo `relax_mild`: backbone + low-E rows stay tight; angular
    // (stage 2) and near (stage 3) rows relax by energy toward 3 MeV, where the
    // cap-limited nodes dominate CPU. Cost falls ~(base/p)^2 on precision-
    // limited nodes; hpge's two lowest-E near rows are <=150 keV so they stay
    // at `base` automatically.
    return [base](uint32_t stage, double E_keV) -> double {
        if (stage <= 1) return base;          // backbone: cheap + load-bearing
        if (E_keV <= 150.0) return base;      // low-E rows kept tight
        if (E_keV <= 500.0) return 0.0045;
        if (E_keV <= 1500.0) return 0.006;
        return 0.008;
    };
}

// ---------------------------------------------------------------------------

int ResponseGenerator::estimated_node_count(const GeometryDescriptor& gd,
                                            const GenerationOptions& o) {
    // D-b: nodes_total is not known upfront (refinement is data-driven), so the
    // estimate is initial-grid + an expected-refine + certificate allowance.
    if (o.closed_loop) {
        GenerationOptions init = o;
        init.closed_loop = false;
        if (o.initial_grid == InitialGrid::Coarse) {
            init.n_shape_energies = std::min(init.n_shape_energies, 5);
            init.n_cos_theta_scan = std::min(init.n_cos_theta_scan, 6);
            init.n_energy_scan = std::min(init.n_energy_scan, 16);
        }
        const int initial = estimated_node_count(gd, init);
        // Rough per-iteration refine: ~one ct column (n_energy_nodes MC) + one
        // near d-plane (n_shape x 9) + a few energy inserts; plus the final
        // certificate probes (random + ~4 structured families x 8).
        const int per_iter = o.n_energy_nodes + o.n_shape_energies * 9 + 4;
        const int cert = o.n_cert_probes + 4 * 8;
        return initial + std::max(0, o.max_refine_iters) * per_iter + cert;
    }
    const bool box = gd.shape == DetectorShape::Box;
    const int n_phi = box ? o.n_phi_nodes : 1;
    const int backbone =
        o.n_energy_scan + 2 * static_cast<int>(
            crystal_edges(gd, o.e_min_keV, o.e_max_keV).size());
    // EFFTRAN transfer mode: only the on-axis backbone (+ a few forced angle
    // anchors when n_anchor_angles > 1); no near-field tensor.
    if (o.transfer_mode) {
        const int anchors = (o.n_anchor_angles > 1)
            ? o.n_shape_energies * std::max(3, o.n_anchor_angles) : 0;
        return backbone + anchors;
    }
    const int angular = o.n_shape_energies * o.n_cos_theta_scan * n_phi;
    int near = 0;
    if (o.profile != ResponseProfile::FarField) {
        // 9 cos_theta nodes x (8 or 9) MC distance nodes per shape energy
        // (the outermost distance node is a no-MC ln N = 0 anchor).
        const int nd = (o.profile == ResponseProfile::Contact) ? 9 : 8;
        near = o.n_shape_energies * nd * 9;
    }
    return backbone + angular + near;
}

std::shared_ptr<DetectorResponse> ResponseGenerator::generate(
    const GeometryDescriptor& gd, const GenerationOptions& opts_in) {
    // D-b: opt-in closed-loop refinement. FALSE (default) falls through to the
    // fixed-grid path below VERBATIM (goldens bit-identical); the loop only runs
    // when explicitly enabled.
    if (opts_in.closed_loop)
        return generate_closed_loop(gd, opts_in);

    GenerationOptions opts = opts_in;
    // Graded precision: an explicit node_precision wins; otherwise a RelaxMild
    // profile installs the D0 map. Uniform + null leaves node_precision null =>
    // the exact historical scalar path (goldens bit-identical).
    if (!opts.node_precision &&
        opts.precision_profile == PrecisionProfile::RelaxMild)
        opts.node_precision = relax_mild_precision_map(opts.node_fep_precision);
    // EFFTRAN transfer mode: a far-field product (no near-field MC tensor). The
    // few-anchor variant runs a small forced cos-theta scan; the flat variant
    // (n_anchor_angles <= 1) skips the angular scan entirely (handled below).
    if (opts.transfer_mode) {
        opts.profile = ResponseProfile::FarField;
        if (opts.n_anchor_angles > 1) {
            opts.n_cos_theta_scan = std::max(3, opts.n_anchor_angles);
            opts.n_cos_theta_nodes = opts.n_anchor_angles;
        }
    }
    if (opts.stats_out)
        *opts.stats_out = GenerationStats{};

    Runner run(gd, opts);
    run.nodes_total = estimated_node_count(gd, opts);

    auto resp = std::make_shared<DetectorResponse>();
    resp->descriptor = gd;
    for (size_t i = 0; i < gd.materials.size(); ++i)
        resp->mu_tables.push_back(
            MuTable::sample(*run.mat(static_cast<int>(i)), static_cast<int>(i)));
    resp->provenance.method = opts.transfer_mode ? ProductionMethod::QuickMcTransfer
                                                 : ProductionMethod::FullMc;
    resp->provenance.profile = opts.profile;
    resp->provenance.node_fep_precision = opts.node_fep_precision;
    resp->provenance.generation_seed = opts.base_seed;
    resp->provenance.kernel_n_rays = opts.kernel_n_rays;
    resp->provenance.detector_name = opts.detector_name;
    resp->provenance.valid_e_min_keV = opts.e_min_keV;
    resp->provenance.valid_e_max_keV = opts.e_max_keV;
    resp->finalize();  // geometry + mu lookups; tables still empty

    const double a = gd.transverse_half_extent();
    const double d_far = std::max(opts.far_distance_a * a, 10.0);
    const bool box = gd.shape == DetectorShape::Box;
    const double phi_span =
        box ? ((gd.dimensions_cm.size() > 1 &&
                std::fabs(gd.dimensions_cm[0] - gd.dimensions_cm[1]) < 1e-9)
                   ? 45.0 : 90.0)
            : 0.0;

    const std::vector<double> edges =
        crystal_edges(gd, opts.e_min_keV, opts.e_max_keV);
    resp->eta_fep.edges_keV = edges;

    // ---- 1. energy backbone (on-axis, far field) ---------------------------
    std::vector<double> scan_E;
    for (int i = 0; i < opts.n_energy_scan; ++i)
        scan_E.push_back(opts.e_min_keV *
                         std::pow(opts.e_max_keV / opts.e_min_keV,
                                  double(i) / (opts.n_energy_scan - 1)));
    for (const double e : edges) {
        scan_E.push_back(e * (1.0 - 1e-3));
        scan_E.push_back(e * (1.0 + 1e-3));
    }
    std::sort(scan_E.begin(), scan_E.end());
    scan_E.erase(std::unique(scan_E.begin(), scan_E.end(),
                             [](double x, double y) { return y - x < x * 1e-6; }),
                 scan_E.end());

    const Eigen::Vector3d far_pos(0.0, 0.0, -d_far);
    const ApertureQuadrature q_far = resp->make_quadrature(far_pos);

    std::vector<double> bb_ln_eta(scan_E.size()), bb_sig(scan_E.size());
    std::vector<double> bb_tot(scan_E.size()), bb_tot_sig(scan_E.size());
    for (size_t i = 0; i < scan_E.size(); ++i) {
        const EfficiencyResult r =
            run.run_node(scan_E[i], far_pos, /*stage=*/1, static_cast<uint32_t>(i));
        const double K = resp->kernel_K(scan_E[i], q_far, MuChoice::Total);
        if (r.full_energy_peak_efficiency <= 0.0 || K <= 0.0)
            throw std::runtime_error("ResponseGenerator: zero FEP at backbone node");
        bb_ln_eta[i] = std::log(r.full_energy_peak_efficiency / K);
        bb_sig[i] = r.fep_uncertainty / r.full_energy_peak_efficiency;
        bb_tot[i] = r.total_efficiency;
        bb_tot_sig[i] = r.total_uncertainty / std::max(r.total_efficiency, 1e-300);
        run.tick("Energy backbone");
    }

    const std::vector<size_t> e_nodes = greedy_energy_nodes(
        scan_E, bb_ln_eta, bb_sig, edges, static_cast<size_t>(opts.n_energy_nodes),
        opts.node_interp_tol, opts.node_noise_k);

    // Backbone eps_tot points -- needed by the tier auto-choice in BOTH the
    // transfer-flat path (which returns before the angular scan) and the full
    // path (whose angular loop appends more below).
    std::vector<TotPoint> tot_points;
    for (size_t i = 0; i < scan_E.size(); ++i) {
        const double K_nors =
            resp->kernel_K(scan_E[i], q_far, MuChoice::NoRayleigh);
        const double K_tot = resp->kernel_K(scan_E[i], q_far, MuChoice::Total);
        tot_points.push_back({scan_E[i], bb_tot[i], bb_tot_sig[i], K_nors, K_tot});
    }

    // ---- EFFTRAN transfer, angle-flat variant ------------------------------
    // MC only the on-axis backbone; two identical sentinel cos-theta nodes make
    // eta angle-flat, and the query-time kernel K carries all angular/distance
    // transfer. Skips the angular + near-field MC entirely.
    if (opts.transfer_mode && opts.n_anchor_angles <= 1) {
        EtaTable& t = resp->eta_fep;
        for (const size_t n : e_nodes) t.energies_keV.push_back(scan_E[n]);
        t.cos_thetas = {opts.cos_theta_min, 1.0};
        t.ln_eta.assign(t.energies_keV.size() * 2, 0.0);
        t.frac_sigma.assign(t.ln_eta.size(), 0.0);
        for (size_t e = 0; e < t.energies_keV.size(); ++e) {
            const double bb = smoothed_node_value(e_nodes[e], e_nodes, scan_E,
                                                  bb_ln_eta, bb_sig, edges);
            const double bs = bb_sig[e_nodes[e]];
            for (size_t c = 0; c < 2; ++c) {
                t.ln_eta[t.index(e, c, 0)] = bb;
                t.frac_sigma[t.index(e, c, 0)] = bs;
            }
        }
        if (decide_tot_tier(*resp, tot_points, opts) == TotDecision::NeedTable)
            resp->tot_eff.tier = TotEffTier::BCurve;  // no angular data -> b(E)
        resp->provenance.min_distance_cm = 2.0 * a;  // far-field validity floor
        resp->model_transfer = SigmaTransferModel{};  // honest off-axis/near sigma
        resp->scatter_in_recapture = kTotalScatterInRecapture;  // total near-field
        resp->finalize();
        if (opts.progress) opts.progress(1.0, "Done");
        return resp;
    }

    // ---- 2. angular shape scans --------------------------------------------
    // Shape energies: log-spaced across the range, snapped to selected nodes.
    std::vector<double> shape_E;
    {
        std::set<double> uniq;
        for (int i = 0; i < opts.n_shape_energies; ++i) {
            const double target = opts.e_min_keV *
                std::pow(opts.e_max_keV / opts.e_min_keV,
                         double(i) / std::max(1, opts.n_shape_energies - 1));
            double best = scan_E[e_nodes[0]];
            for (const size_t n : e_nodes)
                if (std::fabs(std::log(scan_E[n] / target)) <
                    std::fabs(std::log(best / target)))
                    best = scan_E[n];
            uniq.insert(best);
        }
        shape_E.assign(uniq.begin(), uniq.end());
    }

    std::vector<double> scan_ct;
    for (int i = 0; i < opts.n_cos_theta_scan; ++i)
        scan_ct.push_back(opts.cos_theta_min +
                          (1.0 - opts.cos_theta_min) * double(i) /
                              (opts.n_cos_theta_scan - 1));

    std::vector<double> phi_nodes;
    if (box)
        for (int i = 0; i < opts.n_phi_nodes; ++i)
            phi_nodes.push_back(phi_span * double(i) /
                                std::max(1, opts.n_phi_nodes - 1));
    const size_t np = box ? phi_nodes.size() : 1;

    // ang_ln_eta[iE][ict][ip], ang_sig same; the eps_tot `tot_points` seeded
    // from the backbone above gets the angular samples appended in the loop.
    const size_t nE_s = shape_E.size();
    std::vector<std::vector<double>> ang_ln_eta(nE_s), ang_sig(nE_s);
    // The eps_tot residual (eta_tot = eps_tot / K) measured on the same scan;
    // used when the tier auto-choice lands on the eta_tot table (HPGe-class -
    // its ANGULAR shape differs from the FEP one via dead-layer scatter-in).
    std::vector<std::vector<double>> ang_ln_tot(nE_s), ang_tot_sig(nE_s);
    uint32_t ang_node_id = 0;
    for (size_t ie = 0; ie < nE_s; ++ie) {
        ang_ln_eta[ie].resize(scan_ct.size() * np);
        ang_sig[ie].resize(scan_ct.size() * np);
        ang_ln_tot[ie].resize(scan_ct.size() * np);
        ang_tot_sig[ie].resize(scan_ct.size() * np);
        for (size_t ic = 0; ic < scan_ct.size(); ++ic) {
            for (size_t ip = 0; ip < np; ++ip) {
                const double phi_rad =
                    box ? phi_nodes[ip] * kPi / 180.0 : 0.0;
                const Eigen::Vector3d src =
                    source_position(d_far, scan_ct[ic], phi_rad);
                const EfficiencyResult r =
                    run.run_node(shape_E[ie], src, /*stage=*/2, ang_node_id++);
                const ApertureQuadrature q = resp->make_quadrature(src);
                const double K =
                    resp->kernel_K(shape_E[ie], q, MuChoice::Total);
                const double eps = std::max(r.full_energy_peak_efficiency, 1e-300);
                ang_ln_eta[ie][ic * np + ip] = std::log(eps / std::max(K, 1e-300));
                ang_sig[ie][ic * np + ip] =
                    r.fep_uncertainty / eps;
                const double tot = std::max(r.total_efficiency, 1e-300);
                ang_ln_tot[ie][ic * np + ip] = std::log(tot / std::max(K, 1e-300));
                ang_tot_sig[ie][ic * np + ip] = r.total_uncertainty / tot;
                if (r.total_efficiency > 0.0) {
                    const double K_nors =
                        resp->kernel_K(shape_E[ie], q, MuChoice::NoRayleigh);
                    tot_points.push_back(
                        {shape_E[ie], r.total_efficiency,
                         r.total_uncertainty / r.total_efficiency, K_nors, K});
                }
                run.tick("Angular response");
            }
        }
    }

    // Greedy shared cos-theta nodes (phi = 0 curves).
    std::vector<std::vector<double>> ct_curves, ct_sig_curves;
    for (size_t ie = 0; ie < nE_s; ++ie) {
        std::vector<double> c(scan_ct.size()), s(scan_ct.size());
        for (size_t ic = 0; ic < scan_ct.size(); ++ic) {
            c[ic] = ang_ln_eta[ie][ic * np + 0];
            s[ic] = ang_sig[ie][ic * np + 0];
        }
        ct_curves.push_back(std::move(c));
        ct_sig_curves.push_back(std::move(s));
    }
    const std::vector<size_t> ct_nodes = greedy_ct_nodes(
        scan_ct, ct_curves, ct_sig_curves, static_cast<size_t>(opts.n_cos_theta_nodes),
        opts.node_interp_tol, opts.node_noise_k);

    // ---- assemble the eta_fep tensor ---------------------------------------
    ShapeInterp shape_interp;
    shape_interp.set(shape_E);

    EtaTable& t = resp->eta_fep;
    for (const size_t n : e_nodes) t.energies_keV.push_back(scan_E[n]);
    for (const size_t c : ct_nodes) t.cos_thetas.push_back(scan_ct[c]);
    if (box) t.phis_deg = phi_nodes;
    const size_t tnc = t.cos_thetas.size();
    const size_t tnp = box ? phi_nodes.size() : 1;
    t.ln_eta.assign(t.energies_keV.size() * tnc * tnp, 0.0);
    t.frac_sigma.assign(t.ln_eta.size(), 0.0);

    // Angular shape relative to (ct = 1, phi = 0), per shape energy.
    auto shape_at = [&](size_t ie, size_t ic_scan, size_t ip) {
        return ang_ln_eta[ie][ic_scan * np + ip] -
               ang_ln_eta[ie][(scan_ct.size() - 1) * np + 0];
    };

    for (size_t e = 0; e < t.energies_keV.size(); ++e) {
        const double E = t.energies_keV[e];
        const double bb = smoothed_node_value(e_nodes[e], e_nodes, scan_E,
                                              bb_ln_eta, bb_sig, edges);
        const double bb_s = bb_sig[e_nodes[e]];
        size_t i0, i1;
        double w;
        shape_interp.bracket(E, i0, i1, w);
        for (size_t c = 0; c < tnc; ++c) {
            const size_t ic_scan = ct_nodes[c];
            for (size_t p = 0; p < tnp; ++p) {
                const double sh = (1.0 - w) * shape_at(i0, ic_scan, p) +
                                  w * shape_at(i1, ic_scan, p);
                const double sh_s =
                    std::max(ang_sig[i0][ic_scan * np + p],
                             ang_sig[i1][ic_scan * np + p]);
                t.ln_eta[t.index(e, c, p)] = bb + sh;
                t.frac_sigma[t.index(e, c, p)] =
                    std::sqrt(bb_s * bb_s + sh_s * sh_s);
            }
        }
    }

    // ---- 3. near-field model ------------------------------------------------
    if (opts.profile != ResponseProfile::FarField) {
        // Tabulated ln N on a direct (cos_theta x distance) tensor grid per
        // shape energy (PCHIP-interpolated at eval). The prior (z_over_a,theta)
        // scan never sampled grazing-at-small-d -- z = d*cos(theta) collapses,
        // so at contact only near-grazing angles were reached -- which is
        // exactly where the worst probes lived (35-60 keV, 70-89 deg, d<10cm).
        // Angular nodes are dense toward grazing and reach the same edge as the
        // eta grid (cos_theta 0.02 ~ 88.85 deg); distance nodes are geometric
        // (suits PCHIP in ln d) with an outer ln N = 0 anchor row (no MC) so N
        // fades smoothly to 1 above the breakpoint.
        const std::vector<double> near_cts{0.02, 0.06, 0.105, 0.208, 0.342,
                                           0.53, 0.766, 0.94, 1.0};
        std::vector<double> d_over_a{0.12, 0.3, 0.6, 1.0, 1.4, 1.8, 2.4, 3.0};
        if (opts.profile == ResponseProfile::Contact)
            d_over_a.insert(d_over_a.begin(), 0.06);
        const double d_anchor_over_a = 4.0;  // ln N = 0 anchor row, no MC

        // Re-finalize the (now full) eta table so the far prediction works.
        resp->eta_fep.finalize();

        NearFieldModel& nf = resp->near_field;
        nf.energies_keV = shape_E;
        nf.cos_thetas = near_cts;
        for (const double dr : d_over_a) nf.dists_cm.push_back(dr * a);
        nf.dists_cm.push_back(d_anchor_over_a * a);
        const size_t nnc = nf.cos_thetas.size();
        const size_t nnd = nf.dists_cm.size();
        nf.ln_n.assign(nE_s * nnc * nnd, 0.0);      // anchor row stays 0
        nf.frac_sigma.assign(nf.ln_n.size(), 0.002);  // anchor sigma

        uint32_t near_id = 0;
        for (size_t ie = 0; ie < nE_s; ++ie) {
            const double E = shape_E[ie];
            for (size_t ic = 0; ic < nnc; ++ic) {
                const double ct = nf.cos_thetas[ic];
                for (size_t id = 0; id + 1 < nnd; ++id) {  // skip anchor row
                    const double d = nf.dists_cm[id];
                    const Eigen::Vector3d src = source_position(d, ct);
                    // Stage-3 nodes consume ONLY eps_fep (eps_tot below is
                    // never read), so fep_only transport would be free CPU --
                    // but it is DISABLED: bench_fep_only_stage3 (2026-07)
                    // measured a -5..-11% eps_fep bias on bore-hole detectors
                    // (kill-on-scoring-exit kills photons crossing the vacuum
                    // bore that re-enter with full energy). Re-enable per
                    // detector class only after that kill rule handles
                    // void-gap re-entry.
                    const EfficiencyResult r =
                        run.run_node(E, src, /*stage=*/3, near_id++,
                                     /*fep_only=*/false);
                    const ApertureQuadrature q = resp->make_quadrature(src);
                    const double K = resp->kernel_K(E, q, MuChoice::Total);
                    const size_t idx = nf.index(ie, ic, id);
                    if (r.full_energy_peak_efficiency <= 0.0 || K <= 0.0) {
                        // A table cannot skip a node: hold the previous
                        // (smaller-d) value so PCHIP is not kinked; inflate its
                        // sigma honestly.
                        nf.ln_n[idx] =
                            (id > 0) ? nf.ln_n[nf.index(ie, ic, id - 1)] : 0.0;
                        nf.frac_sigma[idx] = 0.75;
                        run.tick("Near field");
                        continue;
                    }
                    bool clamped = false;
                    const double ln_eta_far =
                        resp->eta_fep.eval_ln(E, ct, 0.0, clamped);
                    nf.ln_n[idx] =
                        std::log(r.full_energy_peak_efficiency / K) - ln_eta_far;
                    nf.frac_sigma[idx] =
                        r.fep_uncertainty / r.full_energy_peak_efficiency;
                    run.tick("Near field");
                }
            }
        }
        nf.finalize();

        // Gate: apply N for all d < the anchor distance (the outermost distance
        // node, where ln N is pinned to 0). The table is measured -- not
        // extrapolated -- everywhere inside and fades smoothly to 1 at the
        // anchor, so a single distance gate suffices. (A per-(E,ct) MEASURED
        // 1%-breakpoint scan was tried and reverted: it only sampled ct in
        // {0.5,1.0}, so extreme grazing clamped to the ct=0.5 breakpoint and
        // switched N off far too early -- an ~9% over-prediction where the
        // grazing near-field correction was still large -- and the |lnN|>1%
        // crossing was bimodal/noisy across energy.) break_cos_thetas is kept
        // as a 2-node grid so breakpoint_d_cm's bilinear lookup returns the
        // constant anchor distance for any cos_theta.
        nf.break_cos_thetas = {0.02, 1.0};
        nf.break_d_cm.assign(nf.energies_keV.size() * 2, nf.dists_cm.back());
    } else {
        // Far-field profile: validity floor at the S1-measured kernel
        // breakpoint envelope (1.1-1.8 a) -- conservative 2 a.
        resp->provenance.min_distance_cm = 2.0 * a;
    }

    // ---- 4. eps_tot tier auto-choice ----------------------------------------
    if (decide_tot_tier(*resp, tot_points, opts) == TotDecision::NeedTable) {
        // eta_tot table: same assembly as eta_fep but for eps_tot / K.
        resp->tot_eff.tier = TotEffTier::EtaTotTable;
        {
            EtaTable& tt = resp->tot_eff.eta_tot;
            tt.energies_keV = t.energies_keV;
            tt.cos_thetas = t.cos_thetas;
            tt.phis_deg = t.phis_deg;
            tt.edges_keV = t.edges_keV;
            tt.ln_eta.assign(t.ln_eta.size(), 0.0);
            tt.frac_sigma.assign(t.ln_eta.size(), 0.0);

            // Backbone eta_tot (on-axis, from the dense energy scan) +
            // MEASURED eps_tot angular shape from the same scans as the FEP
            // shape (the eps_tot shape differs - dead-layer scatter-in is
            // angle-dependent for HPGe-class detectors).
            std::vector<double> bb_ln_tot(scan_E.size());
            for (size_t i = 0; i < scan_E.size(); ++i) {
                const double K =
                    resp->kernel_K(scan_E[i], q_far, MuChoice::Total);
                bb_ln_tot[i] = std::log(std::max(bb_tot[i], 1e-300) /
                                        std::max(K, 1e-300));
            }

            auto tot_shape_at = [&](size_t ie, size_t ic_scan, size_t ip) {
                return ang_ln_tot[ie][ic_scan * np + ip] -
                       ang_ln_tot[ie][(scan_ct.size() - 1) * np + 0];
            };

            for (size_t e = 0; e < tt.energies_keV.size(); ++e) {
                size_t i0, i1;
                double w;
                shape_interp.bracket(tt.energies_keV[e], i0, i1, w);
                for (size_t c = 0; c < tnc; ++c) {
                    const size_t ic_scan = ct_nodes[c];
                    for (size_t p = 0; p < tnp; ++p) {
                        const double sh = (1.0 - w) * tot_shape_at(i0, ic_scan, p) +
                                          w * tot_shape_at(i1, ic_scan, p);
                        const double sh_s =
                            std::max(ang_tot_sig[i0][ic_scan * np + p],
                                     ang_tot_sig[i1][ic_scan * np + p]);
                        tt.ln_eta[tt.index(e, c, p)] =
                            bb_ln_tot[e_nodes[e]] + sh;
                        tt.frac_sigma[tt.index(e, c, p)] = std::sqrt(
                            bb_tot_sig[e_nodes[e]] * bb_tot_sig[e_nodes[e]] +
                            sh_s * sh_s);
                    }
                }
            }
        }
    }

    // EFFTRAN few-anchor transfer: the coarse eta(E,theta) still under-resolves
    // the angular residual, so attach the honest off-axis/near sigma envelope.
    if (opts.transfer_mode) {
        resp->model_transfer = SigmaTransferModel{};
        resp->scatter_in_recapture = kTotalScatterInRecapture;
    }

    resp->finalize();
    if (opts.progress) opts.progress(1.0, "Done");
    return resp;
}

// ---------------------------------------------------------------------------

void ResponseGenerator::ground_to_points(DetectorResponse& resp,
                                         std::vector<GroundingPoint> points,
                                         bool curve_derived) {
    if (points.empty())
        throw std::runtime_error("ground_to_points: no points");

    // Model efficiencies with any existing grounding REMOVED (k is a ratio
    // to the ungrounded model).
    resp.grounding = GroundingBlock();
    for (GroundingPoint& p : points) {
        if (p.model_eff > 0.0) continue;
        const double theta = std::acos(std::min(std::max(p.cos_theta, -1.0), 1.0));
        const EffResult r = resp.eps_fep(p.energy_keV, theta,
                                         p.phi_deg * kPi / 180.0, p.distance_cm);
        p.model_eff = r.value;
    }

    // Knots: 2-4 hats across the measured ln-E span (n scales with the
    // number of distinct energies; a near-flat ratio needs few parameters).
    std::vector<double> ln_es;
    for (const GroundingPoint& p : points) ln_es.push_back(std::log(p.energy_keV));
    std::sort(ln_es.begin(), ln_es.end());
    ln_es.erase(std::unique(ln_es.begin(), ln_es.end(),
                            [](double x, double y) { return y - x < 0.01; }),
                ln_es.end());
    const int n_knots =
        std::max(1, std::min({4, static_cast<int>(ln_es.size()),
                              static_cast<int>(points.size())}));
    std::vector<double> knots;
    if (n_knots == 1) {
        knots.push_back(ln_es.front());
    } else {
        for (int i = 0; i < n_knots; ++i)
            knots.push_back(ln_es.front() +
                            (ln_es.back() - ln_es.front()) * double(i) /
                                (n_knots - 1));
    }

    const size_t n = points.size();
    Eigen::MatrixXd X(n, n_knots);
    Eigen::VectorXd y(n);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
    for (size_t i = 0; i < n; ++i) {
        const GroundingPoint& p = points[i];
        if (p.model_eff <= 0.0 || p.measured_eff <= 0.0)
            throw std::runtime_error("ground_to_points: non-positive efficiency");
        const std::vector<double> B = hat_basis(std::log(p.energy_keV), knots);
        for (int j = 0; j < n_knots; ++j) X(i, j) = B[j];
        y(i) = std::log(p.measured_eff / p.model_eff);
        V(i, i) = p.frac_stat_sigma * p.frac_stat_sigma + 1e-8;
    }
    // Per-source certificate blocks: fully correlated within a source.
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (!points[i].source_key.empty() &&
                points[i].source_key == points[j].source_key)
                V(i, j) += points[i].frac_cert_sigma * points[j].frac_cert_sigma;

    const Eigen::MatrixXd Vinv = V.llt().solve(Eigen::MatrixXd::Identity(n, n));
    const Eigen::MatrixXd XtVi = X.transpose() * Vinv;
    const Eigen::MatrixXd C = (XtVi * X).ldlt().solve(
        Eigen::MatrixXd::Identity(n_knots, n_knots));
    const Eigen::VectorXd c = C * (XtVi * y);

    GroundingBlock& g = resp.grounding;
    g.curve_derived = curve_derived;
    g.points = std::move(points);
    g.knot_ln_energies = knots;
    g.ln_k.assign(c.data(), c.data() + n_knots);
    g.cov.resize(static_cast<size_t>(n_knots) * n_knots);
    for (int i = 0; i < n_knots; ++i)
        for (int j = 0; j < n_knots; ++j)
            g.cov[static_cast<size_t>(i) * n_knots + j] = C(i, j);
}

// ---------------------------------------------------------------------------

std::vector<ProbePoint> ResponseGenerator::probe_bank(
    const GeometryDescriptor& gd, const GenerationOptions& opts, int n_points,
    int start_index, double d_min_cm, double d_max_cm) {
    Runner run(gd, opts);
    run.nodes_total = n_points;
    const bool box = gd.shape == DetectorShape::Box;
    const double phi_span =
        box ? ((gd.dimensions_cm.size() > 1 &&
                std::fabs(gd.dimensions_cm[0] - gd.dimensions_cm[1]) < 1e-9)
                   ? 45.0 : 90.0)
            : 0.0;
    const std::vector<double> edges =
        crystal_edges(gd, opts.e_min_keV, opts.e_max_keV);

    std::vector<ProbePoint> bank;
    bank.reserve(static_cast<size_t>(n_points));
    for (int p0 = 0; p0 < n_points; ++p0) {
        const int p = p0 + start_index;
        const double uE = halton(static_cast<uint64_t>(p), 2);
        const double ud = halton(static_cast<uint64_t>(p), 3);
        const double uc = halton(static_cast<uint64_t>(p), 5);
        const double up = halton(static_cast<uint64_t>(p), 7);

        double E = opts.e_min_keV *
                   std::pow(opts.e_max_keV / opts.e_min_keV, uE);
        if (p % 10 == 9 && !edges.empty()) {
            const double edge = edges[static_cast<size_t>(p / 10) % edges.size()];
            const double off = ((p / 10) % 2 == 0) ? -3.0 : 3.0;
            E = std::max(opts.e_min_keV, edge + off);
        }
        const double d = d_min_cm * std::pow(d_max_cm / d_min_cm, ud);
        const double ct = opts.cos_theta_min + (1.0 - opts.cos_theta_min) * uc;
        const double phi_deg = box ? phi_span * up : 0.0;
        const Eigen::Vector3d src = source_position(d, ct, phi_deg * kPi / 180.0);

        ProbePoint row;
        row.energy_keV = E;
        row.d_cm = d;
        row.cos_theta = ct;
        row.phi_deg = phi_deg;
        row.seed = node_seed(opts.base_seed, /*stage=*/9,
                             static_cast<uint32_t>(p));
        EfficiencyCalculator calc;
        run.configure(calc);
        calc.set_point_source(src);
        SimulationConfig cfg;
        cfg.energy_keV = E;
        // Probe/validation banks always run at the uniform scalar precision --
        // never the generation graded map (they measure, not build).
        apply_node_budget(cfg, opts, opts.node_fep_precision);
        cfg.seed = row.seed;
        const EfficiencyResult r = calc.compute(cfg);
        row.eps_fep = r.full_energy_peak_efficiency;
        row.fep_unc = r.fep_uncertainty;
        row.eps_tot = r.total_efficiency;
        row.tot_unc = r.total_uncertainty;
        bank.push_back(row);
        // Optional cost accounting (e.g. the certificate probe MC cost). Stage 9
        // marks probe-bank nodes; only active when the caller supplies stats_out.
        if (opts.stats_out) {
            NodeStat ns;
            ns.stage = 9;
            ns.energy_keV = E;
            ns.d_cm = d;
            ns.cos_theta = ct;
            ns.events = r.num_events_simulated;
            ns.cpu_s = r.cpu_time_seconds;
            ns.wall_s = r.wall_time_seconds;
            ns.stop = static_cast<uint8_t>(r.stop_reason);
            ns.fep_rel_prec = (r.full_energy_peak_efficiency > 0.0)
                ? r.fep_uncertainty / r.full_energy_peak_efficiency : 0.0;
            opts.stats_out->add(ns);
        }
        run.tick("Probe bank");
    }
    return bank;
}

std::vector<ResponseGenerator::StructuredProbe>
ResponseGenerator::plan_structured_probes(
    const GeometryDescriptor& gd, const GenerationOptions& opts,
    const DetectorResponse& resp, ProbeFamilyMask families, int n_per_family,
    int start_index, double d_min_cm, double d_max_cm) {
    std::vector<StructuredProbe> out;
    if (n_per_family <= 0) return out;

    const bool box = gd.shape == DetectorShape::Box;
    const double phi_span =
        box ? ((gd.dimensions_cm.size() > 1 &&
                std::fabs(gd.dimensions_cm[0] - gd.dimensions_cm[1]) < 1e-9)
                   ? 45.0 : 90.0)
            : 0.0;
    const double phi_mid = box ? 0.5 * phi_span : 0.0;

    const double a = gd.transverse_half_extent();
    const double d_far = std::max(opts.far_distance_a * a, 10.0);
    const double d_near06 = std::max(0.6 * a, std::max(d_min_cm, 0.5));

    const std::vector<double>& En = resp.eta_fep.energies_keV;
    const std::vector<double>& Cn = resp.eta_fep.cos_thetas;
    const std::vector<double>& edges = resp.eta_fep.edges_keV;

    // EnergyGap probes run at NODE precision so a failing probe's MC is directly
    // reusable as the inserted backbone node; the rest measure at the scalar.
    const double eg_prec = std::min(opts.node_fep_precision, 0.003);
    const double st_prec = opts.node_fep_precision;

    // Shape (angular-scan carrier) energies: the near-field table's E-axis when
    // present, else a few backbone samples (far-field profile has no near table).
    std::vector<double> shapeE = resp.near_field.energies_keV;
    if (shapeE.empty() && !En.empty()) {
        const int ns = std::min<int>(5, static_cast<int>(En.size()));
        for (int i = 0; i < ns; ++i)
            shapeE.push_back(En[static_cast<size_t>(
                (En.size() - 1) * i / std::max(1, ns - 1))]);
        std::sort(shapeE.begin(), shapeE.end());
        shapeE.erase(std::unique(shapeE.begin(), shapeE.end()), shapeE.end());
    }

    auto want = [&](ProbeFamily f) {
        return (families & probe_family_bit(f)) != 0;
    };
    struct Gap { double v, w; int kind; };  // value, gap-width, kind
    auto widest = [](std::vector<Gap>& g) {
        std::sort(g.begin(), g.end(),
                  [](const Gap& a, const Gap& b) { return a.w > b.w; });
    };

    // ---- EnergyGap: geom-mean midpoints of adjacent backbone nodes (in-segment)
    if (want(ProbeFamily::EnergyGap) && En.size() >= 2) {
        std::vector<Gap> g;
        for (size_t i = 0; i + 1 < En.size(); ++i) {
            if (!edges.empty() &&
                segment_of(En[i], edges) != segment_of(En[i + 1], edges))
                continue;
            g.push_back({std::sqrt(En[i] * En[i + 1]),
                         std::log(En[i + 1] / En[i]), 0});
        }
        widest(g);
        const int n = std::min<int>(n_per_family, static_cast<int>(g.size()));
        for (int i = 0; i < n; ++i)
            out.push_back({g[i].v, d_far, 1.0, 0.0, ProbeFamily::EnergyGap,
                           eg_prec});
    }

    // ---- CtGap: cos-theta midpoints at shape-E node energies, far field
    if (want(ProbeFamily::CtGap) && Cn.size() >= 2 && !shapeE.empty()) {
        std::vector<Gap> g;
        for (size_t i = 0; i + 1 < Cn.size(); ++i)
            g.push_back({0.5 * (Cn[i] + Cn[i + 1]), Cn[i + 1] - Cn[i], 0});
        widest(g);
        const int n = std::min<int>(n_per_family, static_cast<int>(g.size()));
        for (int i = 0; i < n; ++i)
            out.push_back({shapeE[static_cast<size_t>(i) % shapeE.size()], d_far,
                           g[i].v, phi_mid, ProbeFamily::CtGap, st_prec});
    }

    // ---- ShapeEGap: energies between adjacent shape rows, mid-ct, far + d~0.6a
    if (want(ProbeFamily::ShapeEGap) && shapeE.size() >= 2) {
        std::vector<Gap> g;
        for (size_t i = 0; i + 1 < shapeE.size(); ++i) {
            if (!edges.empty() &&
                segment_of(shapeE[i], edges) != segment_of(shapeE[i + 1], edges))
                continue;
            g.push_back({std::sqrt(shapeE[i] * shapeE[i + 1]),
                         std::log(shapeE[i + 1] / shapeE[i]), 0});
        }
        widest(g);
        int pushed = 0;
        for (size_t i = 0; i < g.size() && pushed < n_per_family; ++i) {
            for (const double d : {d_far, d_near06}) {
                if (pushed >= n_per_family) break;
                out.push_back({g[i].v, d, 0.4, phi_mid, ProbeFamily::ShapeEGap,
                               st_prec});
                ++pushed;
            }
        }
    }

    // ---- NearCell: near-field d-gap / ct-gap midpoints (other axes at a node)
    if (want(ProbeFamily::NearCell) && !resp.near_field.empty()) {
        const NearFieldModel& nf = resp.near_field;
        const std::vector<double>& D = nf.dists_cm;
        const std::vector<double>& NC = nf.cos_thetas;
        const std::vector<double>& NE = nf.energies_keV;
        const double ct_rep = NC.empty() ? 1.0 : NC[NC.size() / 2];
        const double E_rep = NE.empty()
            ? (En.empty() ? 662.0 : En[En.size() / 2]) : NE[NE.size() / 2];
        const double d_rep = D.empty() ? d_near06 : D[D.size() / 2];
        std::vector<Gap> g;
        for (size_t i = 0; i + 1 < D.size(); ++i)  // kind 0 = d-gap
            g.push_back({std::sqrt(D[i] * D[i + 1]),
                         std::log(D[i + 1] / std::max(D[i], 1e-9)), 0});
        for (size_t i = 0; i + 1 < NC.size(); ++i)  // kind 1 = ct-gap
            g.push_back({0.5 * (NC[i] + NC[i + 1]), NC[i + 1] - NC[i], 1});
        widest(g);
        const int n = std::min<int>(n_per_family, static_cast<int>(g.size()));
        for (int i = 0; i < n; ++i) {
            if (g[i].kind == 0)
                out.push_back({E_rep, g[i].v, ct_rep, phi_mid,
                               ProbeFamily::NearCell, st_prec});
            else
                out.push_back({E_rep, d_rep, g[i].v, phi_mid,
                               ProbeFamily::NearCell, st_prec});
        }
    }

    // ---- Random: Halton canaries (mirrors probe_bank's sampler)
    if (want(ProbeFamily::Random)) {
        const std::vector<double> pedges =
            crystal_edges(gd, opts.e_min_keV, opts.e_max_keV);
        for (int i = 0; i < n_per_family; ++i) {
            const int pi = start_index + i;
            const double uE = halton(static_cast<uint64_t>(pi), 2);
            const double ud = halton(static_cast<uint64_t>(pi), 3);
            const double uc = halton(static_cast<uint64_t>(pi), 5);
            const double up = halton(static_cast<uint64_t>(pi), 7);
            double E = opts.e_min_keV *
                       std::pow(opts.e_max_keV / opts.e_min_keV, uE);
            if (pi % 10 == 9 && !pedges.empty()) {
                const double edge =
                    pedges[static_cast<size_t>(pi / 10) % pedges.size()];
                const double off = ((pi / 10) % 2 == 0) ? -3.0 : 3.0;
                E = std::max(opts.e_min_keV, edge + off);
            }
            const double d = d_min_cm * std::pow(d_max_cm / d_min_cm, ud);
            const double ct = opts.cos_theta_min + (1.0 - opts.cos_theta_min) * uc;
            const double phi = box ? phi_span * up : 0.0;
            out.push_back({E, d, ct, phi, ProbeFamily::Random, st_prec});
        }
    }
    return out;
}

std::vector<ProbePoint> ResponseGenerator::structured_probe_bank(
    const GeometryDescriptor& gd, const GenerationOptions& opts,
    const DetectorResponse& resp, ProbeFamilyMask families, int n_per_family,
    int start_index, double d_min_cm, double d_max_cm) {
    const std::vector<StructuredProbe> plan = plan_structured_probes(
        gd, opts, resp, families, n_per_family, start_index, d_min_cm, d_max_cm);
    Runner run(gd, opts);
    run.nodes_total = std::max<size_t>(1, plan.size());
    std::vector<ProbePoint> bank;
    bank.reserve(plan.size());
    for (size_t i = 0; i < plan.size(); ++i) {
        const StructuredProbe& sp = plan[i];
        // Deterministic seed: stage 10 (probe bank), node = start_index + i so
        // distinct start_index per iteration gives a disjoint, reproducible set.
        const uint64_t seed = node_seed(opts.base_seed, /*stage=*/10,
                                        static_cast<uint32_t>(start_index +
                                                              static_cast<int>(i)));
        bank.push_back(run_probe_point(run, opts, sp.energy_keV, sp.d_cm,
                                       sp.cos_theta, sp.phi_deg, sp.node_precision,
                                       seed, static_cast<uint8_t>(sp.family),
                                       /*stats stage=*/10));
        run.tick("Structured probes");
    }
    return bank;
}

std::vector<ProbePoint> ResponseGenerator::adversarial_probe_bank(
    const GeometryDescriptor& gd, const GenerationOptions& opts,
    const DetectorResponse& resp, int n_points, int start_index,
    double d_min_cm, double d_max_cm) {
    // Thin wrapper: the EnergyGap/CtGap/NearCell families reproduce the
    // interpolation-gap stresses this used to place by hand. Truncate to
    // n_points total (planner emits up to n_points per family).
    const ProbeFamilyMask fam = probe_family_bit(ProbeFamily::EnergyGap) |
                                probe_family_bit(ProbeFamily::CtGap) |
                                probe_family_bit(ProbeFamily::NearCell);
    std::vector<ProbePoint> bank = structured_probe_bank(
        gd, opts, resp, fam, n_points, start_index, d_min_cm, d_max_cm);
    if (static_cast<int>(bank.size()) > n_points)
        bank.resize(static_cast<size_t>(n_points));
    return bank;
}

void ResponseGenerator::certify(DetectorResponse& resp,
                                const GeometryDescriptor& gd,
                                const GenerationOptions& opts, int n_probes,
                                int seed_offset, ProbeFamilyMask cert_families) {
    // Fresh quasi-random probe bank at a fixed UNIFORM precision (0.005),
    // Halton offset seed_offset (7000 -- disjoint from generation's 5000 and
    // never parity-split). Uniform so the certificate MC never inherits the
    // graded generation map.
    GenerationOptions probe_opts = opts;
    probe_opts.node_precision = nullptr;
    probe_opts.precision_profile = PrecisionProfile::Uniform;
    probe_opts.node_fep_precision = 0.005;
    probe_opts.progress = nullptr;
    probe_opts.cancel = nullptr;
    GenerationStats stats;
    probe_opts.stats_out = &stats;

    std::vector<ProbePoint> bank =
        probe_bank(gd, probe_opts, n_probes, seed_offset);

    // D-b: append the structured families so the certificate's p95/max reflect
    // the worst-case interpolation gaps a random bank under-samples. Each row
    // carries its ProbeFamily tag. 0 (D-a default) keeps this random-only.
    const ProbeFamilyMask sfam = cert_families & kAllStructuredFamilies;
    if (sfam) {
        const std::vector<ProbePoint> sbank = structured_probe_bank(
            gd, probe_opts, resp, sfam, /*n_per_family=*/8, seed_offset);
        bank.insert(bank.end(), sbank.begin(), sbank.end());
    }

    AccuracyCertificate cert;
    cert.probe_seed_base = probe_opts.base_seed;
    cert.iterations = 0;               // caller (closed loop) overrides
    cert.cpu_seconds = stats.total_cpu_s;
    // `converged` is set AFTER scoring, below: it honestly means the certificate
    // met the accuracy gate (fep_p95 <= cert_tol), NOT merely that the probe
    // pass ran. The closed loop further ANDs its own no-actionable-failures
    // condition. A response that certifies above the gate reports converged =
    // false with its true p95 (never a misleading pass).

    std::vector<double> fep_errs, tot_errs;
    for (const ProbePoint& pp : bank) {
        const Eigen::Vector3d src =
            source_position(pp.d_cm, pp.cos_theta, pp.phi_deg * kPi / 180.0);
        const EffResult m = resp.eps_fep_at(pp.energy_keV, src);

        AccuracyCertificate::Row row;
        row.E_keV = pp.energy_keV;
        row.d_cm = pp.d_cm;
        row.cos_theta = pp.cos_theta;
        row.phi_deg = pp.phi_deg;
        row.mc = pp.eps_fep;
        row.mc_sig = pp.fep_unc;
        row.model = m.value;
        row.model_sig = m.sigma;
        row.tag = pp.tag;              // ProbeFamily (0 random; D-b structured)

        // Score only converged probes (mc_sig/mc <= 5%); a noisy MC truth is
        // not a verdict on the model.
        double fep_err = 0.0;
        const bool fep_ok =
            pp.eps_fep > 0.0 && pp.fep_unc / pp.eps_fep <= 0.05;
        if (fep_ok) {
            fep_err = std::fabs(m.value / pp.eps_fep - 1.0);
            fep_errs.push_back(fep_err);
        }
        // Noise-aware pass: within tol OR within 2 sigma of the combined
        // MC + model noise.
        const double noise = 2.0 * std::sqrt(pp.fep_unc * pp.fep_unc +
                                             m.sigma * m.sigma);
        row.pass = fep_ok && (std::fabs(m.value - pp.eps_fep) <=
                              std::max(0.01 * pp.eps_fep, noise));
        cert.rows.push_back(row);

        if (pp.eps_tot > 0.0 && pp.tot_unc / pp.eps_tot <= 0.05) {
            const EffResult mt = resp.eps_total_at(pp.energy_keV, src);
            tot_errs.push_back(std::fabs(mt.value / pp.eps_tot - 1.0));
        }
    }

    auto percentile = [](std::vector<double> v, double q) -> double {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[static_cast<size_t>(q * (v.size() - 1))];
    };
    cert.fep_median = percentile(fep_errs, 0.50);
    cert.fep_p95 = percentile(fep_errs, 0.95);
    cert.fep_max = fep_errs.empty()
        ? 0.0 : *std::max_element(fep_errs.begin(), fep_errs.end());
    cert.tot_median = percentile(tot_errs, 0.50);
    cert.tot_p95 = percentile(tot_errs, 0.95);

    // Honest convergence: the certificate met the accuracy gate. Uses the
    // caller's cert_tol (GenerationOptions default 0.012); a response scored
    // above it reports converged = false with its true percentiles.
    cert.converged = !fep_errs.empty() && (cert.fep_p95 <= opts.cert_tol);

    resp.certificate = std::move(cert);
}

void ResponseGenerator::configure_calculator(
    EfficiencyCalculator& calc, const GeometryDescriptor& gd,
    std::vector<std::unique_ptr<Material>>& owned_materials) {
    // Validate before instantiating anything, so a bad descriptor cannot leave
    // half-built materials in the caller's vector (same ordering rule as
    // GeometryDescriptor::build_geometry).
    const std::vector<GeometryProblem> probs = gd.problems();
    if (!probs.empty())
        throw std::runtime_error(std::string("configure_calculator: ")
                                 + to_string(probs.front()));

    // Instantiate the descriptor's material table; the calculator keeps raw
    // pointers, so the caller owns their lifetime via owned_materials.
    const size_t base = owned_materials.size();
    owned_materials.reserve(base + gd.materials.size());
    for (const MaterialSpec& spec : gd.materials)
        owned_materials.push_back(std::make_unique<Material>(spec.to_material()));
    auto mat = [&](int idx) -> const Material* {
        if (idx < 0 || idx >= static_cast<int>(gd.materials.size()))
            throw std::runtime_error("configure_calculator: bad material index");
        return owned_materials[base + static_cast<size_t>(idx)].get();
    };

    calc.set_detector(gd.shape, mat(gd.crystal_material_index),
                      gd.dimensions_cm);
    // set_detector() clears the fillet/bore/dead layer, so declare them after
    // it; fillet first, so bore_fits() sees the final crystal profile.
    if (gd.bullet_radius_cm > 0.0) calc.set_bullet_radius(gd.bullet_radius_cm);
    if (gd.bore)
        calc.set_bore_hole(gd.bore->radius, gd.bore->depth,
                           gd.bore->rounded_tip);
    if (gd.dead_layer)
        calc.set_dead_layer(gd.dead_layer->front, gd.dead_layer->side,
                            gd.dead_layer->back);
    for (const LayerSpec& l : gd.layers)
        calc.add_attenuator(mat(l.material_index), l.front_thickness_cm,
                            l.side_thickness_cm, l.z_start_cm, l.z_end_cm);
    if (gd.collimator)
        calc.add_collimator(mat(gd.collimator->material_index),
                            gd.collimator->side_thickness_cm,
                            gd.collimator->z_start_cm, gd.collimator->z_end_cm);
}

} // namespace ceelo
