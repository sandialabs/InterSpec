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

#include "transport/FluenceTally.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace ceelo {

namespace {

double standard_error(uint64_t n_histories, double sum_w, double sum_w2)
{
    if (n_histories < 2) return 0.0;
    const double n = static_cast<double>(n_histories);
    const double centered = std::max(0.0, sum_w2 - sum_w * sum_w / n);
    return std::sqrt(centered / (n * (n - 1.0)));
}

} // namespace

int FluenceAxes::energy_bin(double e_keV) const
{
    if (!edges.empty()) {
        // Explicit-edge (landmark) axis: binary search. Out of range is
        // dropped exactly as in the log-uniform case, so `dropped_weight_`
        // keeps its meaning.
        if (e_keV < edges.front() || e_keV > edges.back()) return -1;
        const auto it = std::upper_bound(edges.begin(), edges.end(), e_keV);
        int i = static_cast<int>(it - edges.begin()) - 1;
        if (i < 0) i = 0;
        if (i >= n_energy) i = n_energy - 1;
        return i;
    }
    if (e_keV < e_min_keV || e_keV > e_max_keV) return -1;
    const double u = std::log(e_keV / e_min_keV) / std::log(e_max_keV / e_min_keV);
    int i = static_cast<int>(u * n_energy);
    if (i < 0) i = 0;
    if (i >= n_energy) i = n_energy - 1;
    return i;
}

bool FluenceAxes::same_axis(const FluenceAxes& o) const
{
    if (n_energy != o.n_energy || n_mu != o.n_mu) return false;
    if (edges.size() != o.edges.size()) return false;
    for (std::size_t i = 0; i < edges.size(); ++i)
        if (edges[i] != o.edges[i]) return false;
    if (edges.empty())
        return e_min_keV == o.e_min_keV && e_max_keV == o.e_max_keV;
    return true;
}

FluenceAxes make_registered_landmark_axes(double e0_keV, int n_energy_target,
                                         int n_mu, double e_min_keV)
{
    FluenceAxes ax;
    ax.n_mu = n_mu;
    ax.e_min_keV = e_min_keV;
    ax.e_max_keV = e0_keV;

    const double e_bs = e0_keV / (1.0 + 2.0 * e0_keV / 510.998950);
    const int half = std::max(2, (std::max(4, n_energy_target) + 1) / 2);

    if (!(e0_keV > e_min_keV) || !(e_bs > e_min_keV) || !(e0_keV > e_bs)) {
        // Below ~34 keV the backscatter edge falls under the floor and the
        // two-block structure degenerates; fall back to a plain log axis
        // rather than emit an unusable one.
        ax.n_energy = std::max(2, n_energy_target);
        return ax;
    }

    ax.edges.reserve(2 * half + 1);
    ax.edges.push_back(e_min_keV);
    for (int i = 1; i <= half; ++i)
        ax.edges.push_back(e_min_keV *
                           std::pow(e_bs / e_min_keV, double(i) / half));
    for (int i = 1; i <= half; ++i)
        ax.edges.push_back(e_bs * std::pow(e0_keV / e_bs, double(i) / half));
    // Pin the landmarks exactly: a consumer testing "is this bin below the
    // backscatter edge" must get a clean answer, and the top edge must equal
    // E0 exactly for the container's truncation invariant.
    ax.edges[half] = e_bs;
    ax.edges[2 * half] = e0_keV;
    ax.n_energy = 2 * half;
    return ax;
}

std::vector<FluenceAxes> make_fit_component_axes(double e0_keV,
                                                 int n_energy_target,
                                                 int n_mu,
                                                 double e_min_keV)
{
    std::vector<FluenceAxes> axes(static_cast<int>(FluxOrigin::Count));
    axes[static_cast<int>(FluxOrigin::Primary)] =
        make_registered_landmark_axes(e0_keV, n_energy_target, n_mu,
                                      e_min_keV);
    axes[static_cast<int>(FluxOrigin::Annihilation)] =
        make_registered_landmark_axes(510.998950, n_energy_target, n_mu,
                                      e_min_keV);

    FluenceAxes brems;
    brems.e_min_keV = e_min_keV;
    brems.e_max_keV = e0_keV;
    brems.n_energy = std::max(2, n_energy_target);
    brems.n_mu = n_mu;
    axes[static_cast<int>(FluxOrigin::Brems)] = std::move(brems);
    return axes;
}

FluenceAxes make_landmark_axes(double e0_keV, int n_energy_target, int n_mu,
                              double e_min_keV, double albedo_boost)
{
    FluenceAxes ax;
    ax.n_mu = n_mu;
    ax.e_min_keV = e_min_keV;
    ax.e_max_keV = e0_keV;

    // Block boundaries: e_min, backscatter edge, 511, E0 -- keeping only the
    // ones that are strictly inside the range and ordered.
    const double e_bs = e0_keV / (1.0 + 2.0 * e0_keV / 510.998950);
    std::vector<double> bounds{e_min_keV};
    auto add_bound = [&](double v) {
        if (v > bounds.back() * 1.02 && v < e0_keV * 0.98) bounds.push_back(v);
    };
    add_bound(e_bs);
    add_bound(510.998950);
    bounds.push_back(e0_keV);

    const int n_blocks = static_cast<int>(bounds.size()) - 1;
    if (n_blocks < 1 || !(e0_keV > e_min_keV)) {
        // Degenerate (E0 at or below the floor): fall back to the log-uniform
        // rule rather than emit an unusable axis.
        ax.e_max_keV = std::max(e0_keV, e_min_keV * 1.5);
        ax.n_energy = std::max(2, n_energy_target);
        return ax;
    }

    // Bins per block ~ log width, with the sub-backscatter-edge block boosted:
    // that decade holds the albedo peak and had the worst measured
    // misplacement, and it is where a detector does the least smearing of its
    // own (so the incident grid IS the delivered resolution).
    std::vector<double> weight(n_blocks);
    double wsum = 0.0;
    for (int b = 0; b < n_blocks; ++b) {
        weight[b] = std::log(bounds[b + 1] / bounds[b]);
        if (b == 0) weight[b] *= albedo_boost;
        wsum += weight[b];
    }

    std::vector<int> counts(n_blocks, 2);
    int assigned = 0;
    for (int b = 0; b < n_blocks; ++b) {
        counts[b] = std::max(2, static_cast<int>(
            std::lround(n_energy_target * weight[b] / wsum)));
        assigned += counts[b];
    }
    // Nudge the largest block so the total lands on the target exactly.
    if (assigned != n_energy_target) {
        int big = 0;
        for (int b = 1; b < n_blocks; ++b)
            if (counts[b] > counts[big]) big = b;
        counts[big] = std::max(2, counts[big] + (n_energy_target - assigned));
    }

    ax.edges.push_back(bounds.front());
    for (int b = 0; b < n_blocks; ++b) {
        const double lo = bounds[b], hi = bounds[b + 1];
        for (int i = 1; i <= counts[b]; ++i)
            ax.edges.push_back(
                lo * std::pow(hi / lo, static_cast<double>(i) / counts[b]));
    }
    // Kill floating-point drift so the landmark edges are EXACT: a consumer
    // that tests "is this bin below the backscatter edge" must get a clean
    // answer, and the last edge must equal E0 exactly for the truncation
    // invariant in the binary container.
    int at = 0;
    for (int b = 0; b < n_blocks; ++b) {
        at += counts[b];
        ax.edges[at] = bounds[b + 1];
    }
    ax.n_energy = static_cast<int>(ax.edges.size()) - 1;
    return ax;
}

int FluenceAxes::mu_bin(double mu) const
{
    if (mu < -1.0) mu = -1.0;
    if (mu > 1.0) mu = 1.0;
    int j = static_cast<int>((mu + 1.0) * 0.5 * n_mu);
    if (j < 0) j = 0;
    if (j >= n_mu) j = n_mu - 1;
    return j;
}

double FluenceAxes::energy_low(int i) const
{
    if (!edges.empty()) {
        if (i < 0) return edges.front();
        if (i >= static_cast<int>(edges.size())) return edges.back();
        return edges[i];
    }
    return e_min_keV * std::pow(e_max_keV / e_min_keV,
                                static_cast<double>(i) / n_energy);
}

double FluenceAxes::energy_high(int i) const { return energy_low(i + 1); }

double FluenceAxes::mu_low(int j) const
{
    return -1.0 + 2.0 * static_cast<double>(j) / n_mu;
}

double FluenceAxes::mu_high(int j) const { return mu_low(j + 1); }

FluenceTally::FluenceTally(const FluenceAxes& axes, std::vector<LineDef> lines,
                           std::vector<FluenceAxes> component_axes)
    : axes_(axes),
      lines_(std::move(lines)),
      sum_w_(axes.n_bins(), 0.0),
      sum_w2_(axes.n_bins(), 0.0),
      class_sum_w_(axes.n_bins() * static_cast<int>(ChainClass::Count), 0.0),
      class_sum_w2_(axes.n_bins() * static_cast<int>(ChainClass::Count), 0.0),
      origin_sum_w_(axes.n_bins() * static_cast<int>(FluxOrigin::Count), 0.0),
      hist_w_(axes.n_bins(), 0.0),
      hist_class_w_(axes.n_bins() * static_cast<int>(ChainClass::Count), 0.0)
{
    const int no = static_cast<int>(FluxOrigin::Count);
    if (!component_axes.empty() && static_cast<int>(component_axes.size()) != no)
        throw std::invalid_argument(
            "FluenceTally: component_axes must be empty or have one axis per "
            "FluxOrigin");
    for (int o = 0; o < no; ++o) {
        component_axes_[o] = component_axes.empty() ? axes : component_axes[o];
        if (component_axes_[o].n_mu != axes.n_mu)
            throw std::invalid_argument(
                "FluenceTally: all component axes must share the common mu axis");
    }

    std::size_t component_values = 0;
    for (int c = 0; c < kFluxComponentCount; ++c) {
        component_offsets_[c] = component_values;
        component_values += static_cast<std::size_t>(
            component_axes_[static_cast<int>(flux_component_origin(c))].n_bins());
    }
    component_offsets_[kFluxComponentCount] = component_values;
    component_sum_w_.assign(component_values, 0.0);
    component_sum_w2_.assign(component_values, 0.0);
    hist_component_w_.assign(component_values, 0.0);
    component_total_w_.assign(kFluxComponentCount, 0.0);
    component_total_w2_.assign(kFluxComponentCount, 0.0);
    hist_component_total_.assign(kFluxComponentCount, 0.0);

    const int nl = static_cast<int>(lines_.size());
    const int nc = static_cast<int>(ChainClass::Count);
    line_sum_w_.assign(nl * axes.n_mu, 0.0);
    line_sum_w2_.assign(nl * axes.n_mu, 0.0);
    line_class_sum_w_.assign(nl * axes.n_mu * nc, 0.0);
    line_class_sum_w2_.assign(nl * axes.n_mu * nc, 0.0);
    line_total_w_.assign(nl, 0.0);
    line_total_w2_.assign(nl, 0.0);
    hist_line_w_.assign(nl * axes.n_mu, 0.0);
    hist_line_class_w_.assign(nl * axes.n_mu * nc, 0.0);
    hist_line_total_.assign(nl, 0.0);

    const std::size_t line_component_values =
        static_cast<std::size_t>(kFluxComponentCount) * nl * axes.n_mu;
    line_component_sum_w_.assign(line_component_values, 0.0);
    line_component_sum_w2_.assign(line_component_values, 0.0);
    hist_line_component_w_.assign(line_component_values, 0.0);
}

void FluenceTally::add_line(int line, double mu_in, double weight, ChainClass cls,
                            FluxOrigin origin, PostCreationOrder order)
{
    if (!(weight > 0.0) || !std::isfinite(weight)) return;
    if (line < 0 || line >= static_cast<int>(lines_.size())) return;
    const int o = static_cast<int>(origin);
    const int ord = static_cast<int>(order);
    if (o < 0 || o >= static_cast<int>(FluxOrigin::Count) ||
        ord < 0 || ord >= static_cast<int>(PostCreationOrder::Count))
        return;
    const int imu = axes_.mu_bin(mu_in);
    const int idx = line * axes_.n_mu + imu;

    if (hist_line_w_[idx] == 0.0) hist_line_touched_.push_back(idx);
    hist_line_w_[idx] += weight;
    hist_line_total_[line] += weight;

    const int cidx = idx * static_cast<int>(ChainClass::Count) +
                     static_cast<int>(cls);
    if (hist_line_class_w_[cidx] == 0.0) hist_line_class_touched_.push_back(cidx);
    hist_line_class_w_[cidx] += weight;

    hist_origin_total_[o] += weight;

    const int component = flux_component_index(origin, order);
    const int line_cells = static_cast<int>(lines_.size()) * axes_.n_mu;
    const int component_idx = component * line_cells + idx;
    if (hist_line_component_w_[component_idx] == 0.0)
        hist_line_component_touched_.push_back(component_idx);
    hist_line_component_w_[component_idx] += weight;
    hist_component_total_[component] += weight;

    const int ie = axes_.energy_bin(lines_[line].energy_keV);
    add_coarse(ie, mu_in, weight);

    // Line scores belong to the same history total as grid scores, so
    // total_mean/total_sigma stay representation-independent.
    hist_total_ += weight;
}

void FluenceTally::add(double e_keV, double mu_in, double weight, ChainClass cls,
                       FluxOrigin origin, PostCreationOrder order)
{
    if (!(weight > 0.0) || !std::isfinite(weight)) return;
    const int o = static_cast<int>(origin);
    const int ord = static_cast<int>(order);
    if (o < 0 || o >= static_cast<int>(FluxOrigin::Count) ||
        ord < 0 || ord >= static_cast<int>(PostCreationOrder::Count))
        return;
    const int ie = axes_.energy_bin(e_keV);
    if (ie < 0) {
        dropped_weight_ += weight;
        hist_dropped_weight_ += weight;
        return;
    }
    const int imu = axes_.mu_bin(mu_in);
    const int idx = axes_.index(ie, imu);

    if (hist_w_[idx] == 0.0) hist_touched_.push_back(idx);
    hist_w_[idx] += weight;
    hist_total_ += weight;

    const int cidx = idx * static_cast<int>(ChainClass::Count) +
                     static_cast<int>(cls);
    if (hist_class_w_[cidx] == 0.0) hist_class_touched_.push_back(cidx);
    hist_class_w_[cidx] += weight;

    // Mean-only sums accumulate directly (buffering only matters for w2).
    origin_sum_w_[idx * static_cast<int>(FluxOrigin::Count) +
                  o] += weight;
    hist_origin_total_[o] += weight;

    const int component = flux_component_index(origin, order);
    const FluenceAxes& cax = component_axes_[o];
    const int cie = cax.energy_bin(e_keV);
    if (cie >= 0) {
        const int cimu = cax.mu_bin(mu_in);
        const std::size_t component_idx = component_offsets_[component] +
            static_cast<std::size_t>(cax.index(cie, cimu));
        if (hist_component_w_[component_idx] == 0.0)
            hist_component_touched_.push_back(
                static_cast<int>(component_idx));
        hist_component_w_[component_idx] += weight;
    } else {
        component_dropped_weight_[component] += weight;
    }
    hist_component_total_[component] += weight;
    add_coarse(ie, mu_in, weight);
}

void FluenceTally::add_coarse(int energy_bin, double mu_in, double weight)
{
    if (energy_bin >= 0 && axes_.n_energy > 0) {
        const int band = std::min(
            kNumEnergyBands - 1,
            energy_bin * kNumEnergyBands / axes_.n_energy);
        hist_coarse_[band] += weight;
    }

    const double mu = std::max(-1.0, std::min(1.0, mu_in));
    const double mu2 = mu * mu;
    const double legendre[4] = {
        mu,
        0.5 * (3.0 * mu2 - 1.0),
        0.5 * (5.0 * mu2 * mu - 3.0 * mu),
        0.125 * (35.0 * mu2 * mu2 - 30.0 * mu2 + 3.0)
    };
    for (int l = 0; l < 4; ++l)
        hist_coarse_[kNumEnergyBands + l] += weight * legendre[l];
}

void FluenceTally::end_history()
{
    end_history_impl(static_cast<int>(n_histories_ % kNumLogicalReplicas));
}

void FluenceTally::end_history(int logical_replica)
{
    if (logical_replica < 0 || logical_replica >= kNumLogicalReplicas)
        throw std::out_of_range(
            "FluenceTally::end_history: logical replica must be in [0, 31]");
    end_history_impl(logical_replica);
}

void FluenceTally::end_history_impl(int logical_replica)
{
    ++n_histories_;
    ++replica_histories_[logical_replica];
    // Fold the per-history sums in ONCE: multiple vertex scores of one primary
    // are correlated, so the square must be taken of the history total, not of
    // the individual scores.
    std::vector<float>& replica_continuum =
        replica_continuum_sum_w_[logical_replica];
    if (!hist_touched_.empty() && replica_continuum.empty())
        replica_continuum.assign(sum_w_.size(), 0.0f);
    for (int idx : hist_touched_) {
        const double w = hist_w_[idx];
        sum_w_[idx] += w;
        sum_w2_[idx] += w * w;
        replica_continuum[idx] += static_cast<float>(w);
        hist_w_[idx] = 0.0;
    }
    hist_touched_.clear();
    for (int cidx : hist_class_touched_) {
        const double w = hist_class_w_[cidx];
        class_sum_w_[cidx] += w;
        class_sum_w2_[cidx] += w * w;
        hist_class_w_[cidx] = 0.0;
    }
    hist_class_touched_.clear();

    std::vector<float>& replica_components =
        replica_component_sum_w_[logical_replica];
    if (!hist_component_touched_.empty() && replica_components.empty())
        replica_components.assign(component_sum_w_.size(), 0.0f);
    for (int component_idx : hist_component_touched_) {
        const double w = hist_component_w_[component_idx];
        component_sum_w_[component_idx] += w;
        component_sum_w2_[component_idx] += w * w;
        replica_components[component_idx] += static_cast<float>(w);
        hist_component_w_[component_idx] = 0.0;
    }
    hist_component_touched_.clear();

    std::vector<float>& replica_line_total =
        replica_line_sum_w_[logical_replica];
    if (!hist_line_touched_.empty() && replica_line_total.empty())
        replica_line_total.assign(line_sum_w_.size(), 0.0f);
    for (int idx : hist_line_touched_) {
        const double w = hist_line_w_[idx];
        line_sum_w_[idx] += w;
        line_sum_w2_[idx] += w * w;
        replica_line_total[idx] += static_cast<float>(w);
        hist_line_w_[idx] = 0.0;
    }
    hist_line_touched_.clear();
    for (int cidx : hist_line_class_touched_) {
        const double w = hist_line_class_w_[cidx];
        line_class_sum_w_[cidx] += w;
        line_class_sum_w2_[cidx] += w * w;
        hist_line_class_w_[cidx] = 0.0;
    }
    hist_line_class_touched_.clear();

    std::vector<float>& replica_lines =
        replica_line_component_sum_w_[logical_replica];
    if (!hist_line_component_touched_.empty() && replica_lines.empty())
        replica_lines.assign(line_component_sum_w_.size(), 0.0f);
    for (int component_idx : hist_line_component_touched_) {
        const double w = hist_line_component_w_[component_idx];
        line_component_sum_w_[component_idx] += w;
        line_component_sum_w2_[component_idx] += w * w;
        replica_lines[component_idx] += static_cast<float>(w);
        hist_line_component_w_[component_idx] = 0.0;
    }
    hist_line_component_touched_.clear();

    for (std::size_t l = 0; l < hist_line_total_.size(); ++l) {
        const double w = hist_line_total_[l];
        if (w == 0.0) continue;
        line_total_w_[l] += w;
        line_total_w2_[l] += w * w;
        hist_line_total_[l] = 0.0;
    }

    for (int o = 0; o < static_cast<int>(FluxOrigin::Count); ++o) {
        const double w = hist_origin_total_[o];
        if (w == 0.0) continue;
        origin_total_w_[o] += w;
        origin_total_w2_[o] += w * w;
        hist_origin_total_[o] = 0.0;
    }

    for (int c = 0; c < kFluxComponentCount; ++c) {
        const double w = hist_component_total_[c];
        component_total_w_[c] += w;
        component_total_w2_[c] += w * w;
        hist_component_total_[c] = 0.0;
    }

    bool coarse_nonzero = false;
    for (double x : hist_coarse_)
        coarse_nonzero = coarse_nonzero || x != 0.0;
    if (coarse_nonzero) {
        for (int i = 0; i < kCoarseDimension; ++i) {
            coarse_sum_w_[i] += hist_coarse_[i];
            for (int j = 0; j < kCoarseDimension; ++j)
                coarse_outer_sum_w_[i * kCoarseDimension + j] +=
                    hist_coarse_[i] * hist_coarse_[j];
        }
    }
    hist_coarse_.fill(0.0);

    // The integral's variance must come from the per-history total, not from
    // summing per-bin variances -- see total_sigma().
    sum_total_w_ += hist_total_;
    sum_total_w2_ += hist_total_ * hist_total_;
    max_history_total_ = std::max(max_history_total_, hist_total_);
    hist_total_ = 0.0;

    dropped_weight_sum_w2_ += hist_dropped_weight_ * hist_dropped_weight_;
    hist_dropped_weight_ = 0.0;
}

void FluenceTally::merge(const FluenceTally& other)
{
    if (!axes_.same_axis(other.axes_))
        throw std::invalid_argument(
            "FluenceTally::merge: energy/mu axes differ");
    for (int o = 0; o < static_cast<int>(FluxOrigin::Count); ++o) {
        if (!component_axes_[o].same_axis(other.component_axes_[o]))
            throw std::invalid_argument(
                "FluenceTally::merge: component energy/mu axes differ");
    }
    if (other.lines_.size() != lines_.size())
        throw std::invalid_argument(
            "FluenceTally::merge: line-channel sets differ");
    for (std::size_t l = 0; l < lines_.size(); ++l) {
        if (other.lines_[l].energy_keV != lines_[l].energy_keV ||
            other.lines_[l].origin != lines_[l].origin)
            throw std::invalid_argument(
                "FluenceTally::merge: line-channel definitions differ");
    }
    n_histories_ += other.n_histories_;
    for (std::size_t i = 0; i < sum_w_.size(); ++i) {
        sum_w_[i] += other.sum_w_[i];
        sum_w2_[i] += other.sum_w2_[i];
    }
    for (std::size_t i = 0; i < class_sum_w_.size(); ++i) {
        class_sum_w_[i] += other.class_sum_w_[i];
        class_sum_w2_[i] += other.class_sum_w2_[i];
    }
    for (std::size_t i = 0; i < origin_sum_w_.size(); ++i)
        origin_sum_w_[i] += other.origin_sum_w_[i];
    for (std::size_t i = 0; i < component_sum_w_.size(); ++i) {
        component_sum_w_[i] += other.component_sum_w_[i];
        component_sum_w2_[i] += other.component_sum_w2_[i];
    }
    for (int c = 0; c < kFluxComponentCount; ++c) {
        component_total_w_[c] += other.component_total_w_[c];
        component_total_w2_[c] += other.component_total_w2_[c];
        component_dropped_weight_[c] += other.component_dropped_weight_[c];
    }
    for (int o = 0; o < static_cast<int>(FluxOrigin::Count); ++o) {
        origin_total_w_[o] += other.origin_total_w_[o];
        origin_total_w2_[o] += other.origin_total_w2_[o];
    }
    for (std::size_t i = 0; i < line_sum_w_.size(); ++i) {
        line_sum_w_[i] += other.line_sum_w_[i];
        line_sum_w2_[i] += other.line_sum_w2_[i];
    }
    for (std::size_t i = 0; i < line_class_sum_w_.size(); ++i) {
        line_class_sum_w_[i] += other.line_class_sum_w_[i];
        line_class_sum_w2_[i] += other.line_class_sum_w2_[i];
    }
    for (std::size_t i = 0; i < line_component_sum_w_.size(); ++i) {
        line_component_sum_w_[i] += other.line_component_sum_w_[i];
        line_component_sum_w2_[i] += other.line_component_sum_w2_[i];
    }
    for (std::size_t l = 0; l < line_total_w_.size(); ++l) {
        line_total_w_[l] += other.line_total_w_[l];
        line_total_w2_[l] += other.line_total_w2_[l];
    }
    sum_total_w_ += other.sum_total_w_;
    sum_total_w2_ += other.sum_total_w2_;
    dropped_weight_ += other.dropped_weight_;
    dropped_weight_sum_w2_ += other.dropped_weight_sum_w2_;
    max_history_total_ = std::max(max_history_total_, other.max_history_total_);
    for (int r = 0; r < kNumLogicalReplicas; ++r)
        replica_histories_[r] += other.replica_histories_[r];
    for (int r = 0; r < kNumLogicalReplicas; ++r) {
        const std::vector<float>& other_continuum =
            other.replica_continuum_sum_w_[r];
        if (!other_continuum.empty()) {
            std::vector<float>& continuum = replica_continuum_sum_w_[r];
            if (continuum.empty()) continuum.assign(sum_w_.size(), 0.0f);
            for (std::size_t i = 0; i < continuum.size(); ++i)
                continuum[i] += other_continuum[i];
        }
        const std::vector<float>& other_components =
            other.replica_component_sum_w_[r];
        if (!other_components.empty()) {
            std::vector<float>& components = replica_component_sum_w_[r];
            if (components.empty())
                components.assign(component_sum_w_.size(), 0.0f);
            for (std::size_t i = 0; i < components.size(); ++i)
                components[i] += other_components[i];
        }
        const std::vector<float>& other_line_total =
            other.replica_line_sum_w_[r];
        if (!other_line_total.empty()) {
            std::vector<float>& line_total = replica_line_sum_w_[r];
            if (line_total.empty()) line_total.assign(line_sum_w_.size(), 0.0f);
            for (std::size_t i = 0; i < line_total.size(); ++i)
                line_total[i] += other_line_total[i];
        }
        const std::vector<float>& other_lines =
            other.replica_line_component_sum_w_[r];
        if (!other_lines.empty()) {
            std::vector<float>& lines = replica_line_component_sum_w_[r];
            if (lines.empty())
                lines.assign(line_component_sum_w_.size(), 0.0f);
            for (std::size_t i = 0; i < lines.size(); ++i)
                lines[i] += other_lines[i];
        }
    }
    for (int i = 0; i < kCoarseDimension; ++i)
        coarse_sum_w_[i] += other.coarse_sum_w_[i];
    for (int i = 0; i < kCoarseDimension * kCoarseDimension; ++i)
        coarse_outer_sum_w_[i] += other.coarse_outer_sum_w_[i];
}

double FluenceTally::mean(int ie, int imu) const
{
    if (n_histories_ == 0) return 0.0;
    return sum_w_[axes_.index(ie, imu)] / static_cast<double>(n_histories_);
}

double FluenceTally::sigma(int ie, int imu) const
{
    const int idx = axes_.index(ie, imu);
    return standard_error(n_histories_, sum_w_[idx], sum_w2_[idx]);
}

double FluenceTally::mean_class(int ie, int imu, ChainClass cls) const
{
    if (n_histories_ == 0) return 0.0;
    const int cidx = axes_.index(ie, imu) * static_cast<int>(ChainClass::Count) +
                     static_cast<int>(cls);
    return class_sum_w_[cidx] / static_cast<double>(n_histories_);
}

double FluenceTally::total_mean() const
{
    if (n_histories_ == 0) return 0.0;
    return sum_total_w_ / static_cast<double>(n_histories_);
}

double FluenceTally::total_sigma() const
{
    // Computed from the PER-HISTORY TOTAL, which is exact.
    //
    // Summing the per-bin variances in quadrature would NOT be conservative,
    // even though it looks it (as if a history deposited into one bin RATHER
    // THAN another, a negative multinomial-style covariance). For this estimator:
    // a single history scores at MANY vertices into MANY bins simultaneously,
    // and a history whose walk happens to pass near the detector scores large in
    // all of them at once. The covariance is therefore POSITIVE, and quadrature
    // UNDER-states the uncertainty -- the unsafe direction. A batch-means probe
    // measured the shortfall at a factor 1.485 for a ground geometry with a
    // detector-sized scoring sphere, i.e. exactly the production regime.
    return standard_error(n_histories_, sum_total_w_, sum_total_w2_);
}

double FluenceTally::total_mean_class(ChainClass cls) const
{
    if (n_histories_ == 0) return 0.0;
    double s = 0.0;
    const int nc = static_cast<int>(ChainClass::Count);
    for (std::size_t i = static_cast<int>(cls); i < class_sum_w_.size(); i += nc)
        s += class_sum_w_[i];
    return s / static_cast<double>(n_histories_);
}

double FluenceTally::mean_origin(int ie, int imu, FluxOrigin origin) const
{
    if (n_histories_ == 0) return 0.0;
    const int oidx = axes_.index(ie, imu) * static_cast<int>(FluxOrigin::Count) +
                     static_cast<int>(origin);
    return origin_sum_w_[oidx] / static_cast<double>(n_histories_);
}

double FluenceTally::total_mean_origin(FluxOrigin origin) const
{
    if (n_histories_ == 0) return 0.0;
    return origin_total_w_[static_cast<int>(origin)] /
           static_cast<double>(n_histories_);
}

double FluenceTally::total_sigma_origin(FluxOrigin origin) const
{
    // From the per-history origin total, exact -- same reasoning as
    // total_sigma().
    const int o = static_cast<int>(origin);
    return standard_error(n_histories_, origin_total_w_[o],
                          origin_total_w2_[o]);
}

double FluenceTally::mean_component(int ie, int imu, FluxOrigin origin,
                                    PostCreationOrder order) const
{
    if (n_histories_ == 0) return 0.0;
    const int component = flux_component_index(origin, order);
    const FluenceAxes& ax = component_axes_[static_cast<int>(origin)];
    const std::size_t idx = component_offsets_[component] +
        static_cast<std::size_t>(ax.index(ie, imu));
    return component_sum_w_[idx] / static_cast<double>(n_histories_);
}

double FluenceTally::sigma_component(int ie, int imu, FluxOrigin origin,
                                     PostCreationOrder order) const
{
    const int component = flux_component_index(origin, order);
    const FluenceAxes& ax = component_axes_[static_cast<int>(origin)];
    const std::size_t idx = component_offsets_[component] +
        static_cast<std::size_t>(ax.index(ie, imu));
    return standard_error(n_histories_, component_sum_w_[idx],
                          component_sum_w2_[idx]);
}

double FluenceTally::line_mean_component(int l, int imu, FluxOrigin origin,
                                         PostCreationOrder order) const
{
    if (n_histories_ == 0 || l < 0 || l >= n_lines()) return 0.0;
    const int line_cells = n_lines() * axes_.n_mu;
    const int idx = flux_component_index(origin, order) * line_cells +
                    l * axes_.n_mu + imu;
    return line_component_sum_w_[idx] / static_cast<double>(n_histories_);
}

double FluenceTally::line_sigma_component(int l, int imu, FluxOrigin origin,
                                          PostCreationOrder order) const
{
    if (n_histories_ < 2 || l < 0 || l >= n_lines()) return 0.0;
    const int line_cells = n_lines() * axes_.n_mu;
    const int idx = flux_component_index(origin, order) * line_cells +
                    l * axes_.n_mu + imu;
    return standard_error(n_histories_, line_component_sum_w_[idx],
                          line_component_sum_w2_[idx]);
}

double FluenceTally::total_mean_component(FluxOrigin origin,
                                          PostCreationOrder order) const
{
    if (n_histories_ == 0) return 0.0;
    const int component = flux_component_index(origin, order);
    return component_total_w_[component] / static_cast<double>(n_histories_);
}

double FluenceTally::total_sigma_component(FluxOrigin origin,
                                           PostCreationOrder order) const
{
    const int component = flux_component_index(origin, order);
    return standard_error(n_histories_, component_total_w_[component],
                          component_total_w2_[component]);
}

double FluenceTally::line_mean(int l, int imu) const
{
    if (n_histories_ == 0 || l < 0 || l >= n_lines()) return 0.0;
    return line_sum_w_[l * axes_.n_mu + imu] / static_cast<double>(n_histories_);
}

double FluenceTally::line_sigma(int l, int imu) const
{
    if (n_histories_ < 2 || l < 0 || l >= n_lines()) return 0.0;
    const int idx = l * axes_.n_mu + imu;
    return standard_error(n_histories_, line_sum_w_[idx], line_sum_w2_[idx]);
}

double FluenceTally::line_mean_class(int l, int imu, ChainClass cls) const
{
    if (n_histories_ == 0 || l < 0 || l >= n_lines()) return 0.0;
    const int cidx = (l * axes_.n_mu + imu) * static_cast<int>(ChainClass::Count) +
                     static_cast<int>(cls);
    return line_class_sum_w_[cidx] / static_cast<double>(n_histories_);
}

double FluenceTally::line_total_mean(int l) const
{
    if (n_histories_ == 0 || l < 0 || l >= n_lines()) return 0.0;
    return line_total_w_[l] / static_cast<double>(n_histories_);
}

double FluenceTally::line_total_sigma(int l) const
{
    if (n_histories_ < 2 || l < 0 || l >= n_lines()) return 0.0;
    // From the per-history line total, exact -- same reasoning as total_sigma().
    return standard_error(n_histories_, line_total_w_[l], line_total_w2_[l]);
}

} // namespace ceelo
