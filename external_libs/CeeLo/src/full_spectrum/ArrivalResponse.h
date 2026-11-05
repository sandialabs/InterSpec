#pragma once
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

/// @file ArrivalResponse.h
/// @brief The detector response averaged over a band of arrival directions.
///
/// Scattered photons reach the detector from every side, so folding a
/// scattered field needs the response averaged over arrival direction, not
/// at one beam angle. An arrival band is an isotropic field restricted to
/// direction cosines [mu_lo, mu_hi] (mu = +1 enters the front face head on,
/// -1 the back face). Per envelope photon, its response is the fixed-angle
/// response weighted by the envelope's projected area A(mu): a uniform field
/// puts photons on the envelope in proportion to the area it presents.
///
/// Two sides of the same average live here, so the bake's truth and the
/// model it is scored against cannot drift apart:
///   * ArrivalBeam draws Monte Carlo photons from a band (the detector
///     campaign's `--incidence`);
///   * arrival_quadrature / evaluate_arrival integrate a fitted payload over
///     the band by Gauss-Legendre in the polar angle, per hemisphere.
///
/// ArrivalResponseTable tabulates the same payload on (energy, mu) nodes so
/// that an arbitrary incident field -- a scattered continuum with its own
/// angular distribution -- folds by one contraction, Jet-safe in the field.

#include "full_spectrum/FullSpectrumModel.h"

#include <Eigen/Core>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace ceelo {
class Geometry;
}

namespace ceelo::full_spectrum {

/// Direction cosines [mu_lo, mu_hi], filled uniformly in mu and azimuth.
struct ArrivalBand {
    double mu_lo = -1.0;
    double mu_hi = 1.0;
};

/// The named arrival classes: "iso" [-1, 1], "front" [0, 1], "back" [-1, 0],
/// "side" [-0.3, 0.3], or an explicit "mu-band:LO,HI" with
/// -1 <= LO < HI <= 1. Throws std::invalid_argument for anything else.
ArrivalBand arrival_band(const std::string& name);

/// Projected area of a cylindrical envelope seen from direction cosine mu:
/// pi R^2 |mu| + 2 R L sqrt(1 - mu^2).
double envelope_projected_area_cm2(double envelope_radius_cm,
                                   double envelope_length_cm, double mu);

/// Nodes in mu and weights summing to one for averaging a fixed-angle
/// response over a band, the projected area folded into the weights.
/// Gauss-Legendre in the polar angle, `nodes_per_hemisphere` points on each
/// side of mu = 0 the band covers (the area has a kink there, which a rule
/// spanning it resolves slowly). `mean_projected_area_cm2` is the band's
/// average of A(mu) over uniform mu: photons per unit fluence per envelope.
struct ArrivalQuadrature {
    std::vector<double> mu;
    std::vector<double> weight;
    double mean_projected_area_cm2 = 0.0;
};

ArrivalQuadrature arrival_quadrature(const ArrivalBand& band,
                                     double envelope_radius_cm,
                                     double envelope_length_cm,
                                     std::size_t nodes_per_hemisphere = 16);

/// The payload's response averaged over `band`, per envelope photon, on
/// `deposit_bins_keV`: atoms merged by (component, energy), continuum by
/// (component, bin). Throws std::invalid_argument if the payload carries no
/// envelope (a hand-set pilot).
SpectralMeasure<double> evaluate_arrival(
    const DetectorPulseResponseModel& model, double energy_keV,
    const ArrivalBand& band, const std::vector<double>& deposit_bins_keV,
    std::size_t nodes_per_hemisphere = 16);

/// Monte Carlo photons arriving from a band. Each draw takes mu uniformly in
/// the band, azimuth uniformly, and an origin uniformly on the disk of the
/// geometry's bounding sphere seen from that direction; direction AND origin
/// are redrawn until the ray hits the outer boundary. Accepted photons are
/// therefore distributed in mu as A(mu), and a tally divided by the accepted
/// count is per envelope photon, as for a fixed-angle beam.
class ArrivalBeam {
public:
    ArrivalBeam(const Geometry& geometry, const ArrivalBand& band);

    void sample(std::mt19937_64& rng, Eigen::Vector3d& origin,
                Eigen::Vector3d& direction);

    std::uint64_t accepted() const { return accepted_; }
    std::uint64_t tries() const { return tries_; }

    /// The band's mean projected area of the outer boundary, estimated as
    /// the bounding disk's area times the acceptance, and its binomial sigma.
    double mean_projected_area_cm2() const;
    double mean_projected_area_sigma_cm2() const;

private:
    const Geometry& geometry_;
    ArrivalBand band_;
    Eigen::Vector3d center_;
    double bound_cm_ = 0.0;
    std::uint64_t accepted_ = 0;
    std::uint64_t tries_ = 0;
};

/// An incident photon field at the detector, as fluence (photons per cm², or
/// per cm² per source photon): a continuum binned in energy and arrival
/// cosine, and lines with an angular distribution. Within a bin the photons
/// are uniform in energy and in mu (mu = +1 is head-on into the front face).
/// Only the fluence values carry Scalar, so a fit's parameters reach a fold
/// through them alone.
template <typename Scalar>
struct IncidentField {
    std::vector<double> mu_edges{-1.0, 1.0};  ///< ascending, within [-1, 1]
    std::vector<double> energy_edges_keV;     ///< continuum cells, ascending
    std::vector<Scalar> continuum;            ///< [energy cell][mu bin]
    struct Line {
        double energy_keV = 0.0;
        std::vector<Scalar> fluence;          ///< per mu bin
    };
    std::vector<Line> lines;
};

struct ArrivalTableOptions {
    /// Incident-energy nodes, log-spaced with ratio at most 1 + energy_step;
    /// a field's continuum and lines must lie within [min, max].
    double energy_min_keV = 10.0;
    double energy_max_keV = 3000.0;
    double energy_step = 0.02;
    std::size_t nodes_per_hemisphere = 16;
    /// Line energies tabulated exactly (the direct lines a composition folds
    /// again and again); other lines interpolate between energy nodes.
    std::vector<double> exact_line_energies_keV;
    /// Build threads; 0 = the hardware's.
    unsigned threads = 0;
};

/// A detector payload's response tabulated for folding incident fields: per
/// arrival-cosine node (Gauss-Legendre in theta per hemisphere) and incident
/// energy node, the continuum and fixed-energy atoms on the caller's deposit
/// grid, and the yields of the atoms that track the incident energy (full
/// energy, pair full energy, single and double escape, K escape).
///
/// A fold maps each field cell onto the nodes -- a hat basis in energy (the
/// cell's mass is spread uniformly over it), and in mu the response as the
/// polynomial in theta through each hemisphere's nodes, integrated against the
/// envelope's exact projected area -- and contracts with the table;
/// the tracking atoms are transferred exactly instead, each cell's atom mass
/// landing on its own deposit interval [E_lo - offset, E_hi - offset], so a
/// smooth incident continuum never renders as a comb of node lines. Counts
/// per deposit bin, unbroadened, zero deposits left out.
class ArrivalResponseTable {
public:
    ArrivalResponseTable(const DetectorPulseResponseModel& model,
                         std::vector<double> deposit_edges_keV,
                         ArrivalTableOptions options = {});

    const std::vector<double>& deposit_edges_keV() const
    {
        return deposit_edges_;
    }
    const std::vector<double>& energy_nodes_keV() const { return energies_; }
    const std::vector<double>& mu_nodes() const { return mu_; }
    std::size_t table_bytes() const;

    /// Expected counts per deposit bin for `field`. Throws
    /// std::invalid_argument for a malformed field or one outside the table's
    /// energy range.
    template <typename Scalar>
    std::vector<Scalar> fold(const IncidentField<Scalar>& field) const;

private:
    struct Weight {
        std::size_t index;
        double weight;
    };
    struct Tracked {
        SpectralComponent component;
        double offset_keV;            ///< incident energy - atom energy
        std::vector<double> yield;    ///< [mu node * n_energy + energy node]
    };
    struct ExactLine {
        double energy_keV;
        std::size_t length;           ///< deposit bins below the line
        std::vector<float> rows;      ///< [mu node][length]
        std::vector<double> yield;    ///< [tracked * n_mu + mu node]
    };

    std::vector<Weight> energy_cell_weights(double lo, double hi) const;
    std::vector<Weight> energy_point_weights(double energy) const;
    std::vector<Weight> mu_bin_weights(double lo, double hi) const;
    int exact_line(double energy) const;
    template <typename Scalar>
    void deposit(std::vector<Scalar>& out, double lo, double hi,
                 const Scalar& mass) const;

    std::vector<double> deposit_edges_;
    std::vector<double> energies_;
    std::vector<double> mu_;
    double envelope_radius_cm_ = 0.0;
    double envelope_length_cm_ = 0.0;
    /// Each hemisphere's mu nodes in ascending order, and their indices.
    std::vector<double> front_mu_, back_mu_;
    std::vector<std::size_t> front_index_, back_index_;
    std::vector<std::size_t> length_;     ///< per energy node
    std::vector<std::size_t> offset_;     ///< per (mu node, energy node)
    std::vector<float> rows_;             ///< per envelope photon
    std::vector<Tracked> tracked_;
    std::vector<ExactLine> lines_;
};

template <typename Scalar>
void ArrivalResponseTable::deposit(std::vector<Scalar>& out, double lo,
                                   double hi, const Scalar& mass) const
{
    const auto& e = deposit_edges_;
    if (!(hi > lo)) {
        // A point: the bin holding it, if any.
        if (lo < e.front() || lo >= e.back())
            return;
        const auto bin = static_cast<std::size_t>(
            std::upper_bound(e.begin(), e.end(), lo) - e.begin() - 1);
        out[bin] += mass;
        return;
    }
    const double width = hi - lo;
    auto first = std::upper_bound(e.begin(), e.end(), lo);
    std::size_t bin = first == e.begin()
        ? 0 : static_cast<std::size_t>(first - e.begin() - 1);
    for (; bin + 1 < e.size() && e[bin] < hi; ++bin) {
        const double overlap = std::min(hi, e[bin + 1]) - std::max(lo, e[bin]);
        if (overlap > 0.0)
            out[bin] += mass * (overlap / width);
    }
}

template <typename Scalar>
std::vector<Scalar> ArrivalResponseTable::fold(
    const IncidentField<Scalar>& field) const
{
    const std::size_t n_bins = field.mu_edges.size() >= 2
        ? field.mu_edges.size() - 1 : 0;
    if (n_bins == 0 || field.mu_edges.front() < -1.0
        || field.mu_edges.back() > 1.0)
        throw std::invalid_argument("fold: mu edges must lie in [-1, 1]");
    for (std::size_t m = 0; m < n_bins; ++m) {
        if (!(field.mu_edges[m + 1] > field.mu_edges[m]))
            throw std::invalid_argument("fold: mu edges must ascend");
    }
    const std::size_t n_cells = field.energy_edges_keV.size() >= 2
        ? field.energy_edges_keV.size() - 1 : 0;
    if (field.continuum.size() != n_cells * n_bins)
        throw std::invalid_argument(
            "fold: continuum must hold energy cells x mu bins values");

    const std::size_t n_mu = mu_.size();
    const std::size_t n_energy = energies_.size();
    std::vector<Scalar> out(deposit_edges_.size() - 1, Scalar(0));
    std::vector<std::vector<Weight>> mu_weights(n_bins);
    for (std::size_t m = 0; m < n_bins; ++m)
        mu_weights[m] = mu_bin_weights(field.mu_edges[m], field.mu_edges[m + 1]);

    // Node weights of the whole field, and which nodes it touches.
    std::vector<Scalar> node(n_mu * n_energy, Scalar(0));
    std::vector<char> touched(n_mu * n_energy, 0);
    std::vector<double> coefficient(tracked_.size() * n_bins);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const double lo = field.energy_edges_keV[c];
        const double hi = field.energy_edges_keV[c + 1];
        if (!(hi > lo))
            throw std::invalid_argument("fold: energy edges must ascend");
        const auto energy_weights = energy_cell_weights(lo, hi);
        std::fill(coefficient.begin(), coefficient.end(), 0.0);
        for (std::size_t m = 0; m < n_bins; ++m) {
            const Scalar& fluence = field.continuum[c * n_bins + m];
            for (const Weight& mw : mu_weights[m]) {
                for (const Weight& ew : energy_weights) {
                    const std::size_t index = mw.index * n_energy + ew.index;
                    const double w = mw.weight * ew.weight;
                    node[index] += fluence * w;
                    touched[index] = 1;
                    for (std::size_t t = 0; t < tracked_.size(); ++t)
                        coefficient[t * n_bins + m] += w * tracked_[t].yield[index];
                }
            }
        }
        for (std::size_t t = 0; t < tracked_.size(); ++t) {
            Scalar mass = Scalar(0);
            bool any = false;
            for (std::size_t m = 0; m < n_bins; ++m) {
                const double k = coefficient[t * n_bins + m];
                if (k != 0.0) {
                    mass += field.continuum[c * n_bins + m] * k;
                    any = true;
                }
            }
            if (any)
                deposit(out, lo - tracked_[t].offset_keV,
                        hi - tracked_[t].offset_keV, mass);
        }
    }

    std::vector<Scalar> line_node(lines_.size() * n_mu, Scalar(0));
    for (const auto& line : field.lines) {
        if (line.fluence.size() != n_bins)
            throw std::invalid_argument("fold: a line needs one value per mu bin");
        const int exact = exact_line(line.energy_keV);
        const auto energy_weights = exact >= 0
            ? std::vector<Weight>{} : energy_point_weights(line.energy_keV);
        std::fill(coefficient.begin(), coefficient.end(), 0.0);
        for (std::size_t m = 0; m < n_bins; ++m) {
            for (const Weight& mw : mu_weights[m]) {
                if (exact >= 0) {
                    line_node[static_cast<std::size_t>(exact) * n_mu + mw.index] +=
                        line.fluence[m] * mw.weight;
                    for (std::size_t t = 0; t < tracked_.size(); ++t)
                        coefficient[t * n_bins + m] += mw.weight
                            * lines_[static_cast<std::size_t>(exact)]
                                  .yield[t * n_mu + mw.index];
                    continue;
                }
                for (const Weight& ew : energy_weights) {
                    const std::size_t index = mw.index * n_energy + ew.index;
                    const double w = mw.weight * ew.weight;
                    node[index] += line.fluence[m] * w;
                    touched[index] = 1;
                    for (std::size_t t = 0; t < tracked_.size(); ++t)
                        coefficient[t * n_bins + m] += w * tracked_[t].yield[index];
                }
            }
        }
        for (std::size_t t = 0; t < tracked_.size(); ++t) {
            Scalar mass = Scalar(0);
            bool any = false;
            for (std::size_t m = 0; m < n_bins; ++m) {
                const double k = coefficient[t * n_bins + m];
                if (k != 0.0) {
                    mass += line.fluence[m] * k;
                    any = true;
                }
            }
            if (any) {
                const double at = line.energy_keV - tracked_[t].offset_keV;
                deposit(out, at, at, mass);
            }
        }
    }

    for (std::size_t k = 0; k < n_mu; ++k) {
        for (std::size_t j = 0; j < n_energy; ++j) {
            const std::size_t index = k * n_energy + j;
            if (!touched[index])
                continue;
            const float* row = rows_.data() + offset_[index];
            const Scalar& w = node[index];
            for (std::size_t d = 0; d < length_[j]; ++d)
                out[d] += w * static_cast<double>(row[d]);
        }
    }
    for (std::size_t l = 0; l < lines_.size(); ++l) {
        for (std::size_t k = 0; k < n_mu; ++k) {
            const Scalar& w = line_node[l * n_mu + k];
            const float* row = lines_[l].rows.data() + k * lines_[l].length;
            for (std::size_t d = 0; d < lines_[l].length; ++d)
                out[d] += w * static_cast<double>(row[d]);
        }
    }
    return out;
}

} // namespace ceelo::full_spectrum
