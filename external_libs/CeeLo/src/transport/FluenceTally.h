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

/// @file FluenceTally.h
/// @brief Binned incident-fluence tally, (E', mu_in), with correlated-history
///        variance accounting and scatter-chain class labels.
///
/// THE BINDING CONSTRAINT. This is a LINEAR fluence tally, not a detector
/// estimator. The peak transmission T(E0) and the scattered continuum are the
/// mutually exclusive outcomes of one primary photon
/// (T + sum(scatter fractions) + absorbed = 1), so they must never be summed
/// into one history and passed through a non-linear indicator. CeeLo's existing
/// event tallies do exactly that (all branch deposits of an event are summed and
/// the FEP/total indicators applied once), which is correct for THEM and wrong
/// here. Detector response is applied downstream, separately, as the linear
/// combination
///     result = FEP(E0)*T + integral FEP(E', mu_in) dphi_scat
/// Folding a detector efficiency into this tally would reintroduce exactly the
/// non-linearity hazard CLAUDE.md flags under "Branch-level Russian roulette is
/// fundamentally unsafe for FEP".
///
/// VARIANCE. A single primary history scores at every vertex of its scatter
/// walk, so the scores within one history are CORRELATED. Accumulating them
/// directly into sum_w2 would understate the variance. Scores are therefore
/// buffered per history and folded into the sum/sum-of-squares once, when
/// `end_history()` is called.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ceelo {

/// How a scored chain reached the detector. Used to decompose a full-room run
/// into the pieces the kernel superposition model assumes, and to MEASURE the
/// inter-surface bounce term that the model neglects rather than assuming it.
/// Chain classes are keyed on the number of SLAB INTERACTIONS (bounces), not
/// merely on how many distinct surfaces were touched, because the room
/// decomposition needs both concepts and they are not the same:
///
///   * superposing per-surface kernels can only reproduce chains with EXACTLY
///     ONE bounce -- an isolated single-surface bake also contains chains that
///     bounce repeatedly off that same surface, and summing those over six
///     surfaces double counts (measured: +47% on the one-surface term, and a
///     +31% overshoot once a bath was added on top);
///   * the diffuse-bath term must therefore carry ALL >= 2-bounce chains,
///     same-surface and cross-surface alike.
///
/// `MultiDistinct` is kept separate because it is the piece that is genuinely
/// impossible in any single-surface geometry, and so is the honest measure of
/// what a superposition model omits.
enum class ChainClass : uint8_t {
    AirOnly = 0,        ///< No slab interaction at all.
    OneBounce = 1,      ///< Exactly one slab interaction -- the superposable part.
    MultiSameSlab = 2,  ///< Two or more bounces, all on the SAME slab.
    MultiDistinct = 3,  ///< Two or more bounces spanning DIFFERENT slabs.
    Count = 4
};

/// Where a discrete LINE in the incident spectrum comes from.
///
/// Some scattered flux arrives at an EXACT energy rather than spread over a
/// continuum, and binning it into a log energy bin destroys that distinction.
/// It matters downstream: a line makes a photopeak in the detector, continuum in
/// the same bin does not. The kernel already keeps the primary's uncollided term
/// out of the continuum for exactly this reason; these channels extend the same
/// treatment to the two lines the environment itself produces.
enum class LineOrigin : uint8_t {
    /// Chains that underwent ONLY coherent (Rayleigh) scatters. Rayleigh is
    /// elastic, so these arrive at exactly E0. Measured share of the scattered
    /// flux in a soil field: 3.55% at
    /// 60 keV -- where it is HALF the content of its energy bin -- 0.47% at
    /// 122 keV, and negligible above. Worth its own channel because 60 keV is
    /// squarely in the regime the low-energy work targets.
    ElasticE0 = 0,
    /// The 511 keV photons from pair production, both their uncollided arrival
    /// and any subsequent all-elastic chain.
    Annihilation = 1
};

struct LineDef {
    double energy_keV = 0.0;
    LineOrigin origin = LineOrigin::ElasticE0;
};

/// How the scored photon was CREATED, mirroring the GEANT4 harness's
/// creator-process ORIGIN DECOMPOSITION so the two decompose the scattered
/// flux identically: the tag propagates through Compton/Rayleigh scatters
/// (which keep the same track) and is OVERWRITTEN at each secondary-photon
/// creation -- a 511 from a pair vertex on a brems chain is Annihilation, a
/// brems photon from an annihilation-chain electron is Brems. Chosen over
/// with/without-channel differencing for the same reason the G4 side
/// abandoned differencing: a few-percent channel differenced between two
/// independent runs carries the sum of both variances and every bin is noise,
/// while the tag is exact in one run.
enum class FluxOrigin : uint8_t {
    Primary = 0,       ///< The primary photon and its own scatter chain.
    Annihilation = 1,  ///< 511s from pair production, and their chains.
    Brems = 2,         ///< Bremsstrahlung photons, and their chains.
    Count = 3
};

/// Number of interactions undergone by a photon since its creation.
///
/// This is deliberately local to the current photon, rather than the whole
/// primary history.  Creating an annihilation or bremsstrahlung photon resets
/// the order to Zero; its first Compton/Rayleigh product is One.  TwoPlus is a
/// saturated bucket because the fit needs the exact single-scatter term while
/// retaining all higher-order transport as one additive remainder.
enum class PostCreationOrder : uint8_t {
    Zero = 0,
    One = 1,
    TwoPlus = 2,
    Count = 3
};

constexpr int kFluxComponentCount =
    static_cast<int>(FluxOrigin::Count) *
    static_cast<int>(PostCreationOrder::Count);

/// Component-major index used by all raw fit arrays.
constexpr int flux_component_index(FluxOrigin origin,
                                   PostCreationOrder order)
{
    return static_cast<int>(origin) *
               static_cast<int>(PostCreationOrder::Count) +
           static_cast<int>(order);
}

constexpr FluxOrigin flux_component_origin(int component)
{
    return static_cast<FluxOrigin>(
        component / static_cast<int>(PostCreationOrder::Count));
}

constexpr PostCreationOrder flux_component_order(int component)
{
    return static_cast<PostCreationOrder>(
        component % static_cast<int>(PostCreationOrder::Count));
}

/// Axis definition for the tally: log-spaced energy bins, uniform mu bins.
struct FluenceAxes {
    /// The walk's `e_floor_keV` matches this, so sub-30 keV flux is not merely
    /// untallied, it is not transported. That is policy, not a data limit (the
    /// cross-section tables start at 10 keV, the angular norms at 1 keV), and
    /// lowering it is provably additive: Compton strictly loses energy and
    /// Rayleigh preserves it, so nothing can re-enter the >= 30 keV band from
    /// below. Measured at 20 keV: the
    /// 20-30 keV band is <= 0.02% of the scattered flux even for an Am-241
    /// source, for no measurable run time. Free to move, but it buys nothing
    /// until environmental fluorescence is emitted or a planar/LEGe detector is
    /// in scope.
    double e_min_keV = 30.0;
    double e_max_keV = 3000.0;
    /// 48 is enough for INTEGRATED counts to 0.005% (energy_binning_study.cpp),
    /// and finer or physics-placed grids buy nothing there. It is NOT enough to
    /// carry spectral SHAPE at HPGe resolution. Folding a real pulse-height
    /// response (hpge_binning_study.cpp) leaves the integral exact -- the error
    /// is purely in where the counts land. Percent of counts in the wrong
    /// channel, reading the table honestly as a histogram:
    ///
    ///     bins     59.5 keV   122 keV   662 keV   2614 keV
    ///       48      17.31      9.70      2.17      1.62
    ///       96      10.94      3.35      1.47      1.13
    ///      192       4.22      1.88      0.95      0.84
    ///      384       2.05      1.00      0.57      0.64
    ///      768       0.52      0.25      0.29      0.40
    ///
    /// (reference MC noise floor 0.34-0.71%). Low source energies are worst:
    /// the whole albedo spectrum is squeezed into a few bins AND the detector
    /// does almost no smearing of its own there, so the incident grid IS the
    /// output resolution. ~768-bin density reaches the noise floor everywhere.
    ///
    /// CONSUMERS MUST NOT treat a bin as a delta at its midpoint energy. Doing
    /// so leaves the integral exact while putting 58-74% of the counts in the
    /// wrong channel -- a comb of fabricated lines. Spread the bin's content
    /// across its width.
    int n_energy = 48;
    /// Uniform over mu_in in [-1, 1]. 32 was chosen by measurement, not taste:
    /// the mu axis exists so the next phase can fold a detector response
    /// eps(E,mu), so the bin count was set by the discretization error in
    /// counts = sum_mu phi(E,mu) A_proj(mu) eps(E,mu). Measured worst case over
    /// a soil field and a concrete room, for both a 3"x3" NaI and an extreme
    /// thin-slab aspect ratio:
    ///     n_mu   8     12     16     24     32     64
    ///     err  1.95%  1.20%  0.56%  0.26%  0.10%  0.13%
    /// Extra mu bins cost NO additional MC time (same histories, finer binning),
    /// only table size, so 32 buys ~20x accuracy over 8 for 4x storage. Drop to
    /// 16 (0.56%) only if table size becomes the binding constraint.
    int n_mu = 32;

    /// OPTIONAL explicit energy bin edges, ascending, size n_energy + 1.
    ///
    /// When empty (the default) the axis is log-uniform over
    /// [e_min_keV, e_max_keV], which is what the GEANT4 comparison grid
    /// uses. When populated it defines the bins directly, so
    /// block boundaries can be PINNED at physical landmarks (see
    /// make_landmark_axes). The distinction matters because the dominant
    /// binning error is not resolution but STRADDLE: a bin that spans a
    /// spectral edge smears counts across it no matter how the consumer
    /// re-spreads the bin, and the worst measured case (9.07% misplaced at
    /// 122 keV) was exactly a backscatter edge sitting mid-bin.
    std::vector<double> edges;

    int energy_bin(double e_keV) const;
    int mu_bin(double mu) const;
    double energy_low(int i) const;
    double energy_high(int i) const;
    double mu_low(int j) const;
    double mu_high(int j) const;
    int n_bins() const { return n_energy * n_mu; }
    int index(int ie, int imu) const { return ie * n_mu + imu; }

    /// True when this axis carries explicit edges rather than a log-uniform rule.
    bool has_edges() const { return !edges.empty(); }

    /// Structural equality (an axis mismatch between two tallies is a bug, not
    /// a mergeable difference).
    bool same_axis(const FluenceAxes& o) const;
};

/// Build a landmark-blocked energy axis for a node of primary energy `e0_keV`.
///
/// The scattered spectrum has three physically fixed features, and a bin that
/// straddles any of them destroys information the downstream fit needs:
///   * the **Compton backscatter edge** E0/(1 + 2 E0/511) -- the sharpest
///     feature in an albedo spectrum, and the reason 122 keV was the worst
///     binning case (its edge is at 82.6 keV, mid-bin on a log-uniform grid);
///   * **511 keV**, where the annihilation continuum from pair-produced
///     positrons begins (only for E0 above threshold);
///   * **E0 itself**, the hard upper limit: nothing scatters up.
/// Edges are pinned at each of those, the axis stops at E0 (so above-E0 bins
/// never exist rather than being stored as structural zeros), and bins are
/// distributed log-uniformly WITHIN each block.
///
/// `n_energy_target` is the total bin count; blocks receive bins in proportion
/// to their log width, with `albedo_boost` extra weight on the block below the
/// backscatter edge (the 25-255 keV albedo decade carries the spectral peak
/// and is where measured misplacement was worst). Blocks always get >= 2 bins.
FluenceAxes make_landmark_axes(double e0_keV, int n_energy_target = 192,
                               int n_mu = 32, double e_min_keV = 30.0,
                               double albedo_boost = 2.0);

/// Landmark axis with REGISTERED bins: two blocks ([e_min, E_bs] and
/// [E_bs, E0]) each carrying the SAME bin count for every E0.
///
/// This is the axis to use when E0 is an interpolation axis. `make_landmark_axes`
/// pins the block EDGES at physical landmarks but allocates bins between blocks
/// in proportion to their log-widths, which vary with E0 -- so bin k lands at a
/// different place in the spectrum for different E0 (at 192 bins: 180/12 split
/// at 59 keV versus 139/53 at 2614 keV). Interpolating bin k against bin k then
/// silently mixes unlike parts of the spectrum, which is exactly the
/// feature-misregistration the landmark idea exists to prevent.
///
/// With a fixed split, bin k sits at the same fractional position
///     u = (k + 1/2) / n_energy
/// of the same block for EVERY node, so the backscatter edge is always the
/// block boundary at u = 0.5 and E0 is always the top edge. The scaled
/// coordinate is then built into the STORAGE: interpolating across E0 is a
/// plain bin-index-to-bin-index operation with features automatically
/// registered, and no separate coordinate transform is needed at query time.
///
/// 511 keV is deliberately NOT a block boundary here: the annihilation
/// continuum exists only above the 1022 keV pair threshold, so a 511 block
/// would change the number of blocks partway along the E0 axis and destroy the
/// very registration this function provides. The annihilation component is
/// handled separately, by origin (its shape is E0-independent where it matters).
///
/// `n_energy_target` is rounded up to an even number so the two blocks match.
FluenceAxes make_registered_landmark_axes(double e0_keV,
                                          int n_energy_target = 192,
                                          int n_mu = 32,
                                          double e_min_keV = 10.0);

/// Native energy axes for the three fit origins, ordered by FluxOrigin:
/// Primary uses the E0/backscatter registered coordinate, Annihilation uses a
/// fixed 511-keV/backscatter registered coordinate, and Brems uses a log-relative
/// coordinate ending at E0. The common total tally should normally use the
/// Primary axis so total-shape validation never requires component rebinning.
std::vector<FluenceAxes> make_fit_component_axes(
    double e0_keV, int n_energy_target = 384, int n_mu = 64,
    double e_min_keV = 10.0);

/// Binned tally of incident fluence per source photon.
class FluenceTally {
public:
    static constexpr int kNumLogicalReplicas = 32;
    static constexpr int kNumEnergyBands = 12;
    static constexpr int kCoarseDimension = 16;

    explicit FluenceTally(const FluenceAxes& axes = FluenceAxes{},
                          std::vector<LineDef> lines = {},
                          std::vector<FluenceAxes> component_axes = {});

    /// Score a contribution at incident energy `e_keV` and incidence cosine
    /// `mu_in`. Out-of-range energies are dropped (and counted in
    /// `dropped_weight()` so truncation is never silent).
    void add(double e_keV, double mu_in, double weight, ChainClass cls,
             FluxOrigin origin = FluxOrigin::Primary,
             PostCreationOrder order = PostCreationOrder::One);

    /// Score into a discrete LINE channel rather than the continuum grid. The
    /// line's energy is implied by the channel, so only the incidence angle is
    /// binned. `line` indexes `lines()`.
    void add_line(int line, double mu_in, double weight, ChainClass cls,
                  FluxOrigin origin = FluxOrigin::Primary,
                  PostCreationOrder order = PostCreationOrder::One);

    /// Close the current primary history: fold the buffered per-history sums
    /// into the running mean/variance accumulators. Must be called once per
    /// primary, INCLUDING primaries that scored nothing (they are real zeros and
    /// belong in the denominator).
    void end_history();

    /// Close the current history into an explicitly selected logical replica.
    /// This is the entry point used by fit blocks, whose independently seeded
    /// streams have stable replica IDs. `end_history()` assigns histories
    /// round-robin to the same 32 replicas.
    void end_history(int logical_replica);

    /// Merge another tally (per-thread accumulation).
    void merge(const FluenceTally& other);

    const FluenceAxes& axes() const { return axes_; }
    /// Per-origin fit axis. With the default constructor every origin uses
    /// `axes()`. A fit block may instead supply exactly FluxOrigin::Count
    /// axes (Primary registered, Annihilation fixed-511, Brems log-relative).
    /// All component axes must share the common tally's mu binning.
    const FluenceAxes& component_axes(FluxOrigin origin) const {
        return component_axes_[static_cast<int>(origin)];
    }
    uint64_t num_histories() const { return n_histories_; }

    int n_lines() const { return static_cast<int>(lines_.size()); }
    const std::vector<LineDef>& lines() const { return lines_; }
    const LineDef& line(int l) const { return lines_.at(l); }

    /// Mean line fluence per source photon in incidence bin `imu`.
    double line_mean(int l, int imu) const;
    double line_sigma(int l, int imu) const;
    double line_mean_class(int l, int imu, ChainClass cls) const;
    /// Integrated over incidence angle.
    double line_total_mean(int l) const;
    double line_total_sigma(int l) const;

    /// Mean scored fluence per source photon in bin (ie, imu), all classes.
    double mean(int ie, int imu) const;
    /// 1-sigma statistical uncertainty on `mean`.
    double sigma(int ie, int imu) const;
    /// Mean restricted to one chain class.
    double mean_class(int ie, int imu, ChainClass cls) const;

    /// Integral over all bins (mean total scored fluence per source photon).
    /// INCLUDES the line channels, so "total scattered fluence" is independent
    /// of how the flux is split between lines and continuum -- every existing
    /// comparison against it is unaffected by the split.
    double total_mean() const;
    double total_sigma() const;
    /// Integral restricted to one chain class.
    double total_mean_class(ChainClass cls) const;

    /// Continuum-grid mean restricted to one creation origin (mean-only,
    /// like the chain-class decomposition).
    double mean_origin(int ie, int imu, FluxOrigin origin) const;
    /// Total (grid + lines) restricted to one origin, with an exact
    /// correlated-history sigma -- the numbers compared against the G4
    /// harness's ORIGIN DECOMPOSITION block. The origins PARTITION the
    /// total: sum over origins == total_mean() exactly.
    double total_mean_origin(FluxOrigin origin) const;
    double total_sigma_origin(FluxOrigin origin) const;

    /// Fit component accessors. Components are the Cartesian product of
    /// creation origin and post-creation interaction order and exactly
    /// partition every accepted continuum/line score.
    double mean_component(int ie, int imu, FluxOrigin origin,
                          PostCreationOrder order) const;
    double sigma_component(int ie, int imu, FluxOrigin origin,
                           PostCreationOrder order) const;
    double line_mean_component(int l, int imu, FluxOrigin origin,
                               PostCreationOrder order) const;
    double line_sigma_component(int l, int imu, FluxOrigin origin,
                                PostCreationOrder order) const;
    double total_mean_component(FluxOrigin origin,
                                PostCreationOrder order) const;
    double total_sigma_component(FluxOrigin origin,
                                 PostCreationOrder order) const;

    /// Raw mergeable sufficient statistics for the fit archive. Arrays are
    /// deliberately exposed read-only; their layouts are stable and documented
    /// here so serialization does not reconstruct them from rounded means.
    ///
    /// Continuum totals: [energy][mu].
    const std::vector<double>& continuum_sum_w() const { return sum_w_; }
    const std::vector<double>& continuum_sum_w2() const { return sum_w2_; }
    /// Slab-chain diagnostics: [energy][mu][ChainClass].
    const std::vector<double>& class_sum_w() const { return class_sum_w_; }
    const std::vector<double>& class_sum_w2() const { return class_sum_w2_; }
    /// Continuum components are concatenated component-major. A component's
    /// [energy][mu] slice is
    /// [component_offsets()[c], component_offsets()[c+1]); its shape comes
    /// from component_axes(flux_component_origin(c)). This permits native
    /// Primary/Annihilation/Brems energy coordinates while the total remains on
    /// the common `axes()` grid.
    const std::vector<double>& component_sum_w() const {
        return component_sum_w_;
    }
    const std::vector<double>& component_sum_w2() const {
        return component_sum_w2_;
    }
    const std::array<std::size_t, kFluxComponentCount + 1>&
    component_offsets() const { return component_offsets_; }
    /// Line totals: [line][mu].
    const std::vector<double>& line_sum_w() const { return line_sum_w_; }
    const std::vector<double>& line_sum_w2() const { return line_sum_w2_; }
    /// Exact-line slab-chain diagnostics: [line][mu][ChainClass].
    const std::vector<double>& line_class_sum_w() const {
        return line_class_sum_w_;
    }
    const std::vector<double>& line_class_sum_w2() const {
        return line_class_sum_w2_;
    }
    /// Line components: [component][line][mu].
    const std::vector<double>& line_component_sum_w() const {
        return line_component_sum_w_;
    }
    const std::vector<double>& line_component_sum_w2() const {
        return line_component_sum_w2_;
    }
    /// Integrated raw moments, including continuum and exact lines.
    double total_sum_w() const { return sum_total_w_; }
    double total_sum_w2() const { return sum_total_w2_; }
    const std::vector<double>& component_total_sum_w() const {
        return component_total_w_;
    }
    const std::vector<double>& component_total_sum_w2() const {
        return component_total_w2_;
    }
    const std::vector<double>& line_total_sum_w() const {
        return line_total_w_;
    }
    const std::vector<double>& line_total_sum_w2() const {
        return line_total_w2_;
    }

    /// Replica sums are stored as float32 by design: they estimate covariance,
    /// while the authoritative S1/S2 arrays above remain float64. Layouts are
    /// [replica][concatenated component cells] (using component_offsets) and
    /// [replica][component][line][mu].
    const std::array<uint64_t, kNumLogicalReplicas>& replica_histories() const {
        return replica_histories_;
    }
    /// Common-grid total continuum maps: [replica][energy][mu].
    const std::array<std::vector<float>, kNumLogicalReplicas>&
    replica_continuum_sum_w() const { return replica_continuum_sum_w_; }
    /// A replica vector is allocated lazily on its first nonzero score; an
    /// empty vector means an all-zero map for that replica. This avoids an
    /// accidental 32x memory multiplier when a block runs one tally per
    /// independently seeded replica.
    const std::array<std::vector<float>, kNumLogicalReplicas>&
    replica_component_sum_w() const {
        return replica_component_sum_w_;
    }
    /// Total exact-line maps: [replica][line][mu].
    const std::array<std::vector<float>, kNumLogicalReplicas>&
    replica_line_sum_w() const { return replica_line_sum_w_; }
    const std::array<std::vector<float>, kNumLogicalReplicas>&
    replica_line_component_sum_w() const {
        return replica_line_component_sum_w_;
    }

    /// Exact moments of the per-history 16-vector used for stopping: entries
    /// 0..11 are registered-energy band masses, entries 12..15 are the total
    /// fluence weighted by Legendre P1..P4(mu). The outer product is row-major.
    const std::array<double, kCoarseDimension>& coarse_sum_w() const {
        return coarse_sum_w_;
    }
    const std::array<double, kCoarseDimension * kCoarseDimension>&
    coarse_outer_sum_w() const { return coarse_outer_sum_w_; }

    /// Weight that fell outside the energy axis (never silently discarded).
    double dropped_weight() const { return dropped_weight_; }
    double dropped_weight_sum_w2() const { return dropped_weight_sum_w2_; }
    double max_history_total() const { return max_history_total_; }
    /// Per-component weight accepted by the common grid but outside its native
    /// component axis. A correctly configured fit tally keeps every entry zero.
    const std::array<double, kFluxComponentCount>&
    component_dropped_weight() const { return component_dropped_weight_; }

private:
    void end_history_impl(int logical_replica);
    void add_coarse(int energy_bin, double mu_in, double weight);

    FluenceAxes axes_;
    std::array<FluenceAxes, static_cast<int>(FluxOrigin::Count)>
        component_axes_;
    std::array<std::size_t, kFluxComponentCount + 1> component_offsets_{};
    uint64_t n_histories_ = 0;

    std::vector<LineDef> lines_;

    // Line channels, laid out exactly like the grid: [n_lines * n_mu] and
    // [n_lines * n_mu * ChainClass::Count], plus per-line history totals so
    // line_total_sigma() is exact rather than a quadrature of bin variances.
    std::vector<double> line_sum_w_;
    std::vector<double> line_sum_w2_;
    std::vector<double> line_class_sum_w_;
    std::vector<double> line_class_sum_w2_;
    std::vector<double> line_total_w_;
    std::vector<double> line_total_w2_;
    std::vector<double> hist_line_w_;
    std::vector<int> hist_line_touched_;
    std::vector<double> hist_line_class_w_;
    std::vector<int> hist_line_class_touched_;
    std::vector<double> hist_line_total_;

    // Running accumulators over histories.
    std::vector<double> sum_w_;    ///< per bin
    std::vector<double> sum_w2_;   ///< per bin, per-history-summed then squared
    // Per chain class, summed over bins (for the decomposition/bounce measure).
    std::vector<double> class_sum_w_;   ///< [n_bins * ChainClass::Count]
    std::vector<double> class_sum_w2_;  ///< per-history-summed then squared
    // Per creation origin: continuum grid sums (mean-only) plus per-history
    // totals so the per-origin sigma is exact under intra-history correlation.
    std::vector<double> origin_sum_w_;  ///< [n_bins * FluxOrigin::Count]
    double origin_total_w_[static_cast<int>(FluxOrigin::Count)] = {};
    double origin_total_w2_[static_cast<int>(FluxOrigin::Count)] = {};
    double hist_origin_total_[static_cast<int>(FluxOrigin::Count)] = {};

    // Fit components, component-major. Unlike the per-origin grid, these carry
    // exact per-cell and integrated S2 after per-history aggregation.
    std::vector<double> component_sum_w_;
    std::vector<double> component_sum_w2_;
    std::vector<double> hist_component_w_;
    std::vector<int> hist_component_touched_;
    std::vector<double> component_total_w_;
    std::vector<double> component_total_w2_;
    std::vector<double> hist_component_total_;

    // Exact-line counterparts, component-major [component][line][mu].
    std::vector<double> line_component_sum_w_;
    std::vector<double> line_component_sum_w2_;
    std::vector<double> hist_line_component_w_;
    std::vector<int> hist_line_component_touched_;

    // Logical-replica aggregate maps for downstream covariance estimates.
    std::array<uint64_t, kNumLogicalReplicas> replica_histories_{};
    std::array<std::vector<float>, kNumLogicalReplicas>
        replica_continuum_sum_w_;
    std::array<std::vector<float>, kNumLogicalReplicas>
        replica_component_sum_w_;
    std::array<std::vector<float>, kNumLogicalReplicas>
        replica_line_sum_w_;
    std::array<std::vector<float>, kNumLogicalReplicas>
        replica_line_component_sum_w_;

    // Per-history coarse observables and their exact raw first/outer moments.
    std::array<double, kCoarseDimension> hist_coarse_{};
    std::array<double, kCoarseDimension> coarse_sum_w_{};
    std::array<double, kCoarseDimension * kCoarseDimension>
        coarse_outer_sum_w_{};

    std::array<double, kFluxComponentCount> component_dropped_weight_{};

    // Per-history TOTAL accumulators. total_sigma() needs the variance of the
    // history total; summing per-bin variances would under-state it, because one
    // history scores into many bins at once with positively correlated
    // magnitudes.
    double sum_total_w_ = 0.0;
    double sum_total_w2_ = 0.0;

    // Buffer for the history in progress.
    double hist_total_ = 0.0;
    std::vector<double> hist_w_;
    std::vector<int> hist_touched_;     ///< indices touched this history
    std::vector<double> hist_class_w_;  ///< [n_bins * Count] contributions
    std::vector<int> hist_class_touched_;

    double dropped_weight_ = 0.0;
    double dropped_weight_sum_w2_ = 0.0;
    double hist_dropped_weight_ = 0.0;
    double max_history_total_ = 0.0;
};

} // namespace ceelo
