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

/// @file EnvScatterKernel.h
/// @brief The product layer: incident flux onto a detector's surfaces, and
///        storage/interpolation of baked environment-scatter kernels.
///
/// SCOPE. This phase delivers the photon FLUX ARRIVING AT the detector, resolved
/// by energy and incidence angle and split by which surface it enters through.
/// It deliberately does NOT apply any detector efficiency — folding in the
/// detector response is the next phase. Everything here is therefore a property
/// of the source, the environment, and the detector's OUTER SHAPE only, and is
/// reusable across crystals.
///
/// RAYLEIGH PAIRING RULE — read before combining a kernel with a peak factor.
/// The peak term and this scattered kernel must agree about who owns coherent
/// scatter, or forward Rayleigh is counted twice:
///
///   * `AirAttenuation::AnalyticNoScatter` already treats Rayleigh out-scatter as
///     compensated by in-scatter (the open-geometry argument), so the peak it
///     produces ALREADY INCLUDES forward coherent in-scatter. Pair it with a
///     kernel baked using `--no-rayleigh`.
///   * `AirAttenuation::AnalyticCollimated` removes Rayleigh from the peak, so
///     coherent scatter is unaccounted for and the kernel must supply it. Pair it
///     with a normal (Rayleigh-on) kernel.
///
/// Mixing `AnalyticNoScatter` with a Rayleigh-on kernel double-counts the
/// near-peak coherent flux. Minor at 662 keV (only ~0.3% of scattered flux lands
/// within ~10 keV of the line) but growing to several percent near 60 keV, where
/// Rayleigh dominates the air cross-section.

#include "transport/FluenceTally.h"

#include <string>
#include <vector>

namespace ceelo {

/// A right-circular-cylinder detector envelope, in the CeeLo detector frame:
/// symmetry axis along +z, front face at z = 0, body extending toward +z.
struct DetectorEnvelope {
    double radius_cm = 3.81;   ///< 3" diameter default (a 3"x3" NaI).
    double length_cm = 7.62;
};

/// Photon flux entering each surface of the envelope, per source photon.
///
/// Convention: mu_in is the cosine between the photon's DIRECTION OF TRAVEL and
/// the detector axis (+z). A photon arriving head-on through the front face has
/// mu_in = +1; one arriving from behind, travelling in -z, enters the back with
/// mu_in = -1; mu_in = 0 is pure side-on.
///
/// For a parallel beam at incidence cosine mu the projected areas are
///     face/back: pi * R^2 * |mu|      (face if mu > 0, back if mu < 0)
///     side:      2 * R * L * sqrt(1 - mu^2)
/// so each surface's incident photon count is the fluence times its projected
/// area, integrated over the angular distribution.
struct SurfaceFlux {
    double face = 0.0;   ///< Photons entering the front face, per source photon.
    double side = 0.0;   ///< Photons entering the curved side.
    double back = 0.0;   ///< Photons entering the rear face.

    double total() const { return face + side + back; }
};

/// Per-surface flux integrated over all energies, CONTINUUM PLUS ALL LINES.
SurfaceFlux surface_flux(const FluenceTally& tally, const DetectorEnvelope& det);

/// Per-surface flux of ONE discrete line channel. Kept separate because a line
/// and a continuum of equal magnitude produce very different detector
/// responses -- the line makes a photopeak. This is the whole reason the line
/// channels exist.
SurfaceFlux line_surface_flux(const FluenceTally& tally,
                              const DetectorEnvelope& det, int line);

/// Per-surface flux restricted to one energy bin (index into tally.axes()).
SurfaceFlux surface_flux_in_bin(const FluenceTally& tally,
                                const DetectorEnvelope& det, int energy_bin);

/// 1-sigma statistical uncertainty on `surface_flux`, propagated from the
/// per-bin tally variances (added in quadrature — slightly conservative, the
/// same convention FluenceTally::total_sigma uses).
SurfaceFlux surface_flux_sigma(const FluenceTally& tally,
                               const DetectorEnvelope& det);

/// Per-surface flux of the UNCOLLIDED beam. The direct term is a single
/// direction, not a distribution, so it is projected analytically rather than
/// binned — and it is kept out of the scattered tally because the two are
/// mutually exclusive outcomes of one primary, never co-tallied.
SurfaceFlux direct_surface_flux(double direct_fluence, double mu_in,
                                const DetectorEnvelope& det);

// ---------------------------------------------------------------------------
// Baked kernel storage
// ---------------------------------------------------------------------------

/// One baked node: the scattered-fluence field at one (geometry, energy) point.
///
/// Values are stored NORMALIZED to the uncollided direct fluence,
/// s_hat = phi_scat / (T(E0) / (4 pi d^2)), which makes them dimensionless,
/// O(1e-3 .. 1), and smooth in every axis — far better behaved for
/// interpolation than the raw fluence, which spans decades with 1/d^2.
/// The peak transmission T(E0, rho*d) is NOT tabulated: it is analytic.
struct EnvKernelNode {
    double e0_keV = 0.0;
    double distance_cm = 0.0;
    double h_src_cm = 0.0;      ///< 0 for the air-only kernel.
    double h_det_cm = 0.0;
    double transverse_extent_cm = 0.0;  ///< Air kernel only; 0 = unbounded.

    double direct_transmission = 0.0;
    double direct_fluence = 0.0;
    uint64_t histories = 0;
    double max_rel_sigma = 0.0;         ///< Worst per-bin relative sigma.

    /// s_hat[energy_bin * n_mu + mu_bin].
    std::vector<double> s_hat;
    /// Same, but only chains that touched a slab (the albedo part). Empty for
    /// the air kernel. Kept separate so a room can be assembled by superposing
    /// per-surface albedo kernels on top of one air kernel without double
    /// counting the pure-air chains.
    std::vector<double> s_hat_slab;
    /// Discrete line channels: flux arriving at an EXACT energy (all-elastic
    /// chains at E0, and the 511 keV annihilation pair). Stored apart from the
    /// continuum grid so a consumer can make a photopeak from them instead of
    /// smearing them across a ~10%-wide log bin.
    struct Line {
        double energy_keV = 0.0;
        uint8_t origin = 0;              ///< LineOrigin
        std::vector<double> s_hat;       ///< [n_mu], normalized to direct fluence
        std::vector<double> s_hat_slab;  ///< [n_mu], slab-touching chains only
    };
    std::vector<Line> lines;

    /// Fraction of the scattered fluence from chains touching TWO OR MORE slabs.
    /// The superposition model neglects these, so this MEASURES the size of what
    /// is being neglected rather than assuming it.
    double multi_slab_fraction = 0.0;

    /// PER-NODE energy axis. A landmark-blocked grid pins its bin edges at
    /// this node's OWN backscatter edge and E0, so the axis is a property of
    /// the node, not of the table. Empty means "use the table-level axes"
    /// (a log-uniform axis).
    std::vector<double> energy_edges;

    /// Number of energy bins this node actually carries (edges.size()-1 when
    /// `energy_edges` is populated, else the table axes' n_energy).
    int n_energy(const FluenceAxes& table_axes) const
    {
        return energy_edges.empty() ? table_axes.n_energy
                                    : static_cast<int>(energy_edges.size()) - 1;
    }
};

/// A baked table: a set of nodes on a structured grid, plus provenance.
///
/// Storage follows the VirtualDepthFit precedent: self-describing plain text,
/// with the axes, the validity ranges, the material identities and the per-node
/// statistics all stamped into the file, so a consumer can tell exactly what it
/// is holding and refuse to extrapolate outside it.
struct EnvKernelTable {
    std::string kind;            ///< "air", "soil", "concrete".
    std::string air_name;
    double air_density = 0.0;
    std::string slab_name;       ///< Empty for the air kernel.
    double slab_density = 0.0;
    double scoring_radius_cm = 0.0;
    FluenceAxes axes;

    std::vector<double> e0_grid_keV;
    std::vector<double> distance_grid_cm;
    std::vector<double> h_src_grid_cm;
    std::vector<double> h_det_grid_cm;
    std::vector<double> extent_grid_cm;

    std::vector<EnvKernelNode> nodes;

    /// Human-readable note recorded in the header (generator, date, caveats).
    std::string provenance;
};

/// Write/read the self-describing text form. Returns false on I/O failure.
/// `load_env_kernel` auto-detects the binary form below by its magic, so
/// consumers have a single load path for either representation.
bool save_env_kernel(const EnvKernelTable& table, const std::string& path);
bool load_env_kernel(EnvKernelTable& table, const std::string& path);

/// Binary float32 form -- the production container (the text form is ~6x the
/// 20-30 MB vendorable budget at the 48x32 grid).
///
///   * Payload arrays (s_hat, s_hat_slab, line channels) are float32: ~1e-7
///     relative, three orders below the ~1% statistical floor of any bake.
///     Scalars/metadata stay full precision.
///   * PER-NODE E0 TRUNCATION: energy bins strictly above the node's E0 are
///     structurally zero (Compton loses energy, Rayleigh preserves it, brems
///     photons are below their electron's KE < E0, and the 511 line is a
///     separate channel), so only bins up to the one containing E0 are
///     stored. The writer VERIFIES the dropped tail is zero and falls back
///     to the full array if not (so a future physics channel that breaks the
///     invariant degrades size, never correctness).
///   * Little-endian, as written by the only supported hosts. The magic
///     doubles as the detector: 'CEEK' + format version.
bool save_env_kernel_binary(const EnvKernelTable& table,
                            const std::string& path);
bool load_env_kernel_binary(EnvKernelTable& table, const std::string& path);

/// Merge a SHARD of a distributed bake into `table`.
///
/// A production bake is split across machines by node index (see the bake
/// driver's --shard/--nshards), so the results come back as several files that
/// must reassemble into one table. Requirements enforced here, because a
/// silently mismatched merge would poison a write-once product:
///   * identical kind/materials/scoring radius and identical axis GRIDS
///     (the geometry ladders and, for log-uniform tables, the energy
///     axis) -- a shard baked from a different grid is refused;
///   * no DUPLICATE nodes (same geometry+E0 coordinates): a duplicate means
///     overlapping shard assignments, i.e. some other node is missing, and
///     averaging them would hide that.
/// Returns false and leaves `table` untouched on any mismatch.
bool merge_env_kernel_shard(EnvKernelTable& table, const EnvKernelTable& shard);

/// True when `table` contains exactly one node per point of its declared
/// grids -- the completeness check to run before shipping a merged bake.
/// `missing` (optional) receives a human-readable list of absent coordinates.
bool env_kernel_is_complete(const EnvKernelTable& table,
                            std::string* missing = nullptr);

} // namespace ceelo
