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

/// @file NextEventEstimator.h
/// @brief Next-event (point-detector) estimator for air and environment scatter.
///
/// WHY. Brute-force MC cannot converge the scattered flux arriving at a small
/// detector far from a source: scatter-in is wide-angle, so cone biasing (which
/// only emits toward the detector) structurally misses it, and isotropic
/// emission almost never hits: a cone-scale series is still climbing at 5 m and
/// nowhere near converged at 50 m (measured on the air benchmark). A next-event
/// estimator removes the "must randomly hit the detector" variance: at EVERY
/// scatter vertex it deterministically scores the attenuated probability of
/// reaching the detector.
///
/// WHAT IS SCORED. A binned incident-fluence field phi(E', mu_in) at the
/// detector location, per source photon, in 1/cm^2. NOT an efficiency: no
/// detector response is applied here (see FluenceTally.h for why folding one in
/// would be a correctness error, and because the detector-response fold is
/// deliberately the NEXT phase of this work). The uncollided/direct term is
/// reported separately as a scalar transmission — the direct and scattered terms
/// are mutually exclusive outcomes of one primary, not co-emitted photons.
///
/// THE BOUNDED ESTIMATOR. A pure point-detector estimator carries a 1/r^2 that
/// diverges for vertices near the detector, giving unbounded variance — and
/// near-detector vertices dominate the contribution, so an ad-hoc cap would bias
/// the dominant term. Instead the detector is scored as a SPHERE of radius
/// `scoring_radius_cm`, and what is scored is the expected TRACK LENGTH inside
/// that sphere:
///
///     phi = w * P_branch * integral_cone p(Omega) * l(Omega) dOmega
///             * exp(-tau) / V
///
/// where the cone is the sphere's exact subtension from the vertex, l(Omega) is
/// the chord that direction cuts through the sphere, and V is its volume. This
/// is the standard track-length flux estimator, and it is bounded (the cone
/// integral of a normalized density cannot exceed 1, and l <= 2R) with no
/// singularity at any r.
///
/// Scoring track length rather than CROSSINGS matters. Converting crossings to
/// fluence needs the mean-chord identity (mean chord = 4V/S = 4R/3), which holds
/// for a parallel or isotropic incident field but not for a nearby point vertex,
/// whose rays clip the sphere at grazing incidence and cut systematically
/// shorter chords. Since near vertices dominate, that approximation biased the
/// result high by ~4% in a ground-plane geometry (measured against the analog
/// truth in tests/test_env_scatter.cpp) before the chord form replaced it. The
/// chord form is exact at any distance, including vertices INSIDE the sphere,
/// where the chord is simply the distance to the exit.
///
/// The cone integral is evaluated by sampling one direction uniformly in the
/// cone — unbiased at any r, essentially exact in the far field where the cone
/// is tiny, and it assigns the incidence angle mu_in correctly by construction.
///
/// MULTIPLE SCATTER IS MANDATORY, not optional: in air at 50 m / 662 keV the
/// optical depth is ~0.46, so double/single scatter is ~23%, and in a dense slab
/// the albedo is intrinsically a multiple-collision quantity. The walk therefore
/// continues analog after each score, and scores again at every subsequent
/// vertex.

#include "geometry/EnvironmentGeometry.h"
#include "transport/FluenceTally.h"

#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace ceelo {

/// Differential scattering coefficient of a material, per steradian, evaluated
/// at one scattering cosine. This is the physics primitive the next-event score
/// is built from:
///
///     dMu_s/dOmega(mu_s) = sum_e mu_{s,e} * p_e(mu_s)
///
/// where p_e is the normalized per-steradian angular density of element e.
/// Summing over the material's elements (rather than sampling one) is
/// deterministic and lower variance, and it matches how Material::macroscopic_xs
/// builds mu. Integrating either component over the full sphere must return the
/// material's mu_cs / mu_rs exactly — the closure identity that guarantees the
/// estimator is normalized consistently with the analog walk it accompanies.
struct DiffScatter {
    double compton = 0.0;   ///< 1/cm/sr
    double rayleigh = 0.0;  ///< 1/cm/sr
};

DiffScatter differential_scatter(const Material& mat, double energy_keV,
                                 double mu_s, bool include_rayleigh = true,
                                 bool free_electron = false);

/// How the primary's first interaction site is produced.
enum class PrimaryMode : uint8_t {
    /// Isotropic 4pi emission followed by an analog flight to the first
    /// collision. Simple and assumption-free; currently the ONLY implemented
    /// mode, and the one every validation result was produced with.
    Analog = 0,
    /// NOT YET IMPLEMENTED -- `run()` throws if this is selected, rather than
    /// silently running Analog.
    ///
    /// Intended design: sample the first collision vertex directly from an
    /// importance density concentrated on the region that actually contributes
    /// (a prolate halo around the source->detector line, weighted toward the
    /// detector end, since both 1/r^2 and exp(-tau) favor near-detector
    /// vertices), with the source->vertex leg becoming a deterministic
    /// attenuated weight w1 = mu_tot*exp(-tau)/(4 pi r^2) / q(x). This is a bake
    /// COST optimization, not a correctness issue: Analog is unbiased but leaves
    /// large variance in open air at long range, where most primaries fly off
    /// and contribute nothing while rare near-detector vertices contribute a lot.
    VertexSampled = 1
};

/// Experimental, unbiased variance-reduction policies.  Both are deliberately
/// OFF by default and remain benchmark-only until the campaign acceptance gates
/// are demonstrated.  They are separate switches because their cost/variance
/// tradeoffs are different and must be measured independently.
enum class ExperimentalPrimarySampling : uint8_t {
    Analog = 0,
    /// One logical primary history transports a direction and its antipode,
    /// each at weight 1/2.  Both walks are closed by ONE tally.end_history(),
    /// so S2 contains the square of their combined contribution.
    AntipodalPair = 1
};

enum class ExperimentalConeSampling : uint8_t {
    Single = 0,
    /// Each detector-cone integral uses two measure-preserving antithetic
    /// points, each at integration weight 1/2.  The policy is applied to
    /// Compton/Rayleigh, annihilation, and brems next-event scores together.
    AntitheticPair = 1
};

struct EnvScatterConfig {
    /// Radius of the spherical scoring body at the detector location. Should be
    /// comparable to the detector's own size (a 3"x3" NaI has R = 3.81 cm).
    double scoring_radius_cm = 4.0;

    /// Walk termination.
    double e_floor_keV = 30.0;   ///< Library validity floor.
    /// Optional half-size of a cubic world centred on the detector, in cm.
    /// Zero (the default) means the air and any slabs are UNBOUNDED, which is
    /// the physical idealization and what production bakes should use.
    ///
    /// It exists because GEANT4 cannot represent an unbounded world: the
    /// exported GDML puts the ground and the air in a finite box (a 3000 cm
    /// cube), so G4 kills any chain that leaves it while CeeLo keeps
    /// integrating. Scattered flux falls off roughly as r^-3 per unit radius,
    /// so the truncated tail is small but not zero, and it biases the
    /// comparison in ONE direction: CeeLo high. Set this to the GDML world's
    /// half-size to make a CeeLo-vs-G4 comparison geometry-matched; leave it at
    /// zero for physics.
    double world_half_size_cm = 0.0;
    int max_vertices = 50;       ///< Matches SourceGeometry's kMaxSourceScatters.
    double weight_floor = 1e-5;  ///< Below this, Russian-roulette the history.
    /// Roulette a history whose remaining optical depth to the nearest air-facing
    /// slab surface exceeds this: it is buried too deep to re-emerge. Unbiased
    /// here because the tally is linear (the FEP non-linearity hazard that makes
    /// branch roulette unsafe in the main engine does not exist in a fluence
    /// tally).
    double roulette_tau = 8.0;

    bool include_rayleigh = true;
    /// Model the 511 keV pair from pair production above threshold: score its
    /// uncollided arrival into the Annihilation LINE channel AND walk both
    /// photons so their own scatter chains reach the continuum.
    ///
    /// Walking them is not optional at high energy. With only the uncollided
    /// term (the original behavior) the kernel measured 8.3% BELOW GEANT4 at
    /// 2614 keV, where pair production is 3.4% of interactions in soil.
    ///
    /// Remaining approximation: annihilation IN FLIGHT is not modeled -- the
    /// positron is always assumed to stop first. That is ~2-3% of positrons at
    /// these energies, reshaping <~0.3% of the total tally, below the GEANT4
    /// campaign's resolution. It is a one-sided approximation (G4 option4 does
    /// model it) and matters only if the 2614 keV target tightens below ~1%.
    bool include_annihilation = true;

    /// Walk the annihilation photons rather than only scoring their uncollided
    /// arrival. A knob purely for A/B quantification; production wants it on.
    bool walk_annihilation = true;

    /// Displace the 511 birth point to where the positron actually stops.
    /// The positron CSDA range at the mean pair energy from a 2614 keV photon
    /// is 0.26 cm in soil but 310 cm in AIR -- and 40% of pair-production
    /// vertices are in air, because primaries traverse a lot of it before
    /// reaching the ground. Annihilating at the vertex therefore misplaces a
    /// substantial part of the channel by a large fraction of the
    /// source-detector distance. `false` restores annihilate-at-vertex for A/B.
    bool positron_range = true;

    /// Model bremsstrahlung from the electrons the photon walk creates
    /// (Compton recoils, photoelectrons, both pair electrons). Measured
    /// against the G4 harness's creator-process origin decomposition, brems
    /// is 4.26% of the scattered flux at 2614 keV -- the entire remaining
    /// deficit there -- and although it is a small fraction of the TOTAL, it
    /// is 100% of the 30-33 keV bin and ~half of everything below 60 keV, so
    /// a spectral fit trained without it learns "no flux" where reality is
    /// twice the model. Each emission is scored next-event through the 2BS
    /// angular density about the electron's DIFFUSED direction at the true
    /// along-walk emission point (metres from the vertex in air), and the
    /// analog photon is pushed so its own scatter chain scores.
    bool include_brems = true;

    /// Electrons below this kinetic energy do not radiate (cost gate). Their
    /// photons are almost all below the 30 keV tally floor anyway: the 2BS/SB
    /// spectrum is ~1/k and an electron's radiative yield scales ~KE^2, so
    /// sub-150 keV electrons contribute negligibly to the >= 30 keV band.
    /// Final value set from the origin-tagged 662/122 keV GEANT4 runs.
    double brems_electron_floor_keV = 150.0;

    /// Walk the brems photons rather than only scoring their uncollided
    /// arrival. A knob purely for A/B quantification; production wants it on.
    bool walk_brems = true;
    /// DIAGNOSTIC ONLY: use the free-electron Klein-Nishina angular
    /// distribution instead of KN x S(x,Z), in BOTH the scoring and the walk.
    /// S(x,Z) suppresses small-angle (low momentum transfer) scatter, so an
    /// error in it redistributes weight between forward and wide angles while
    /// leaving the total mu_cs fixed. Toggling it bounds how much of a
    /// forward/wide-angle discrepancy S can possibly account for. It is NOT a
    /// physical option -- binding is real -- and must stay false in production.
    bool free_electron_compton = false;

    /// Doppler-broaden the scored Compton energy. Off by default: it costs time
    /// and only reshapes bins far narrower than this kernel's energy grid.
    bool doppler = false;

    PrimaryMode primary_mode = PrimaryMode::Analog;

    /// BENCHMARK ONLY.  Production manifests and payload bakes must leave both
    /// at their default values unless a paired common-random-number campaign
    /// has passed the bias, closure, and variance x wall-time gates.
    ExperimentalPrimarySampling experimental_primary_sampling =
        ExperimentalPrimarySampling::Analog;
    ExperimentalConeSampling experimental_cone_sampling =
        ExperimentalConeSampling::Single;

    FluenceAxes axes;
};

/// The discrete line channels a node carries: the elastic-E0 line always, plus
/// the 511 annihilation line when the primary is above the pair threshold.
std::vector<LineDef> env_line_channels(double e0_keV);

/// Exact transport-event counters carried into fit-block diagnostics.  These
/// count events in the analog walk (not next-event scores), so they are useful
/// for closure/performance audits without changing the estimator observable.
enum class EnvTransportCounter : std::uint8_t {
    Photoelectric = 0,
    Compton,
    Rayleigh,
    PairProduction,
    BremsEmission,
    Count
};

struct EnvTransportDiagnostics {
    std::array<std::uint64_t,
               static_cast<std::size_t>(EnvTransportCounter::Count)>
        interaction_counts{};
    /// Photon chains that remained alive after exactly max_vertices steps.
    /// This is an explicit transport truncation, distinct from a campaign's
    /// history cap and from roulette termination.
    std::uint64_t max_vertex_terminations = 0;
    /// Reserved for recoverable transport failures. Fatal exceptions still
    /// fail the worker block rather than being hidden in this counter.
    std::uint64_t failures = 0;

    void add_interaction(EnvTransportCounter counter,
                         std::uint64_t count = 1)
    {
        interaction_counts[static_cast<std::size_t>(counter)] += count;
    }

    void merge(const EnvTransportDiagnostics& other)
    {
        for (std::size_t i = 0; i < interaction_counts.size(); ++i)
            interaction_counts[i] += other.interaction_counts[i];
        max_vertex_terminations += other.max_vertex_terminations;
        failures += other.failures;
    }
};

/// Everything the bake needs from one (geometry, energy) node.
struct EnvScatterResult {
    FluenceTally tally;
    EnvTransportDiagnostics diagnostics;
    /// Uncollided transmission source -> detector at E0 (the peak factor).
    double direct_transmission = 0.0;
    /// Uncollided fluence per source photon, T/(4 pi d^2), in 1/cm^2.
    double direct_fluence = 0.0;
    double source_distance_cm = 0.0;
    uint64_t histories = 0;

    /// `e0_keV` fixes the line channels; run() validates the tally against the
    /// energy it is asked to run, so a result cannot be silently reused across
    /// energies with the wrong elastic line.
    EnvScatterResult(const FluenceAxes& axes, double e0_keV,
                     std::vector<FluenceAxes> component_axes = {})
        : tally(axes, env_line_channels(e0_keV),
                std::move(component_axes)), e0_keV(e0_keV) {}

    double e0_keV = 0.0;
};

class NextEventEstimator {
public:
    /// @param env       The environment (air + slabs).
    /// @param cfg       Estimator configuration.
    /// @param source    Source position, detector frame.
    /// @param detector  Detector reference point (normally the origin).
    /// @param det_axis  Detector symmetry axis; mu_in is measured against it.
    ///                  Default +z, matching CeeLo's detector frame (crystal
    ///                  face at z = 0 extending toward +z), so mu_in = +1 is a
    ///                  photon arriving through the front face and mu_in = -1
    ///                  one arriving through the back.
    NextEventEstimator(const EnvironmentGeometry& env,
                       const EnvScatterConfig& cfg,
                       const Eigen::Vector3d& source,
                       const Eigen::Vector3d& detector = Eigen::Vector3d::Zero(),
                       const Eigen::Vector3d& det_axis = Eigen::Vector3d::UnitZ());

    /// Run `n_histories` primary photons at `e0_keV`, accumulating into `out`.
    /// Thread-safe in the sense that separate instances/tallies may run in
    /// parallel and be merged; a single call is single-threaded.
    void run(double e0_keV, uint64_t n_histories, uint64_t seed,
             EnvScatterResult& out, int logical_replica = -1) const;

    /// INDEPENDENT ANALOG TRUTH — no next-event scoring at all.
    ///
    /// Photons are walked fully analog and the track length inside the scoring
    /// sphere is tallied directly (phi = L / V, the standard track-length flux
    /// estimator). In expectation this equals the crossing form the next-event
    /// estimator produces, since for a convex body crossings * mean_chord =
    /// track length and (4R/3)/(4 pi R^3/3) = 1/(pi R^2).
    ///
    /// This converges far too slowly to be useful at real distances — which is
    /// the entire reason the next-event estimator exists — but at short range it
    /// converges fine, and comparing the two is the unbiasedness gate for the
    /// estimator. Only POST-SCATTER track length is tallied, so it is directly
    /// comparable to the scattered fluence tally.
    void run_brute(double e0_keV, uint64_t n_histories, uint64_t seed,
                   EnvScatterResult& out, int logical_replica = -1) const;

    /// Deterministic direct terms (no MC): uncollided transmission and fluence.
    double direct_transmission(double e0_keV) const;

    const EnvScatterConfig& config() const { return cfg_; }

private:
    /// One photon awaiting its walk: the primary, the 511 pair from a pair-
    /// production vertex, and any brems photons the walk's electrons emit.
    /// Depth stays small in practice -- every secondary is strictly below its
    /// parent's energy, so the generation chain terminates fast.
    ///
    /// NOTE: initialize by explicit member assignment, not positional
    /// aggregate init. A field inserted mid-struct silently shifts positional
    /// initializers (line_e0 landing on slab_hits, elastic_line falling back
    /// to -1), which empties the ElasticE0 line channel while leaving every
    /// total unchanged.
    struct PendingPhoton {
        Eigen::Vector3d pos{0, 0, 0};
        Eigen::Vector3d dir{0, 0, 1};
        double energy = 0.0;
        double weight = 1.0;
        bool touched_slab = false;
        int first_slab = -1;
        bool multi_slab = false;
        /// Number of SLAB interactions so far. Distinct from `multi_slab`:
        /// superposing per-surface kernels can only reproduce chains with
        /// exactly one bounce, so the bounce COUNT is what decides which term
        /// of the room decomposition a chain belongs to.
        int slab_hits = 0;
        /// Line channel this chain still belongs to while it has undergone ONLY
        /// elastic (Rayleigh) scatters; -1 once any Compton has occurred.
        int elastic_line = -1;
        /// Creation origin, inherited by the whole scatter chain and
        /// overwritten at each secondary creation (G4 creator semantics).
        FluxOrigin origin = FluxOrigin::Primary;
        /// Interactions undergone by THIS photon since its creator. Secondary
        /// creation resets it to zero; it saturates at two for fit tagging.
        int post_creation_interactions = 0;
    };

    struct Walker;  // per-history scratch (segment buffers, photon stack)

    void run_history(double e0_keV, std::mt19937_64& rng, Walker& w,
                     FluenceTally& tally,
                     EnvTransportDiagnostics& diagnostics) const;

    /// Score the next-event contribution of one vertex.
    void score_vertex(const Eigen::Vector3d& vertex,
                      const Eigen::Vector3d& dir_in,
                      double energy_keV,
                      double weight,
                      const Material& mat,
                      ChainClass cls,
                      int elastic_line,
                      PostCreationOrder order,
                      std::mt19937_64& rng,
                      FluenceTally& tally,
                      EnvTraceScratch& attenuation_scratch,
                      FluxOrigin origin = FluxOrigin::Primary) const;

    /// Score an isotropic emission of `n_photons` at `energy_keV` from a vertex
    /// (used for the 511 keV annihilation pair).
    void score_isotropic(const Eigen::Vector3d& vertex,
                         double energy_keV,
                         double weight,
                         ChainClass cls,
                         int line,
                         std::mt19937_64& rng,
                         FluenceTally& tally,
                         EnvTraceScratch& attenuation_scratch,
                         FluxOrigin origin = FluxOrigin::Primary) const;

    /// Form a brems-radiating electron at a vertex: sample its emission
    /// events along the diffused geometry-free walk (true emission
    /// positions), next-event score each photon's uncollided arrival through
    /// the 2BS density about the electron's direction AT EMISSION, and push
    /// the analog-sampled photon so its own scatter chain scores. The pushed
    /// photon's uncollided first flight is never tallied by the walk (the
    /// walk scores only at collisions), so it is covered exactly once, by
    /// the emission score here -- the same partition the annihilation
    /// channel uses.
    void score_brems_electron(const Eigen::Vector3d& vertex,
                              const Eigen::Vector3d& electron_dir,
                              double electron_KE_keV,
                              double weight,
                              const Material& mat,
                              bool is_positron,
                              const PendingPhoton& parent,
                              ChainClass cls,
                              std::mt19937_64& rng,
                              Walker& walk,
                              FluenceTally& tally,
                              EnvTransportDiagnostics& diagnostics) const;

    /// Expected UNCOLLIDED track length for a geometric chord `chord_cm` of air
    /// at `energy_keV`: (1 - exp(-mu*c))/mu. Scoring the raw chord would
    /// over-count by ~mu*c/2, which is negligible for a small sphere but reaches
    /// ~0.5% at the 50 cm radius wanted for far-field GEANT4 gating.
    double attenuated_chord(double chord_cm, double energy_keV) const;

    /// Solid angle of the scoring sphere from `vertex`, and its axis.
    /// Returns cos(half-angle); dOmega = 2*pi*(1 - cos_half).
    double scoring_cone(const Eigen::Vector3d& vertex,
                        Eigen::Vector3d& axis_out) const;

    const EnvironmentGeometry& env_;
    EnvScatterConfig cfg_;
    Eigen::Vector3d source_;
    Eigen::Vector3d detector_;
    Eigen::Vector3d det_axis_;
    double scoring_volume_cm3_;  ///< (4/3) pi R^3
};

} // namespace ceelo
