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

#include "transport/NextEventEstimator.h"

#include "transport/AngularNorms.h"
#include "transport/ComptonScatter.h"
#include "transport/RayleighAngular.h"
#include "transport/TransportUtils.h"
#include "transport/PhotonTransport.h"
#include "cross_sections/CrossSectionData.h"
#include "physics/BremsAngular.h"
#include "physics/ElectronCsda.h"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ceelo {

namespace {

/// Where the ray from `origin` along `dir` (within [0, s_max]) meets the sphere.
struct SphereHit {
    double entry = 0.0;   ///< Distance to the sphere surface (0 if inside).
    double chord = 0.0;   ///< Length of the portion inside the sphere.
};

SphereHit sphere_hit(const Eigen::Vector3d& origin, const Eigen::Vector3d& dir,
                     double s_max, const Eigen::Vector3d& center, double radius)
{
    const Eigen::Vector3d oc = origin - center;
    const double b = 2.0 * oc.dot(dir);
    const double c = oc.squaredNorm() - radius * radius;
    const double disc = b * b - 4.0 * c;   // a == 1 for a unit direction
    if (disc <= 0.0) return {};
    const double sq = std::sqrt(disc);
    double t0 = 0.5 * (-b - sq);
    double t1 = 0.5 * (-b + sq);
    t0 = std::max(t0, 0.0);
    t1 = std::min(t1, s_max);
    SphereHit h;
    h.entry = t0;
    h.chord = std::max(0.0, t1 - t0);
    return h;
}

/// Length of the portion of the ray inside the sphere (convenience).
double chord_in_sphere(const Eigen::Vector3d& origin, const Eigen::Vector3d& dir,
                       double s_max, const Eigen::Vector3d& center, double radius)
{
    return sphere_hit(origin, dir, s_max, center, radius).chord;
}

/// Map (bounce count, spanned-distinct-slabs) onto a chain class.
ChainClass chain_class(int slab_hits, bool multi_distinct)
{
    if (slab_hits == 0) return ChainClass::AirOnly;
    if (slab_hits == 1) return ChainClass::OneBounce;
    return multi_distinct ? ChainClass::MultiDistinct : ChainClass::MultiSameSlab;
}

PostCreationOrder order_bucket(int interactions)
{
    if (interactions <= 0) return PostCreationOrder::Zero;
    if (interactions == 1) return PostCreationOrder::One;
    return PostCreationOrder::TwoPlus;
}

EnvTransportCounter diagnostic_counter(InteractionType interaction)
{
    switch (interaction) {
    case InteractionType::Photoelectric:
        return EnvTransportCounter::Photoelectric;
    case InteractionType::Compton:
        return EnvTransportCounter::Compton;
    case InteractionType::Rayleigh:
        return EnvTransportCounter::Rayleigh;
    case InteractionType::PairProduction:
        return EnvTransportCounter::PairProduction;
    }
    return EnvTransportCounter::Photoelectric;
}

constexpr double kTwoPi = 6.283185307179586;
constexpr double kFourPi = 12.566370614359172;
constexpr double kAnnihilationKeV = 510.998950;
constexpr double kPairThresholdKeV = 2.0 * kAnnihilationKeV;

/// Sample a direction uniformly inside a cone of half-angle acos(cos_half)
/// about `axis`.
Eigen::Vector3d sample_in_cone(const Eigen::Vector3d& axis, double cos_half,
                               std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double mu = 1.0 - uni(rng) * (1.0 - cos_half);
    const double phi = kTwoPi * uni(rng);
    return rotate_direction(axis, mu, phi);
}

/// Draw either the one-point cone rule or a two-point antithetic
/// rule.  The second point is the measure-preserving transform
/// (u_mu,u_phi) -> (1-u_mu,(u_phi+1/2) mod 1), so each point is marginally
/// uniform in the cone.  Averaging them is therefore unbiased.  Keeping the
/// one-point path separate preserves the production RNG stream bit-for-bit.
struct ConeDirections {
    std::array<Eigen::Vector3d, 2> direction;
    int count = 1;
};

ConeDirections sample_cone_directions(
    const Eigen::Vector3d& axis, double cos_half, std::mt19937_64& rng,
    ExperimentalConeSampling policy)
{
    ConeDirections result;
    if (policy == ExperimentalConeSampling::Single) {
        result.direction[0] = sample_in_cone(axis, cos_half, rng);
        return result;
    }

    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double u_mu = uni(rng);
    const double u_phi = uni(rng);
    const double width = 1.0 - cos_half;
    const double phi0 = kTwoPi * u_phi;
    const double phi1 = kTwoPi *
        ((u_phi < 0.5) ? (u_phi + 0.5) : (u_phi - 0.5));
    result.direction[0] =
        rotate_direction(axis, 1.0 - u_mu * width, phi0);
    result.direction[1] =
        rotate_direction(axis, 1.0 - (1.0 - u_mu) * width, phi1);
    result.count = 2;
    return result;
}

} // namespace

std::vector<LineDef> env_line_channels(double e0_keV)
{
    std::vector<LineDef> lines;
    lines.push_back({e0_keV, LineOrigin::ElasticE0});
    if (e0_keV >= 2.0 * 510.998950)
        lines.push_back({510.998950, LineOrigin::Annihilation});
    return lines;
}

namespace {

/// Fraction of the CSDA path length that shows up as net displacement. A
/// positron does not travel in a straight line, so its vector range is shorter
/// than its CSDA path; 0.6 is the standard projected-range fraction for ~1 MeV
/// electrons in low-Z media. Named and separated because it is the obvious knob
/// if the high-energy comparison leaves a residual.
constexpr double kPositronDetour = 0.6;

/// Where the positron from a pair-production vertex actually annihilates.
///
/// Consumes range as MASS thickness along the traced environment segments, so a
/// positron born just under the soil surface correctly emerges into air (where
/// it then travels metres) while one born deep in soil stops within millimetres.
/// Direction is the incoming photon's: pair-production secondaries are strongly
/// forward-peaked at these energies.
///
/// The positron kinetic energy is a PARAMETER (the caller draws the flat pair
/// split once and shares it): the same positron that radiates brems along its
/// flight is the one displacing the 511 birth point, so one draw keeps the two
/// coherent within an event.
Eigen::Vector3d annihilation_point(const EnvironmentGeometry& env,
                                   const Eigen::Vector3d& vertex,
                                   const Eigen::Vector3d& photon_dir,
                                   double positron_KE_keV,
                                   std::vector<EnvSegment>& seg_buf,
                                   EnvTraceScratch& trace_scratch)
{
    const double ke = positron_KE_keV;
    if (ke <= 0.0) return vertex;

    env.trace(vertex, photon_dir, seg_buf, trace_scratch);
    if (seg_buf.empty()) return vertex;

    // Track the FRACTION of the range consumed, not an absolute mass budget:
    // the range is material-dependent, so crossing into a new material rescales
    // what is left. A segment of length L in a material whose range is R_cm
    // consumes L/R_cm of the remaining fraction.
    double frac_used = 0.0;
    double travelled = seg_buf.back().t_end_cm;
    for (const EnvSegment& sg : seg_buf) {
        const double rho = sg.material->density();
        if (rho <= 0.0) continue;
        const double range_cm =
            ElectronCsda::instance().range_gcm2_material(
                *sg.material, ke, /*is_positron=*/true) * kPositronDetour / rho;
        if (!(range_cm > 0.0)) continue;

        const double seg_len = sg.t_end_cm - sg.t_start_cm;
        const double remaining_cm = (1.0 - frac_used) * range_cm;
        if (seg_len >= remaining_cm) {
            travelled = sg.t_start_cm + remaining_cm;   // stops inside this segment
            break;
        }
        frac_used += seg_len / range_cm;
        travelled = sg.t_end_cm;
    }
    return vertex + photon_dir * travelled;
}

} // namespace

DiffScatter differential_scatter(const Material& mat, double energy_keV,
                                 double mu_s, bool include_rayleigh,
                                 bool free_electron)
{
    const double e_MeV = energy_keV * 1e-3;
    const CrossSectionData& xs = CrossSectionData::instance();
    DiffScatter d;
    const auto& comp = mat.composition();
    for (std::size_t i = 0; i < comp.size(); ++i) {
        const int Z = comp[i].Z;
        const double n = mat.number_density(i) * kBarnToCm2;  // per barn -> 1/cm
        const auto sig = xs.all_cross_sections(Z, e_MeV);

        const int Zc = free_electron ? 0 : Z;   // Z=0 => free-electron KN
        const double norm_c = angular_norm_interp(ScatterProcess::Compton,
                                                  energy_keV, Zc);
        if (norm_c > 0.0) {
            const double p = compton_angular_pdf_unnorm(energy_keV, Zc, mu_s) /
                             (kTwoPi * norm_c);
            d.compton += n * sig.sigma_cs * p;
        }
        if (include_rayleigh) {
            const double norm_r = angular_norm_interp(ScatterProcess::Rayleigh,
                                                      energy_keV, Z);
            if (norm_r > 0.0) {
                const double p = rayleigh_angular_pdf_unnorm(energy_keV, Z, mu_s) /
                                 (kTwoPi * norm_r);
                d.rayleigh += n * sig.sigma_rs * p;
            }
        }
    }
    return d;
}


/// Per-history scratch, hoisted so the segment vectors are not reallocated on
/// every history.
struct NextEventEstimator::Walker {
    std::vector<EnvSegment> segs;
    std::vector<PendingPhoton> stack;
    EnvTraceScratch transport_trace;
    EnvTraceScratch detector_trace;
};

NextEventEstimator::NextEventEstimator(const EnvironmentGeometry& env,
                                       const EnvScatterConfig& cfg,
                                       const Eigen::Vector3d& source,
                                       const Eigen::Vector3d& detector,
                                       const Eigen::Vector3d& det_axis)
    : env_(env), cfg_(cfg), source_(source), detector_(detector),
      det_axis_(det_axis.normalized()),
      scoring_volume_cm3_((4.0 / 3.0) * M_PI * cfg.scoring_radius_cm *
                          cfg.scoring_radius_cm * cfg.scoring_radius_cm)
{
    // The scoring sphere must lie entirely in the air. If a slab cut through it,
    // chords inside that slab would be scored with the air attenuation (see
    // attenuated_chord) and, worse, the whole track-length-to-fluence relation
    // would mix media -- errors of tens of percent, silently. Catch it at
    // construction rather than let a bake produce quietly wrong numbers.
    for (std::size_t k = 0; k < env_.num_slabs(); ++k) {
        const EnvSlab& s = env_.slabs()[k];
        const double d = (detector_ - s.point).dot(s.normal);
        if (d < cfg_.scoring_radius_cm) {
            throw std::invalid_argument(
                "NextEventEstimator: scoring sphere (radius " +
                std::to_string(cfg_.scoring_radius_cm) + " cm) reaches slab " +
                std::to_string(k) + " at " + std::to_string(d) +
                " cm; reduce scoring_radius_cm or move the detector");
        }
    }
}

double NextEventEstimator::direct_transmission(double e0_keV) const
{
    return env_.transmission(source_, detector_, e0_keV * 1e-3);
}

double NextEventEstimator::scoring_cone(const Eigen::Vector3d& vertex,
                                        Eigen::Vector3d& axis_out) const
{
    const Eigen::Vector3d rel = detector_ - vertex;
    const double r = rel.norm();
    const double R = cfg_.scoring_radius_cm;
    if (r <= R) {
        // The vertex is inside the scoring sphere: the photon crosses the
        // surface whichever way it goes, so the cone is the full sphere.
        axis_out = Eigen::Vector3d::UnitZ();
        return -1.0;
    }
    axis_out = rel / r;
    const double sin_half = R / r;
    return std::sqrt(std::max(0.0, 1.0 - sin_half * sin_half));
}

void NextEventEstimator::score_vertex(const Eigen::Vector3d& vertex,
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
                                      FluxOrigin origin) const
{
    Eigen::Vector3d axis;
    const double cos_half = scoring_cone(vertex, axis);
    const double d_omega = kTwoPi * (1.0 - cos_half);

    const double mu_tot = mat.mu_total(energy_keV * 1e-3);
    if (!(mu_tot > 0.0)) return;

    const ConeDirections samples = sample_cone_directions(
        axis, cos_half, rng, cfg_.experimental_cone_sampling);
    const double integration_weight = 1.0 / samples.count;
    for (int sample = 0; sample < samples.count; ++sample) {
        const Eigen::Vector3d& dir_out = samples.direction[sample];
        const double mu_s = std::max(
            -1.0, std::min(1.0, dir_in.dot(dir_out)));
        const double mu_in = dir_out.dot(det_axis_);

        // Score the expected TRACK LENGTH through the scoring sphere, not the
        // number of crossings.  The sampled chord is exact even in near field.
        const SphereHit hit = sphere_hit(
            vertex, dir_out, std::numeric_limits<double>::max(), detector_,
            cfg_.scoring_radius_cm);
        if (hit.chord <= 0.0) continue;

        // Attenuation follows each sampled direction to its own entry point.
        // This remains necessary for the antithetic rule: averaging two center-
        // ray attenuations would not estimate the original cone integral.
        const Eigen::Vector3d entry_pt = vertex + dir_out * hit.entry;
        const DiffScatter d = differential_scatter(
            mat, energy_keV, mu_s, cfg_.include_rayleigh,
            cfg_.free_electron_compton);

        // Each antithetic point carries half the original quadrature weight.
        // Both contributions remain in the same primary-history buffer.
        const double common = weight * integration_weight * d_omega /
                              (mu_tot * scoring_volume_cm3_);
        const double e_out = compton_scattered_energy(energy_keV, mu_s);
        const bool score_compton =
            d.compton > 0.0 && e_out >= cfg_.e_floor_keV;
        const bool score_rayleigh =
            d.rayleigh > 0.0 && energy_keV >= cfg_.e_floor_keV;

        // Compton and Rayleigh at this cone point share one geometry trace.
        double compton_transmission = 0.0;
        double rayleigh_transmission = 0.0;
        if (score_compton && score_rayleigh) {
            const auto pair = env_.transmission_pair(
                vertex, entry_pt, e_out * 1e-3, energy_keV * 1e-3,
                attenuation_scratch);
            compton_transmission = pair[0];
            rayleigh_transmission = pair[1];
        } else if (score_compton) {
            compton_transmission = env_.transmission(
                vertex, entry_pt, e_out * 1e-3, attenuation_scratch);
        } else if (score_rayleigh) {
            rayleigh_transmission = env_.transmission(
                vertex, entry_pt, energy_keV * 1e-3, attenuation_scratch);
        }

        if (score_compton) {
            tally.add(e_out, mu_in,
                      common * d.compton * compton_transmission *
                          attenuated_chord(hit.chord, e_out),
                      cls, origin, order);
        }

        // Pure-Rayleigh descendants remain exact atoms.
        if (score_rayleigh) {
            const double w = common * d.rayleigh * rayleigh_transmission *
                             attenuated_chord(hit.chord, energy_keV);
            if (elastic_line >= 0)
                tally.add_line(elastic_line, mu_in, w, cls, origin, order);
            else
                tally.add(energy_keV, mu_in, w, cls, origin, order);
        }
    }
}

double NextEventEstimator::attenuated_chord(double chord_cm,
                                            double energy_keV) const
{
    // The analog truth tallies REAL track length, which is cut short when the
    // photon interacts inside the scoring sphere. Scoring the full geometric
    // chord therefore over-counts by ~mu*<chord>/2. That is only ~0.05% for a
    // 5 cm sphere in air, but it scales linearly with R: at the 50 cm sphere
    // wanted for far-field GEANT4 gating it reaches ~+0.5%, which would read as
    // a phantom G4 discrepancy. The exact expected uncollided track length in a
    // homogeneous medium is the closed form
    //     integral_0^c exp(-mu s) ds = (1 - exp(-mu c)) / mu.
    // The scoring sphere is required to lie entirely in the air (enforced in the
    // constructor), so a single mu applies.
    const double mu = env_.air().mu_total(energy_keV * 1e-3);
    const double x = mu * chord_cm;
    if (x < 1e-8) return chord_cm;           // mu -> 0 limit, avoids 0/0
    return (1.0 - std::exp(-x)) / mu;
}

void NextEventEstimator::score_isotropic(const Eigen::Vector3d& vertex,
                                         double energy_keV,
                                         double weight,
                                         ChainClass cls,
                                         int line,
                                         std::mt19937_64& rng,
                                         FluenceTally& tally,
                                         EnvTraceScratch& attenuation_scratch,
                                         FluxOrigin origin) const
{
    if (energy_keV < cfg_.e_floor_keV) return;
    Eigen::Vector3d axis;
    const double cos_half = scoring_cone(vertex, axis);
    const double d_omega = kTwoPi * (1.0 - cos_half);
    const ConeDirections samples = sample_cone_directions(
        axis, cos_half, rng, cfg_.experimental_cone_sampling);
    const double integration_weight = 1.0 / samples.count;
    for (int sample = 0; sample < samples.count; ++sample) {
        const Eigen::Vector3d& dir_out = samples.direction[sample];
        const SphereHit hit = sphere_hit(
            vertex, dir_out, std::numeric_limits<double>::max(), detector_,
            cfg_.scoring_radius_cm);
        if (hit.chord <= 0.0) continue;
        const double t_leg = env_.transmission(
            vertex, vertex + dir_out * hit.entry, energy_keV * 1e-3,
            attenuation_scratch);
        const double w = weight * integration_weight * (d_omega / kFourPi) *
                         attenuated_chord(hit.chord, energy_keV) * t_leg /
                         scoring_volume_cm3_;
        // The annihilation pair arrives at exactly 511 keV: an atom, not a
        // continuum bin.  Both antithetic scores retain that taxonomy.
        if (line >= 0)
            tally.add_line(line, dir_out.dot(det_axis_), w, cls, origin,
                           PostCreationOrder::Zero);
        else
            tally.add(energy_keV, dir_out.dot(det_axis_), w, cls, origin,
                      PostCreationOrder::Zero);
    }
}

void NextEventEstimator::score_brems_electron(const Eigen::Vector3d& vertex,
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
                                              EnvTransportDiagnostics& diagnostics) const
{
    if (!cfg_.include_brems ||
        electron_KE_keV < cfg_.brems_electron_floor_keV)
        return;

    // Geometry-free diffused walk in the vertex's material. Accepted
    // approximation: a track that would cross a material boundary (an
    // air-born electron arcing into the ground) keeps its birth material's
    // stopping power; the emitted photon's attenuation leg below IS traced
    // through the real geometry, which is the part that decides whether a
    // soft photon reaches the detector. The G4 harness's brems-by-material
    // origin split measures how much of the channel rides on this.
    const std::vector<BremsEmission> emissions =
        ElectronCsda::instance().sample_brems_emissions_in_material(
            mat, vertex, electron_dir, electron_KE_keV, rng, is_positron);
    if (emissions.empty()) return;
    diagnostics.add_interaction(
        EnvTransportCounter::BremsEmission,
        static_cast<std::uint64_t>(emissions.size()));

    const double eff_z = ElectronCsda::effective_Z(mat);
    for (const BremsEmission& e : emissions) {
        // Sub-band photons can never re-enter the >= e_floor band (Compton
        // loses energy, Rayleigh preserves it): skip both the score and the
        // walk so dropped_weight() stays a clean diagnostic.
        if (e.photon_keV < cfg_.e_floor_keV) continue;

        // Next-event score of the uncollided arrival, exactly the
        // score_isotropic pattern but with the sampler-faithful 2BS density
        // about the electron's diffused direction AT EMISSION replacing
        // 1/(4pi).
        Eigen::Vector3d axis;
        const double cos_half = scoring_cone(e.position, axis);
        const double d_omega = kTwoPi * (1.0 - cos_half);
        const ConeDirections samples = sample_cone_directions(
            axis, cos_half, rng, cfg_.experimental_cone_sampling);
        const double integration_weight = 1.0 / samples.count;
        for (int sample = 0; sample < samples.count; ++sample) {
            const Eigen::Vector3d& dir_out = samples.direction[sample];
            const SphereHit hit = sphere_hit(
                e.position, dir_out, std::numeric_limits<double>::max(),
                detector_, cfg_.scoring_radius_cm);
            if (hit.chord > 0.0) {
                const double mu_e = std::max(
                    -1.0, std::min(1.0, e.electron_dir.dot(dir_out)));
                const double p_sr = brems_2bs_pdf_sr(
                    e.electron_KE_keV, e.photon_keV, eff_z, mu_e);
                if (p_sr > 0.0) {
                    const double t_leg = env_.transmission(
                        e.position, e.position + dir_out * hit.entry,
                        e.photon_keV * 1e-3, walk.detector_trace);
                    const double w = weight * integration_weight * d_omega *
                                     p_sr * t_leg *
                                     attenuated_chord(hit.chord,
                                                      e.photon_keV) /
                                     scoring_volume_cm3_;
                    tally.add(e.photon_keV, dir_out.dot(det_axis_), w, cls,
                              FluxOrigin::Brems, PostCreationOrder::Zero);
                }
            }
        }

        // Walk the analog photon so its own scatter chain scores. Its
        // uncollided first flight is covered exactly once, by the emission
        // score above -- the walk scores only at collisions (see the
        // annihilation channel for the same partition).
        if (cfg_.walk_brems) {
            PendingPhoton b;
            b.pos = e.position;
            b.dir = e.photon_dir;
            b.energy = e.photon_keV;
            b.weight = weight;
            b.touched_slab = parent.touched_slab;
            b.first_slab = parent.first_slab;
            b.multi_slab = parent.multi_slab;
            b.slab_hits = parent.slab_hits;
            b.elastic_line = -1;   // continuum from birth
            b.origin = FluxOrigin::Brems;
            b.post_creation_interactions = 0;
            walk.stack.push_back(b);
        }
    }
}

void NextEventEstimator::run_history(double e0_keV, std::mt19937_64& rng,
                                     Walker& walk, FluenceTally& tally,
                                     EnvTransportDiagnostics& diagnostics) const
{
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    // A history is a small stack of photons: the primary, plus the two 511 keV
    // photons a pair-production vertex creates. Depth is bounded at one
    // generation because a 511 is below the 1022 keV pair threshold and can
    // never itself pair-produce, so at most 3 photons per history.
    const int line_e0 = 0;                       // ElasticE0 is always channel 0
    const int line_ann = (tally.n_lines() > 1) ? 1 : -1;

    walk.stack.clear();
    {
        // Explicit member assignment, NOT positional aggregate init: a field
        // insertion would shift a positional list (e.g. line_e0 landing on
        // slab_hits and elastic_line falling back to -1), silently emptying
        // the ElasticE0 line channel while every total stays correct.
        const Eigen::Vector3d primary_dir = sample_isotropic_direction(rng);
        PendingPhoton p;
        p.pos = source_;
        p.dir = primary_dir;
        p.energy = e0_keV;
        p.weight = 1.0;
        p.elastic_line = line_e0;
        p.post_creation_interactions = 0;

        if (cfg_.experimental_primary_sampling ==
            ExperimentalPrimarySampling::AntipodalPair) {
            // The stack is LIFO.  Push the antipode first so the +direction
            // branch consumes the same downstream RNG prefix as the analog A/B
            // run.  Both half-weight branches remain open until the caller's
            // single end_history(), which is essential for an exact S2.
            PendingPhoton antipode = p;
            antipode.dir = -primary_dir;
            antipode.weight = 0.5;
            walk.stack.push_back(antipode);
            p.weight = 0.5;
        }
        walk.stack.push_back(p);
    }

    // GEANT4 cannot represent an unbounded world, so a matched comparison has
    // to stop where the GDML box stops. Zero means unbounded (the default).
    const double whs = cfg_.world_half_size_cm;
    auto outside_world = [whs](const Eigen::Vector3d& p) {
        return whs > 0.0 && (std::abs(p.x()) > whs || std::abs(p.y()) > whs ||
                             std::abs(p.z()) > whs);
    };

    while (!walk.stack.empty()) {
        PendingPhoton ph = walk.stack.back();
        walk.stack.pop_back();

        int vtx = 0;
        for (; vtx < cfg_.max_vertices; ++vtx) {
            env_.trace(ph.pos, ph.dir, walk.segs, walk.transport_trace);
            if (walk.segs.empty()) break;          // straight out of the world

            const double target = -std::log(std::max(1e-300, 1.0 - uni(rng)));
            double accumulated = 0.0;
            const Material* hit_mat = nullptr;
            MacroscopicXS hit_xs{0.0, 0.0, 0.0, 0.0};
            int hit_slab = -1;
            double hit_t = 0.0;
            for (const EnvSegment& sg : walk.segs) {
                const MacroscopicXS segment_xs =
                    sg.material->macroscopic_xs(ph.energy * 1e-3);
                const double mu = segment_xs.mu_total();
                const double len = sg.t_end_cm - sg.t_start_cm;
                if (mu <= 0.0) continue;
                if (accumulated + mu * len >= target) {
                    hit_t = sg.t_start_cm + (target - accumulated) / mu;
                    hit_mat = sg.material;
                    hit_xs = segment_xs;
                    hit_slab = sg.slab_index;
                    break;
                }
                accumulated += mu * len;
            }
            if (hit_mat == nullptr) break;         // escaped without interacting

            const Eigen::Vector3d vertex = ph.pos + ph.dir * hit_t;
            // Beyond the emulated world the photon would already have been
            // killed at the boundary, so it neither scores nor continues.
            if (outside_world(vertex)) break;

            if (hit_slab >= 0) {
                ++ph.slab_hits;
                if (ph.touched_slab && hit_slab != ph.first_slab) ph.multi_slab = true;
                if (!ph.touched_slab) { ph.touched_slab = true; ph.first_slab = hit_slab; }
            }
            const ChainClass cls = chain_class(ph.slab_hits, ph.multi_slab);

            // SCORE FIRST, then continue the walk analog.
            score_vertex(vertex, ph.dir, ph.energy, ph.weight, *hit_mat, cls,
                         ph.elastic_line,
                         order_bucket(ph.post_creation_interactions + 1),
                         rng, tally, walk.detector_trace, ph.origin);

            const InteractionType it = select_interaction(hit_xs, uni(rng));
            diagnostics.add_interaction(diagnostic_counter(it));

            if (it == InteractionType::Photoelectric) {
                // The photoelectron carries essentially the full photon
                // energy (binding neglected: the air/soil K-edges are <= 5%
                // of any KE above the brems gate); let it radiate before the
                // photon is absorbed.
                if (cfg_.include_brems &&
                    ph.energy >= cfg_.brems_electron_floor_keV) {
                    const Eigen::Vector3d e_dir =
                        sample_photoelectron_direction(ph.dir, ph.energy, rng);
                    score_brems_electron(vertex, e_dir, ph.energy, ph.weight,
                                         *hit_mat, /*is_positron=*/false, ph,
                                         cls, rng, walk, tally, diagnostics);
                }
                break;   // absorbed
            }

            if (it == InteractionType::PairProduction) {
                // One flat split shared by everything the pair does: the e-
                // and e+ that radiate brems and the e+ whose range displaces
                // the 511 birth point are the same particles.
                const double avail =
                    std::max(0.0, ph.energy - kPairThresholdKeV);
                const double ke_electron = uni(rng) * avail;
                const double ke_positron = avail - ke_electron;
                // Both pair electrons radiate (forward-peaked at these
                // energies, so their initial direction is the photon's).
                score_brems_electron(vertex, ph.dir, ke_electron, ph.weight,
                                     *hit_mat, /*is_positron=*/false, ph, cls,
                                     rng, walk, tally, diagnostics);
                score_brems_electron(vertex, ph.dir, ke_positron, ph.weight,
                                     *hit_mat, /*is_positron=*/true, ph, cls,
                                     rng, walk, tally, diagnostics);
                if (cfg_.include_annihilation &&
                    ph.energy >= kPairThresholdKeV && line_ann >= 0) {
                    // The 511s are born where the POSITRON STOPS, not at the
                    // pair-production vertex -- metres away in air.
                    const Eigen::Vector3d p_ann =
                        cfg_.positron_range
                            ? annihilation_point(env_, vertex, ph.dir,
                                                 ke_positron, walk.segs,
                                                 walk.transport_trace)
                            : vertex;
                    // Uncollided arrival of BOTH photons: a discrete line.
                    score_isotropic(p_ann, kAnnihilationKeV, 2.0 * ph.weight, cls,
                                    line_ann, rng, tally, walk.detector_trace,
                                    FluxOrigin::Annihilation);
                    if (cfg_.walk_annihilation) {
                        // Emit antipodally: physically right, and the pairing is
                        // free antithetic variance reduction. The MEAN needs only
                        // each photon's marginal isotropy, so this choice cannot
                        // bias a linear tally either way.
                        const Eigen::Vector3d w_dir = sample_isotropic_direction(rng);
                        for (int sgn = 0; sgn < 2; ++sgn) {
                            PendingPhoton a;
                            a.pos = p_ann;
                            a.dir = (sgn == 0) ? w_dir : Eigen::Vector3d(-w_dir);
                            a.energy = kAnnihilationKeV;
                            a.weight = ph.weight;
                            a.touched_slab = ph.touched_slab;
                            a.first_slab = ph.first_slab;
                            a.multi_slab = ph.multi_slab;
                            a.slab_hits = ph.slab_hits;
                            a.elastic_line = line_ann;   // still purely elastic
                            a.origin = FluxOrigin::Annihilation;
                            a.post_creation_interactions = 0;
                            walk.stack.push_back(a);
                        }
                    }
                }
                break;   // the parent photon is consumed by pair production
            }

            ph.pos = vertex;
            if (it == InteractionType::Compton) {
                const int Z = hit_mat->select_element(ph.energy * 1e-3, 1, uni(rng));
                const double mu = sample_compton_cos_theta(
                    ph.energy, cfg_.free_electron_compton ? 0 : Z, rng);
                const Eigen::Vector3d dir_in = ph.dir;
                const double e_in = ph.energy;
                ComptonScatterResult r =
                    finish_compton_at_angle(e_in, dir_in, mu, rng, Z, cfg_.doppler);
                // The recoil electron's KE is the dominant brems source.
                if (cfg_.include_brems &&
                    r.deposited_energy_keV >= cfg_.brems_electron_floor_keV) {
                    const Eigen::Vector3d e_dir = compton_recoil_direction(
                        dir_in, r.new_direction, e_in, r.scattered_energy_keV);
                    score_brems_electron(vertex, e_dir, r.deposited_energy_keV,
                                         ph.weight, *hit_mat,
                                         /*is_positron=*/false, ph, cls, rng,
                                         walk, tally, diagnostics);
                }
                ph.energy = r.scattered_energy_keV;
                ph.dir = r.new_direction;
                ph.elastic_line = -1;      // inelastic: the chain is continuum now
            } else {                        // Rayleigh: elastic, line flag survives
                const int Z = hit_mat->select_element(ph.energy * 1e-3, 2, uni(rng));
                const double mu = sample_rayleigh_cos_theta(Z, ph.energy, rng);
                ph.dir = rotate_direction(ph.dir, mu, kTwoPi * uni(rng));
            }
            ph.post_creation_interactions =
                std::min(2, ph.post_creation_interactions + 1);

            if (ph.energy < cfg_.e_floor_keV) break;

            if (hit_slab >= 0) {
                const EnvSlab& sl = env_.slabs()[hit_slab];
                const double depth = -(ph.pos - sl.point).dot(sl.normal);
                if (depth > 0.0) {
                    const double tau = sl.material.mu_total(ph.energy * 1e-3) * depth;
                    if (tau > cfg_.roulette_tau) {
                        if (uni(rng) < 0.5) break;
                        ph.weight *= 2.0;
                    }
                }
            }
            if (ph.weight < cfg_.weight_floor) {
                if (uni(rng) < 0.5) break;
                ph.weight *= 2.0;
            }
        }
        if (vtx == cfg_.max_vertices)
            ++diagnostics.max_vertex_terminations;
    }
}

void NextEventEstimator::run_brute(double e0_keV, uint64_t n_histories,
                                   uint64_t seed, EnvScatterResult& out,
                                   int logical_replica) const
{
    if (logical_replica < -1 ||
        logical_replica >= FluenceTally::kNumLogicalReplicas)
        throw std::out_of_range(
            "NextEventEstimator::run_brute: logical replica must be -1 or in "
            "[0, 31]");
    if (out.tally.n_lines() < 1 ||
        out.tally.line(0).energy_keV != e0_keV)
        throw std::invalid_argument(
            "NextEventEstimator::run_brute: EnvScatterResult was built for a "
            "different primary energy than the one being run");
    if (cfg_.experimental_primary_sampling !=
            ExperimentalPrimarySampling::Analog ||
        cfg_.experimental_cone_sampling !=
            ExperimentalConeSampling::Single)
        throw std::invalid_argument(
            "NextEventEstimator::run_brute: experimental next-event variance "
            "reduction policies are not implemented for analog truth runs");
    std::seed_seq seq{static_cast<uint32_t>(seed),
                      static_cast<uint32_t>(seed >> 32),
                      static_cast<uint32_t>(seed * 2654435761ULL)};
    std::mt19937_64 rng(seq);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    const double R = cfg_.scoring_radius_cm;
    const double volume = scoring_volume_cm3_;

    // The annihilation pair is walked here too, so this channel is gated: the
    // estimator scores the 511s explicitly, and a brute walk that terminated
    // at pair production would differ from it by exactly that term with
    // nothing to catch it.
    struct BrutePhoton {
        Eigen::Vector3d pos{0, 0, 0};
        Eigen::Vector3d dir{0, 0, 1};
        double energy = 0.0;
        bool touched_slab = false;
        int first_slab = -1;
        bool multi_slab = false;
        int slab_hits = 0;
        FluxOrigin origin = FluxOrigin::Primary;
        int post_creation_interactions = 0;
        int elastic_line = -1;
    };

    std::vector<EnvSegment> segs;
    EnvTraceScratch trace_scratch;
    std::vector<BrutePhoton> pending;
    const int line_e0 = 0;
    const int line_ann = (out.tally.n_lines() > 1) ? 1 : -1;
    for (uint64_t h = 0; h < n_histories; ++h) {
      pending.clear();
      {
        BrutePhoton p;
        p.pos = source_;
        p.dir = sample_isotropic_direction(rng);
        p.energy = e0_keV;
        p.elastic_line = line_e0;
        pending.push_back(p);
      }
      while (!pending.empty()) {
        BrutePhoton cur = pending.back();
        pending.pop_back();
        Eigen::Vector3d pos = cur.pos;
        Eigen::Vector3d dir = cur.dir;
        double energy = cur.energy;
        // Annihilation and brems photons are scattered products from birth, so
        // their very FIRST flight already counts -- that flight is the
        // counterpart of the estimator's uncollided emission score.
        bool scattered = cur.origin != FluxOrigin::Primary ||
                         cur.post_creation_interactions > 0;
        bool touched_slab = cur.touched_slab;
        int first_slab = cur.first_slab;
        bool multi_slab = cur.multi_slab;
        int slab_hits = cur.slab_hits;
        const FluxOrigin origin = cur.origin;
        int post_creation_interactions = cur.post_creation_interactions;
        int elastic_line = cur.elastic_line;

        // Grow the same brems channel the estimator scores, so the internal
        // gate stays honest (the un-walked-511 lesson): each electron the
        // walk creates radiates via the SAME emission-recording sampler, and
        // its photons join the analog stack as scattered-from-birth.
        auto emit_brems = [&](const Material& mat, const Eigen::Vector3d& e_dir,
                              double ke, bool is_positron,
                              const Eigen::Vector3d& vtx) {
            if (!cfg_.include_brems || !cfg_.walk_brems ||
                ke < cfg_.brems_electron_floor_keV)
                return;
            const std::vector<BremsEmission> ems =
                ElectronCsda::instance().sample_brems_emissions_in_material(
                    mat, vtx, e_dir, ke, rng, is_positron);
            out.diagnostics.add_interaction(
                EnvTransportCounter::BremsEmission,
                static_cast<std::uint64_t>(ems.size()));
            for (const BremsEmission& e : ems) {
                if (e.photon_keV < cfg_.e_floor_keV) continue;
                BrutePhoton b;
                b.pos = e.position;
                b.dir = e.photon_dir;
                b.energy = e.photon_keV;
                b.touched_slab = touched_slab;
                b.first_slab = first_slab;
                b.multi_slab = multi_slab;
                b.slab_hits = slab_hits;
                b.origin = FluxOrigin::Brems;
                b.post_creation_interactions = 0;
                b.elastic_line = -1;
                pending.push_back(b);
            }
        };

        int vtx = 0;
        for (; vtx < cfg_.max_vertices; ++vtx) {
            env_.trace(pos, dir, segs, trace_scratch);
            if (segs.empty()) break;

            // Sample the flight; if it escapes, the flight runs to the end of
            // the traced path (and the photon then leaves the world).
            const double target = -std::log(std::max(1e-300, 1.0 - uni(rng)));
            double accumulated = 0.0;
            double s_flight = segs.back().t_end_cm;
            const Material* hit_mat = nullptr;
            MacroscopicXS hit_xs{0.0, 0.0, 0.0, 0.0};
            int hit_slab = -1;
            for (const EnvSegment& s : segs) {
                const MacroscopicXS segment_xs =
                    s.material->macroscopic_xs(energy * 1e-3);
                const double mu = segment_xs.mu_total();
                const double len = s.t_end_cm - s.t_start_cm;
                if (mu > 0.0 && accumulated + mu * len >= target) {
                    s_flight = s.t_start_cm + (target - accumulated) / mu;
                    hit_mat = s.material;
                    hit_xs = segment_xs;
                    hit_slab = s.slab_index;
                    break;
                }
                if (mu > 0.0) accumulated += mu * len;
            }

            // Tally the track length inside the scoring sphere for this flight.
            // Only post-scatter flights count, so this is comparable to the
            // estimator's scattered-fluence tally.
            if (scattered) {
                const double L = chord_in_sphere(pos, dir, s_flight, detector_, R);
                if (L > 0.0) {
                    const ChainClass cls = chain_class(slab_hits, multi_slab);
                    if (elastic_line >= 0)
                        out.tally.add_line(
                            elastic_line, dir.dot(det_axis_), L / volume, cls,
                            origin, order_bucket(post_creation_interactions));
                    else
                        out.tally.add(
                            energy, dir.dot(det_axis_), L / volume, cls, origin,
                            order_bucket(post_creation_interactions));
                }
            }

            if (hit_mat == nullptr) break;  // escaped

            pos = pos + dir * s_flight;
            // Same world bound as the estimator, so the two stay comparable
            // when the emulation is switched on.
            if (cfg_.world_half_size_cm > 0.0 &&
                (std::abs(pos.x()) > cfg_.world_half_size_cm ||
                 std::abs(pos.y()) > cfg_.world_half_size_cm ||
                 std::abs(pos.z()) > cfg_.world_half_size_cm)) break;
            if (hit_slab >= 0) {
                ++slab_hits;
                if (touched_slab && hit_slab != first_slab) multi_slab = true;
                if (!touched_slab) { touched_slab = true; first_slab = hit_slab; }
            }

            const InteractionType it = select_interaction(hit_xs, uni(rng));
            out.diagnostics.add_interaction(diagnostic_counter(it));
            if (it == InteractionType::Photoelectric) {
                if (cfg_.include_brems && cfg_.walk_brems &&
                    energy >= cfg_.brems_electron_floor_keV) {
                    const Eigen::Vector3d e_dir =
                        sample_photoelectron_direction(dir, energy, rng);
                    emit_brems(*hit_mat, e_dir, energy,
                               /*is_positron=*/false, pos);
                }
                break;
            }
            if (it == InteractionType::PairProduction) {
                // One flat split shared by the radiating pair electrons and
                // the positron's range displacement (same draw as run()).
                const double avail = std::max(0.0, energy - kPairThresholdKeV);
                const double ke_electron = uni(rng) * avail;
                const double ke_positron = avail - ke_electron;
                emit_brems(*hit_mat, dir, ke_electron,
                           /*is_positron=*/false, pos);
                emit_brems(*hit_mat, dir, ke_positron,
                           /*is_positron=*/true, pos);
                // Emit the annihilation pair so this channel is actually gated.
                // The estimator scores the 511s' UNCOLLIDED arrival explicitly,
                // so here the pair must be walked with scattered = true from
                // birth -- their very first flight has to contribute its real
                // track length, or the two sides differ by exactly that term.
                if (cfg_.include_annihilation && cfg_.walk_annihilation &&
                    energy >= kPairThresholdKeV) {
                    const Eigen::Vector3d p_ann =
                        cfg_.positron_range
                            ? annihilation_point(env_, pos, dir, ke_positron,
                                                 segs, trace_scratch)
                            : pos;
                    const Eigen::Vector3d w_dir = sample_isotropic_direction(rng);
                    for (int sgn = 0; sgn < 2; ++sgn) {
                        BrutePhoton a;
                        a.pos = p_ann;
                        a.dir = (sgn == 0) ? w_dir : Eigen::Vector3d(-w_dir);
                        a.energy = kAnnihilationKeV;
                        a.touched_slab = touched_slab;
                        a.first_slab = first_slab;
                        a.multi_slab = multi_slab;
                        a.slab_hits = slab_hits;
                        a.origin = FluxOrigin::Annihilation;
                        a.post_creation_interactions = 0;
                        a.elastic_line = line_ann;
                        pending.push_back(a);
                    }
                }
                break;
            }
            if (it == InteractionType::Compton) {
                const int Z = hit_mat->select_element(energy * 1e-3, 1, uni(rng));
                const double mu = sample_compton_cos_theta(
                    energy, cfg_.free_electron_compton ? 0 : Z, rng);
                const Eigen::Vector3d dir_in = dir;
                const double e_in = energy;
                ComptonScatterResult r =
                    finish_compton_at_angle(e_in, dir_in, mu, rng, Z, cfg_.doppler);
                if (cfg_.include_brems && cfg_.walk_brems &&
                    r.deposited_energy_keV >= cfg_.brems_electron_floor_keV) {
                    const Eigen::Vector3d e_dir = compton_recoil_direction(
                        dir_in, r.new_direction, e_in, r.scattered_energy_keV);
                    emit_brems(*hit_mat, e_dir, r.deposited_energy_keV,
                               /*is_positron=*/false, pos);
                }
                energy = r.scattered_energy_keV;
                dir = r.new_direction;
                elastic_line = -1;
            } else {
                const int Z = hit_mat->select_element(energy * 1e-3, 2, uni(rng));
                const double mu = sample_rayleigh_cos_theta(Z, energy, rng);
                dir = rotate_direction(dir, mu, kTwoPi * uni(rng));
            }
            post_creation_interactions =
                std::min(2, post_creation_interactions + 1);
            scattered = true;
            if (energy < cfg_.e_floor_keV) break;
        }
        if (vtx == cfg_.max_vertices)
            ++out.diagnostics.max_vertex_terminations;
      }
      if (logical_replica >= 0) out.tally.end_history(logical_replica);
      else                       out.tally.end_history();
    }

    out.histories += n_histories;
    const double d = (detector_ - source_).norm();
    out.source_distance_cm = d;
    out.direct_transmission = direct_transmission(e0_keV);
    out.direct_fluence = (d > 0.0)
        ? out.direct_transmission / (kFourPi * d * d)
        : 0.0;
}

void NextEventEstimator::run(double e0_keV, uint64_t n_histories, uint64_t seed,
                             EnvScatterResult& out, int logical_replica) const
{
    if (logical_replica < -1 ||
        logical_replica >= FluenceTally::kNumLogicalReplicas)
        throw std::out_of_range(
            "NextEventEstimator::run: logical replica must be -1 or in [0, 31]");
    if (out.tally.n_lines() < 1 ||
        out.tally.line(0).energy_keV != e0_keV) {
        // A result carries line channels keyed to ONE primary energy; reusing it
        // at another energy would silently mis-file the elastic line.
        throw std::invalid_argument(
            "NextEventEstimator::run: EnvScatterResult was built for a different "
            "primary energy than the one being run");
    }
    if (cfg_.primary_mode != PrimaryMode::Analog) {
        // Fail loudly rather than silently running Analog: a caller asking for
        // importance-sampled vertices would otherwise get correct-but-slow
        // results and never learn the mode was ignored.
        throw std::invalid_argument(
            "NextEventEstimator: PrimaryMode::VertexSampled is not implemented; "
            "only PrimaryMode::Analog is available");
    }

    std::seed_seq seq{static_cast<uint32_t>(seed),
                      static_cast<uint32_t>(seed >> 32),
                      static_cast<uint32_t>(seed * 2654435761ULL)};
    std::mt19937_64 rng(seq);

    Walker walk;
    for (uint64_t i = 0; i < n_histories; ++i) {
        run_history(e0_keV, rng, walk, out.tally, out.diagnostics);
        if (logical_replica >= 0) out.tally.end_history(logical_replica);
        else                       out.tally.end_history();
    }
    out.histories += n_histories;

    const double d = (detector_ - source_).norm();
    out.source_distance_cm = d;
    out.direct_transmission = direct_transmission(e0_keV);
    out.direct_fluence = (d > 0.0)
        ? out.direct_transmission / (kFourPi * d * d)
        : 0.0;
}

} // namespace ceelo
