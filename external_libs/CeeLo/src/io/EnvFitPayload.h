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

/// @file EnvFitPayload.h
/// @brief Compact, positive runtime representation of fitted environment flux.
///
/// This is deliberately not an MC archive.  It contains only fitted spline
/// coefficients, normalized joint-shape archetypes, validity domains, and
/// validation metadata.  Exact-energy atoms remain separate from continuum
/// cells at evaluation time.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace ceelo {

constexpr std::uint32_t kEnvFitPayloadFormatVersion = 1;
constexpr double kEnvFitElectronRestEnergyKeV = 510.998950;

enum class EnvFitAxisKind : std::uint32_t {
    Primary = 1,
    Annihilation = 2,
    Brems = 3,
    LineE0 = 4,
    Line511 = 5
};

/// The fitted equation is dimensionless s_hat.  Raw archive fluence is divided
/// by the independently known analytic direct fluence before fitting, and the
/// runtime query supplies that direct fluence to recover incident cm^-2.
enum class EnvFitPayloadNormalization : std::uint32_t {
    SHatToDirectFluence = 1
};

/// One positive dimensionless-s_hat model
/// Y(q) sum_r softmax(g_r(q)) B_r(u,mu).
///
/// The q fields use a four-dimensional tensor product of clamped uniform
/// cubic B-splines.  Levels 1--4 are nested dyadic refinements with 1, 2, 4,
/// and 8 intervals per coordinate, hence 4^4, 5^4, 7^4, and 11^4 features.
/// Coefficients are row-major with the fourth q coordinate fastest.  Archetypes are stored as
/// nonnegative coefficients of individually normalized tensor-product cubic
/// B-spline PDFs.  Consequently every coefficient row is a convex combination
/// and every archetype has exactly unit integral.
struct EnvFitRuntimeModel {
    std::string name;
    EnvFitAxisKind axis_kind = EnvFitAxisKind::Primary;
    double e_cut_keV = 10.0;
    bool has_energy_threshold = false;
    double energy_threshold_keV = 0.0;
    std::uint32_t rank = 0;
    std::uint32_t level = 0;

    /// Clamped cubic knots on registered energy u in [0,1].  Empty for lines.
    std::vector<double> energy_knots;
    /// Clamped cubic knots on t=(mu+1)/2 in [0,1].
    std::vector<double> mu_knots;
    std::array<double, 4> q_min{};
    std::array<double, 4> q_span{};

    /// [n_features]
    std::vector<double> log_yield_coefficients;
    /// [n_features][rank], rank fastest.
    std::vector<double> logit_coefficients;
    /// [rank][n_energy_basis][n_mu_basis], mu fastest.  Line models omit the
    /// energy dimension and store [rank][n_mu_basis].  Each rank row sums one.
    std::vector<double> archetype_coefficients;
};

/// Explicit expected runtime inventory.  The codec rejects any fitted model
/// that is absent here, as well as any inventory entry without a corresponding
/// fitted model.  Names and axes are the canonical physics taxonomy.
struct EnvFitModelInventoryEntry {
    std::string name;
    EnvFitAxisKind axis_kind = EnvFitAxisKind::Primary;
    bool line = false;
};

struct EnvFitPayload {
    std::string config_sha256;
    std::string material_name;
    std::string geometry_tier;
    EnvFitPayloadNormalization normalization =
        EnvFitPayloadNormalization::SHatToDirectFluence;
    /// Inclusive campaign bounds in (E0_keV, distance_m, hsrc_m, hdet_m).
    /// These are the authoritative physical domain; q boxes only delimit each
    /// fitted spline.  The evaluator additionally enforces |hsrc-hdet| < d.
    std::array<double, 4> physical_domain_min{};
    std::array<double, 4> physical_domain_max{};
    /// True only after untouched-final acceptance.  File loading rejects false
    /// by default; validation tooling must opt in explicitly.
    bool accepted_for_publication = false;
    /// Campaign physics setting.  Every brems model must use exactly this
    /// source-energy activation threshold (zero means active throughout).
    double brems_activation_threshold_keV = 0.0;
    /// SHA-256 over configuration, domain, inventory, and fitted model state.
    /// Validation JSON and the publication flag are deliberately excluded.
    /// Newly constructed payloads may leave this empty; save computes it.
    std::string model_content_sha256;
    std::vector<EnvFitModelInventoryEntry> model_inventory;
    /// Canonical JSON generated by the validation pipeline.  The evaluator
    /// treats it as opaque provenance, never as executable model state.
    std::string validation_metadata_json;
    std::vector<EnvFitRuntimeModel> continuum_models;
    std::vector<EnvFitRuntimeModel> line_models;
};

struct EnvFitRuntimeQuery {
    double e0_keV = 0.0;
    /// Three-dimensional source-detector separation.
    double distance_m = 0.0;
    double source_height_m = 0.0;
    double detector_height_m = 0.0;
    /// Analytic T(E0)/(4*pi*d_cm^2), per source photon.  Required at runtime;
    /// deliberately not predicted by, or baked into, the fitted model.
    double direct_fluence_per_cm2 = 0.0;
};

struct EnvFitEvaluatedComponent {
    std::string name;
    /// Dimensionless [requested energy][requested mu] s_hat masses.
    std::vector<double> s_hat_cell_masses;
    /// Incident fluence per source photon [cm^-2] = s_hat * direct fluence.
    std::vector<double> incident_cell_masses_per_cm2;
};

struct EnvFitEvaluatedLine {
    std::string name;
    double energy_keV = 0.0;
    std::vector<double> s_hat_mu_bin_masses;
    std::vector<double> incident_mu_bin_masses_per_cm2;
};

struct EnvFitEvaluation {
    std::uint32_t n_energy = 0;
    std::uint32_t n_mu = 0;
    double direct_fluence_per_cm2 = 0.0;
    std::vector<double> continuum_s_hat_cell_masses;
    std::vector<double> incident_continuum_cell_masses_per_cm2;
    std::vector<EnvFitEvaluatedComponent> components;
    std::vector<EnvFitEvaluatedLine> lines;
};

/// Transform physical coordinates to
/// (log E0, log rho, atan(hs/rho), atan(hd/rho)).
bool env_fit_transform_query(const EnvFitRuntimeQuery& query,
                             std::array<double, 4>& q,
                             std::string* error = nullptr);

/// Validate dimensions, domains, finiteness, exact support taxonomy, and
/// normalized/nonnegative archetypes.
bool validate_env_fit_payload(const EnvFitPayload& payload,
                              std::string* error = nullptr);

/// Stable candidate identity used to bind detector-fold and acceptance
/// reports.  It excludes validation metadata and accepted_for_publication.
std::string env_fit_payload_model_content_sha256(const EnvFitPayload& payload);

/// Evaluate dimensionless s_hat and incident continuum (s_hat times the
/// query's analytic direct fluence) as bin-integrated masses on arbitrary
/// increasing requested axes.  Requested axes may extend beyond or select a
/// subset of physical support; bins outside support receive exactly zero.
/// Exact lines are returned separately and are never smeared into energy bins.
/// Extrapolation in q or physical coordinates is rejected rather than clipped.
bool evaluate_env_fit_payload(const EnvFitPayload& payload,
                              const EnvFitRuntimeQuery& query,
                              const std::vector<double>& energy_edges_keV,
                              const std::vector<double>& mu_edges,
                              EnvFitEvaluation& result,
                              std::string* error = nullptr);

/// Dependency-free, checksummed CEEFIT-PAYLOAD-v1 binary codec.  All integer
/// fields and IEEE-754 coefficients are little-endian.  The 64-byte header is
/// magic[16], version:u32, endian-tag:u32, payload-size:u64, SHA256[32].  The
/// body is config/material/tier strings, normalization:u32, publication:u8,
/// brems threshold:f64, model-content SHA string, 4 physical minima, 4 physical
/// maxima, explicit model inventory, validation JSON, then continuum and line
/// model vectors.  A string/vector is
/// prefixed by a u64 element count.  Each model stores name, axis:u32, Ecut,
/// threshold flag/value, rank/level, energy/mu knots, q domain, yield/logit
/// vectors, and archetype coefficients in that order.  save is atomic within
/// the destination directory.
bool save_env_fit_payload(const EnvFitPayload& payload,
                          const std::string& path,
                          std::string* error = nullptr);
bool load_env_fit_payload(EnvFitPayload& payload,
                          const std::string& path,
                          std::string* error = nullptr,
                          bool allow_candidate = false);

/// Reusable evaluator that validates a payload once.  Prefer this over the
/// free convenience function when evaluating many geometries; it avoids an
/// O(payload-size) coefficient audit on every query.
class EnvFitPayloadEvaluator {
public:
    bool set_payload(EnvFitPayload payload, std::string* error = nullptr,
                     bool allow_candidate = false);
    bool load(const std::string& path, std::string* error = nullptr,
              bool allow_candidate = false);
    bool ready() const { return ready_; }
    const EnvFitPayload& payload() const { return payload_; }

    bool evaluate(const EnvFitRuntimeQuery& query,
                  const std::vector<double>& energy_edges_keV,
                  const std::vector<double>& mu_edges,
                  EnvFitEvaluation& result,
                  std::string* error = nullptr) const;

private:
    EnvFitPayload payload_;
    bool ready_ = false;
};

}  // namespace ceelo
