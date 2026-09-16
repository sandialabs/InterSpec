#ifndef CEELO_IO_DETECTOR_RESPONSE_H
#define CEELO_IO_DETECTOR_RESPONSE_H
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

/// @file DetectorResponse.h
/// @brief The storable, MC-parameterized detector response object.
///
/// Implements the response function of the July 2026 parameterization
/// campaign (its response-function specification
/// in the dev repo -- "the spec" below):
///
///     eps_fep(E, pos) = k(E) * eta(E, theta[, phi]) * N(E, Omega, theta_eff)
///                       * K(E, pos)                                  (Eq. 3)
///
///   K       geometry kernel (io/ResponseKernel.h) evaluated at query time
///           over the exact detector geometry -- never interpolated. K uses
///           the mu-tables CAPTURED AT GENERATION (MuTable below), so eta*K
///           cannot drift if the library's cross-section data is later
///           regenerated; the stored object is self-contained.
///   eta     small residual table eps/K on greedy-placed (E x cos-theta
///           [x phi]) nodes, K-edge segmented log-log PCHIP (Eq. spec 3.1)
///   N       near-field multiplier: tabulated ln N on a (cos-theta x d) grid
///           per E-node, PCHIP-interpolated, gated by measured breakpoints
///           (spec Eq. 6, tabulated form)
///   k       grounding ratio to measured peak efficiencies: piecewise-linear
///           ln k(ln E) hat basis + full covariance (spec Eqs. 7a-8b)
///
/// eps_total is tiered (spec Eq. 4): bare crystal -> kernel-exact
/// K_{mu-mu_RS}; canned scintillator -> x b(E); thick-dead-layer HPGe ->
/// k(E)*eta_tot(E,theta)*K -- the tier is chosen at generation time.
///
/// Every query returns {value, sigma, flag}; sigma combines interpolated
/// node-MC variance, coverage-tuned per-regime model floors, grounding
/// covariance and the empirical grounding-transfer inflation (spec sec 5).
///
/// XML: to_xml_string()/from_xml_string() -- ONE codec shared by the
/// generator and by InterSpec (which embeds the fragment in its DRF XML).

#include "geometry/Geometry.h"
#include "io/Pchip.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ceelo {

// ---------------------------------------------------------------------------
// Query result
// ---------------------------------------------------------------------------

/// Provenance flag returned with every efficiency query (spec sec 4.5).
/// Priority order when several apply: NeedsMc > Shadowed >
/// NearFieldUnmodeled > OutOfRangeClamped > Ok.
enum class ResponseFlag : uint8_t {
    Ok = 0,               ///< inside validated ranges
    OutOfRangeClamped,    ///< E or theta outside node range; value clamped
    NearFieldUnmodeled,   ///< d below measured breakpoint, no near model
                          ///< (far-field profile queried close in)
    Shadowed,             ///< collimator hole fraction s in [0.05, 0.3):
                          ///< scatter/leak-dominated, sigma inflated
    NeedsMc               ///< refuse-grade query (s < 0.05, behind-plane, ...):
                          ///< value is a guess; offer a bespoke MC run
};

const char* to_string(ResponseFlag f);

/// {value, 1-sigma absolute uncertainty, provenance flag}.  `sigma` is the total; `sigma_model`
/// is the part of it that comes from the ad hoc model envelopes alone (see model_sigma and
/// frac_covariance) - the remainder is data-derived (MC node statistics, anchor/grounding fit
/// covariance).  sigma_model <= sigma always.
struct EffResult {
    double value = 0.0;
    double sigma = 0.0;
    ResponseFlag flag = ResponseFlag::Ok;
    double sigma_model = 0.0;
};

// ---------------------------------------------------------------------------
// Shielded-eps_tot build-up seam (Stage E3 A1)
// ---------------------------------------------------------------------------

/// Effective-shield descriptor handed to a build-up model.  The uncollided-only
/// eps_total under shielding misses in-peak (and, for eps_tot, degraded) scatter
/// build-up by 10-49% (report sec 3.3), which matters for cascade summing-out.
/// This is the interface used by the cascade-summing path: the host application
/// (e.g. InterSpec, via computeEffectiveShielding) fills in the effective areal
/// density, effective atomic number and hydrogen fraction, and supplies the
/// build-up model.
struct ShieldContext {
    double areal_density_g_cm2 = 0.0;  ///< effective shield areal density
    double eff_atomic_number = 0.0;    ///< effective Z of the shield stack
    double hydrogen_frac = 0.0;        ///< H mass fraction (scatter softness)
    int geometry_class = 0;            ///< caller's geometry tag (point/slab/...)
};

/// Build-up factor B(E, shield) >= 1 = (uncollided + in-peak scatter) /
/// uncollided.  Applied MULTIPLICATIVELY to eps_total ONLY when a caller passes
/// a ShieldContext (see eps_total_* overloads); RATIO use only -- never an
/// absolute scatter estimate (GADRAS under-predicts detector-collected scatter
/// 1.7-2x).  NOT serialized (default null => zero behavior change for every
/// existing caller).  CeeLo ships the seam empty; InterSpec installs the model.
using BuildupModel = std::function<double(double E_keV, const ShieldContext&)>;

// ---------------------------------------------------------------------------
// Model-envelope uncertainty constants
// ---------------------------------------------------------------------------

/// EVERY ad hoc model-uncertainty constant CeeLo applies at query time lives here, and nowhere
/// else.  They are ENVELOPES - "where the model is not measured it may be wrong by about this
/// much" - not measurements, and DetectorResponse::SigmaBudget keeps them apart from the
/// data-derived terms (MC node sigma, anchor/grounding fit covariance) so a host can tell the two
/// kinds apart (EffResult::sigma_model, frac_covariance's `model_part`).  Each enters a query as a
/// fully-correlated common mode across energies at the query geometry.
///
/// Provenance, per constant.  Anything not listed as measured below still dates from the original
/// CeeLo import (InterSpec commit b26e5a67), whose S1/S7 campaign notes and "spec sec 4" grounding
/// table are not in this repository - those remain placeholders, and say so individually.
///
/// The 2026-09 corpus study (harness: InterSpec's test_CeeLoDrfIntegration, the `--envelope-*`
/// flags on `envelope_corpus_measure`; working notes and raw per-probe CSVs under
/// scratch/20260913_envelope_study/) measured 33 detectors - 4 CeeLo presets, 9 ANGLE imports and
/// 20 GADRAS Detector.dat imports, spanning NaI / HPGe / LaBr3 / CZT and cylinder / box - against
/// fresh, never-fitted MC on a stratified (E, d, theta) grid, generation and probes both at
/// node_fep_precision = 0.003.
///
/// WHICH STATISTIC.  These constants are consumed by a chi-square (InterSpec's
/// compute_efficiency_whitening), so each measured value is the floor that makes the **RMS pull
/// about MC truth equal 1** over every far-field probe in the corpus, with the response's own
/// declared node sigma already in the denominator.  That is deliberately NOT the same as a 68%
/// interval: at these values 75% (FEP) / 84% (total) of probes fall inside one sigma.  The
/// residual is heavy-tailed, so the two calibrations differ by ~1.6x and cannot both be had.  If
/// you re-derive, say which one you targeted.
///
/// NOT INCLUDED, deliberately: CeeLo-vs-GEANT4 model error.  Measured at <= 0.48% (FEP) and
/// <= 0.23% (total) on the bare/lightly-shielded G4 configs {1,2,3,5,6,25,26}, with no new MC, from
/// the committed references in tests/data/{geant4,ceelo}_reference (see
/// scratch/20260913_envelope_study/mc_vs_geant4.py).  It is a BOUND rather than a measurement
/// (the two codes' own counting noise is ~73% of the observed spread), it is code-vs-code rather
/// than against data, and its applicability to an arbitrary user detector is not established.  A
/// careful comparison with measured data is what would justify folding it in.  Do not add it on
/// the strength of the number above.
namespace model_sigma {
    /// theta > 90 degrees: no nodes behind the face plane; the value is a clamped guess.
    constexpr double behind_plane = 0.30;
    /// A far-field profile (no NearFieldModel) queried inside the near-field gate: the
    /// kernel-only near-field error (S1).
    constexpr double near_unmodeled = 0.05;
    /// Collimator shadow, by transmitted hole fraction s: below `shadow_refuse_s` the query is
    /// refuse-grade (sigma ~100%); up to `shadow_ramp_s` the sigma ramps linearly from
    /// `shadow_ramp_max` down to 0.
    constexpr double shadow_refuse = 1.0;
    constexpr double shadow_refuse_s = 0.05;
    constexpr double shadow_ramp_s = 0.30;
    constexpr double shadow_ramp_max = 0.5;
    /// Model-form floor on the (ratio-only) build-up correction of a shielded eps_total, added
    /// in quadrature when a ShieldContext is supplied.
    constexpr double buildup_floor = 0.10;
    /// FAR-FIELD PEAK FLOOR - MEASURED, 2026-09 corpus (33 detectors, far field, p = 0.003).
    /// Calibrated to RMS pull 1 over 5508 probes; achieved coverage 74.8% within one sigma.
    /// One constant serves every detector class: solving per class gives 0.475% excluding CZT and
    /// 0.356% for CZT alone, so the class spread does not justify a split for the peak efficiency
    /// (it does for the total - see tot_far_floor).  The previous 0.014 over-covered by 2.8x: at
    /// that value 97.7% of probes sat inside one sigma and the RMS pull was 0.44.
    /// Landed together with the re-derived transfer envelope, which is the only way it
    /// works: on its own this change made `act_fit_pulls_calibrated` WORSE off axis
    /// (RMS pull 1.23 -> 1.46 at 30 degrees against a gate of 1.3), because the oversized
    /// floor had been compensating for an off-axis term ~4x too small.  Two errors were
    /// cancelling; fixing either alone exposes the other.
    constexpr double fep_far_floor = 0.005;
    /// NEAR-FIELD PEAK FLOOR - NOT re-derived.  The 2026-09 study measured the FAR stratum only;
    /// no near-field probe bank was run, so this keeps its original import value and its original
    /// lack of provenance.  Deriving it needs the near stratum plus the d/a step test that would
    /// say whether `near_regime_a` marks a real step at all.
    constexpr double fep_near_floor = 0.023;
    /// FAR-FIELD TOTAL FLOOR - MEASURED, 2026-09 corpus, EXCLUDING CdTe-class crystals.
    /// Calibrated to RMS pull 1 over the 6588 far-field probes of the 31 non-CdTe detectors;
    /// achieved coverage 83.2%.  The previous 0.016 over-covered by 2.6x.
    constexpr double tot_far_floor = 0.006;
    /// FAR-FIELD TOTAL FLOOR for CdTe/CZT-class crystals - MEASURED, same corpus, 5 detectors,
    /// 1080 probes; 1.971% for RMS pull 1, coverage 81.9%.  Rounded to 0.020.
    ///
    /// This split is EMPIRICAL and its root cause is NOT understood.  Keeping one constant would
    /// mean either under-covering CZT by 3.2x or over-covering everything else by 1.5x, so the
    /// split earns its place - but it is a patch over a modelling gap, not a description of one.
    ///
    /// What is known: the total-efficiency error does NOT track detector size (the corpus's
    /// SMALLEST crystal, a bare 0.5 cm3 CZT, is its most accurate at 0.044%, while a 1.0 cm3 CZT
    /// with one attenuator layer is its worst at 3.490%), and it does NOT track `TotEffTier`
    /// either - `EtaTotTable` spans the entire range from unresolved to 3.490%.  Both hypotheses
    /// were tested and rejected on this corpus.  The leading remaining idea is that an uncollided
    /// kernel plus an angle-flat correction cannot carry scatter-in from the surrounding housing.
    /// See scratch/20260915_nonbare_det_uncert_investigate_prompt.md.
    ///
    /// DELETE THIS CONSTANT once the underlying model is fixed; it should not outlive the defect.
    constexpr double tot_far_floor_cdte = 0.020;
    /// NEAR-FIELD TOTAL FLOOR - NOT re-derived; the 2026-09 study measured the far stratum only.
    constexpr double tot_near_floor = 0.029;
    /// NOT re-derived.  Whether a step at 4 a exists at all is untested: it needs the near
    /// stratum plus a breakpoint scan, neither of which the 2026-09 (far-field) study ran.
    constexpr double near_regime_a = 4.0;
    /// NOT re-derived, and DEAD IN PRACTICE: it is applied only by the opt-in closed loop
    /// (`GenerationOptions::closed_loop`), which nothing in InterSpec ever enables.
    constexpr double generator_floor_inflation = 1.25;
    /// TRANSFER FAR-FIELD ON-AXIS FLOOR - MEASURED, 2026-09 corpus, and it survives.
    /// A transfer anchored on a MEASURED MC curve at one distance and queried on axis at
    /// another runs 0.55% RMS, against this 0.5%.  (The older
    /// curve_transfer_envelope_corpus reported 0.00% here, but it anchored on a golden and
    /// re-queried the same point - that is a tautology, not an accuracy, as its own gate
    /// comment says.  Anchor at one distance and predict at another, or measure nothing.)
    constexpr double transfer_far_onaxis = 0.005;
    /// OFF-AXIS AMPLITUDE - MEASURED, 2026-09 corpus (36 detectors, 7668 held-out points).
    /// The residual SATURATES with angle; see transfer_offaxis_s2_half.  This is the
    /// amplitude the saturating form approaches, not a per-sin^2 slope, so it is not
    /// comparable to the 0.03 it replaces.
    constexpr double transfer_offaxis_mid = 0.037;
    /// EXTRA AMPLITUDE at low energy - MEASURED, same corpus, and far smaller than the
    /// 0.25 it replaces.  Splitting the corpus at mid_e_ref_keV gives A = 3.68% above and
    /// 5.23% below, i.e. the amplitude rises by a factor 1.42, not the factor ~9 the old
    /// value implied.  Enters as (mid + low_e * w^2) with w the ln ramp below.
    constexpr double transfer_offaxis_low_e = 0.016;
    constexpr double transfer_low_e_ref_keV = 45.0;
    constexpr double transfer_mid_e_ref_keV = 150.0;
    /// Half-saturation in sin^2(theta) for the off-axis term - MEASURED, 2026-09 corpus.
    /// 0.153 = sin^2(23 degrees); fitted 0.116 (= sin^2 20 deg) above mid_e_ref_keV and
    /// 0.203 (= sin^2 27 deg) below, so one value serves.
    ///
    /// WHY THE FORM CHANGED.  The error climbs to ~2.4% by 30 degrees and then flattens
    /// (measured 2.43 / 3.64 / 4.01 / 3.04% at 30/45/60/75), while sin^2(theta) grows
    /// without bound - so the old form was ~4x short at 15-30 degrees, where ordinary
    /// measurements are made, and over-covered by ~2x past 60 degrees, where they are not.
    /// The saturating form beat it on HELD-OUT angles in every split tried (fit {15,45,75}
    /// predict {30,60}: RMS 0.49% vs 1.15%; and the reverse: 0.75% vs 1.53%), which is the
    /// bar a change of functional form has to clear.
    ///
    /// <= 0 selects the LEGACY sin^2(theta) form, which is what a stored response written
    /// before this field existed gets on read - old files keep old behaviour exactly.
    ///
    /// HOW MUCH OF THIS IS REDUCIBLE, measured on the same corpus - because the obvious
    /// idea (add angular resolution) buys much less than it looks like it should:
    ///
    ///   off-axis RMS error, no correction                      3.08%
    ///   + a shared saturating shape f(theta)                   2.65%   (free)
    ///   + amplitude linear in crystal aspect ratio L/2R        2.50%   (free)
    ///   + a PERFECT per-(detector, angle) correction           2.28%   <- MC anchor angles
    ///   + a PERFECT per-(detector, angle, ENERGY) correction   0.38%   <- full eta(E,theta)
    ///
    /// Only 45% of the variance is a fixed bias per (detector, angle); the other **55% varies
    /// with ENERGY at fixed detector and angle**.  So angular refinement alone cannot remove
    /// most of it, and a geometry-only correction removes just 19% in RMS - the per-detector
    /// means correlate nicely with aspect ratio (flat crystals near zero, L/2R ~ 1 running
    /// +2 to +4%), but those means are only part of the variance.  What collapses the error is
    /// resolving ENERGY off axis, which is what a full characterization does.
    ///
    /// Predictors TRIED AND REJECTED on this corpus, so they are not retried: crystal
    /// aspect ratio, the kernel's solid-angle-weighted mean chord, and that chord weighted
    /// by interaction probability 1 - exp(-mu(E)L).  All three correlate with the residual
    /// at |r| ~ 0.35-0.38 and all three leave ~1.9% - the interaction weighting buys nothing
    /// over plain geometry even though it does carry real energy dependence (for a 3"x3" NaI
    /// at 45 degrees its ratio runs 0.755 -> 0.871 over 40 keV to 2.5 MeV).  It is the wrong
    /// energy dependence: the entry chord says where a photon first interacts, while
    /// full-energy containment depends on escape FROM that point, which is a different
    /// geometric quantity.  See envelope_chord_predictor in test_CeeLoDrfIntegration.
    constexpr double transfer_offaxis_s2_half = 0.153;
    constexpr double transfer_near_contact = 0.10;
    constexpr double transfer_near_gate_a = 5.0;
}  // namespace model_sigma

// ---------------------------------------------------------------------------
// Geometry descriptor (storable; rebuilds the ray-trace Geometry)
// ---------------------------------------------------------------------------

/// Full material definition so the stored response is self-contained
/// (user-defined materials round-trip; nothing resolves by name at load).
struct MaterialSpec {
    std::string name;
    double density_g_per_cm3 = 0.0;
    std::vector<MaterialComponent> composition;  ///< {Z, mass_fraction}

    static MaterialSpec from(const Material& m);
    Material to_material() const { return Material(name, density_g_per_cm3, composition); }
};

/// One concentric cup-shaped layer (endcap, window, housing, wrap...).
/// Multiple layers are supported, ordered innermost -> outermost; each maps
/// to one Geometry::add_attenuator call. NOTE: CeeLo attenuators cover front
/// + side; dedicated back layers are not modeled (deferred with theta > 90).
struct LayerSpec {
    int material_index = -1;      ///< into GeometryDescriptor::materials
    double front_thickness_cm = 0.0;
    double side_thickness_cm = 0.0;
    double z_start_cm = 0.0;      ///< axial extent (crystal front face = 0)
    double z_end_cm = 0.0;
};

/// Side-only collimator tube (maps to Geometry::add_collimator).
struct CollimatorSpec {
    int material_index = -1;
    double side_thickness_cm = 0.0;
    double z_start_cm = 0.0;      ///< typically negative (extends past the face)
    double z_end_cm = 0.0;
};

enum class ResponseSymmetry : uint8_t {
    Axial,     ///< cylinders: no phi axis
    Quadrant   ///< boxes: phi axis over one quadrant (reflect by symmetry)
};

/// Which plane a caller's `dist_cm` is measured FROM. Distances and the kernel
/// must share one origin -- a mismatch silently corrupts everything (master
/// plan sec 3.2 warning) rather than failing.
///
/// Both data below are planes perpendicular to the axis, and BOTH are measured
/// to the source, so a larger value always means further away:
///
///   CrystalFace  -- the front face of the CRYSTAL SOLID, z = 0 in the crystal
///                   frame. NOT the front of the active volume: with a dead
///                   layer the active volume starts a further t_front behind
///                   this plane (see Geometry::set_dead_layer).
///   EndcapFront  -- the outermost front surface of the attenuator stack, i.e.
///                   the face a user can physically touch. This is InterSpec's
///                   user-facing convention. It sits
///                   endcap_front_offset_cm() in FRONT of CrystalFace.
///
/// The conversion between them is GeometryDescriptor::endcap_front_offset_cm(),
/// applied in DetectorResponse::query_position(). Nothing else should be doing
/// that arithmetic by hand.
enum class ReferencePoint : uint8_t {
    CrystalFace,
    EndcapFront
};

/// Ways a GeometryDescriptor violates a Geometry/RayTrace precondition.
///
/// Those preconditions are asserts, so they vanish in a release build and a
/// violating descriptor traces silent garbage instead of failing. Every
/// descriptor must therefore be run through GeometryDescriptor::problems()
/// before it reaches a Geometry.
enum class GeometryProblem : uint8_t {
    DimensionsMissing,     ///< dimensions_cm too short for the shape
    BulletOnNonCylinder,   ///< fillet requested on a box
    BulletNotFinite,       ///< NaN / inf / negative
    BulletTooWide,         ///< bullet_radius_cm >= crystal radius (rho_c <= 0)
    BulletTooLong,         ///< bullet_radius_cm >= crystal length
    BulletNoDeadLayerRoom, ///< dead layer leaves no active fillet
    DeadLayerTooThick,     ///< dead layer consumes the whole crystal
    BoreOnNonCylinder,
    BoreNotFinite,
    BoreTooWide,           ///< bore radius >= crystal radius
    BoreTooDeep,           ///< bore depth >= crystal length
    BoreTipTooBlunt,       ///< rounded_tip && radius > depth
    BoreOutsideFillet,     ///< bore_fits() fails against the fillet
    BoreInsideDeadLayer    ///< bore radius >= (crystal radius - side dead layer)
};

/// Short English description of `p`; for developer-facing messages (InterSpec
/// maps the enum onto its own localized strings).
const char* to_string(GeometryProblem p);

/// The storable geometry: shape, dimensions, bore, dead layer, layers,
/// collimator, symmetry, reference point + the material table they index.
struct GeometryDescriptor {
    DetectorShape shape = DetectorShape::Cylinder;
    /// The serialized crystal dimensions, in cm: Cylinder {radius, FULL length};
    /// Box {half_x, half_y, FULL length}.  See CRYSTAL DIMENSION CONVENTION in
    /// geometry/Geometry.h -- and prefer set_dimensions() / cylinder_dims() /
    /// box_dims() below, which name the fields, over indexing this by hand.
    ///
    /// This layout is written into every response file; it does not change.
    std::vector<double> dimensions_cm;
    /// Cylinder only: quarter-torus fillet radius on the outer FRONT edge
    /// ("bulletization", ANGLE's `bulletizingRadius`), cm. 0 = a sharp
    /// 90-degree edge -- the default, and bit-for-bit the pre-feature trace.
    /// Must be < radius and < length; see problems().
    double bullet_radius_cm = 0.0;
    int crystal_material_index = -1;
    std::optional<BoreHoleConfig> bore;          ///< coax HPGe finger
    std::optional<DeadLayerConfig> dead_layer;   ///< cm, crystal material
    std::vector<LayerSpec> layers;               ///< innermost -> outermost
    std::optional<CollimatorSpec> collimator;
    ResponseSymmetry symmetry = ResponseSymmetry::Axial;
    ReferencePoint reference_point = ReferencePoint::CrystalFace;
    std::vector<MaterialSpec> materials;

    /// Set `shape` and `dimensions_cm` together, so a Box can never end up
    /// declared with two numbers.  See CRYSTAL DIMENSION CONVENTION.
    void set_dimensions(const CylinderDims& dims);
    void set_dimensions(const BoxDims& dims);

    /// Named read-back of `dimensions_cm`.  Asserts the shape matches and the
    /// vector is long enough; consult problems() first for untrusted input.
    CylinderDims cylinder_dims() const;
    BoxDims box_dims() const;

    /// Build the ray-trace geometry. The returned Geometry references the
    /// Material instances appended to `owned` -- keep them alive as long as
    /// the Geometry (DetectorResponse holds both).
    Geometry build_geometry(std::vector<std::unique_ptr<Material>>& owned) const;

    /// Transverse half-extent `a` -- the near/far regime SCALE PARAMETER
    /// (cylinder: outer radius; box: half-diagonal), not a physical dimension.
    ///
    /// It intentionally includes the side dead layer, even though the dead
    /// layer is internal to the crystal (Geometry::set_dead_layer), so it runs
    /// slightly large. That is harmless and must stay: `a` is computed from the
    /// descriptor at BOTH generation time (ResponseGenerator picks its far
    /// reference distance as far_distance_a * a) and query time (the
    /// near-regime test is d < near_regime_a * a). Changing the value would
    /// desynchronise every already-generated response from the code that made
    /// it. Contrast endcap_front_offset_cm(), which is query-time only and so
    /// could be -- and was -- corrected in place.
    ///
    /// If you need a true physical outer radius, sum the layer side
    /// thicknesses onto dimensions_cm yourself; do not use this.
    double transverse_half_extent() const;

    /// |z_min| of the outermost front surface: the (positive) offset from the
    /// endcap front to the CRYSTAL SOLID's face (z = 0), i.e. the summed front
    /// thicknesses of the attenuator shells. The dead layer is deliberately
    /// excluded -- it is carved out of the inside of the crystal, so it sets
    /// where the ACTIVE volume starts, not where the crystal face is.
    double endcap_front_offset_cm() const;

    /// Crystal K-edges within (e_min, e_max) keV -- the mandatory segment
    /// breaks for any ln-eta(E) interpolation (both flanks become nodes).
    std::vector<double> crystal_k_edges(double e_min_keV,
                                        double e_max_keV) const;

    /// Every Geometry/RayTrace precondition this descriptor violates; empty
    /// means it is safe to build. See GeometryProblem.
    std::vector<GeometryProblem> problems() const;

    /// Standalone XML for JUST the geometry, so a host can store a detector's
    /// shape before (or without) any response having been generated for it.
    ///
    /// The payload is the same <Detector> element DetectorResponse's own codec
    /// writes -- shared code, so the two cannot drift -- wrapped in a
    /// <CeeLoGeometry> root. Round-trips through from_xml_string().
    std::string to_xml_string() const;
    static GeometryDescriptor from_xml_string(const std::string& xml);
};

// ---------------------------------------------------------------------------
// Stored mu tables (generation-time attenuation snapshot)
// ---------------------------------------------------------------------------

/// Per-material macroscopic attenuation table, sampled at generation on a
/// log-energy grid with absorption-edge flank pairs. Evaluation is log-log
/// linear between samples (flank pairs keep interpolation from bridging an
/// edge). Stored so the query-time kernel uses the SAME mu(E) as the MC that
/// produced eta -- eta*K must not drift when cross-section data changes.
struct MuTable {
    int material_index = -1;              ///< into GeometryDescriptor::materials
    std::vector<double> energy_keV;       ///< ascending, with edge flanks
    std::vector<double> mu_pe, mu_cs, mu_rs, mu_pp;  ///< 1/cm (linear values)

    MacroscopicXS eval(double energy_keV) const;

    /// Sample `mat` on a log grid of ~n_per_decade points/decade over
    /// [e_min, e_max] keV plus K/L-edge flank pairs for its elements.
    static MuTable sample(const Material& mat, int material_index,
                          double e_min_keV = 10.0, double e_max_keV = 10000.0,
                          int n_per_decade = 45);
};

// ---------------------------------------------------------------------------
// Residual tables eta(E, cos-theta [, phi])
// ---------------------------------------------------------------------------

/// Tensor table of ln eta = ln(eps / K) with per-node MC sigma.
/// Interpolation (spec sec 3.1): per (ct, phi) node, log-log PCHIP in E
/// SEGMENTED at the crystal K-edges (both flanks are nodes; no interpolant
/// spans an edge); then PCHIP across cos-theta; then linear across phi
/// (boxes only). Outside the node range: clamp (never extrapolate).
class EtaTable {
public:
    std::vector<double> energies_keV;   ///< ascending; edge flanks as nodes
    std::vector<double> cos_thetas;     ///< ascending, e.g. 0.02 .. 1.0
    std::vector<double> phis_deg;       ///< empty (axial) or quadrant nodes
    std::vector<double> ln_eta;         ///< [e][c][p] energy-major flattened
    std::vector<double> frac_sigma;     ///< per-node MC fractional sigma
    std::vector<double> edges_keV;      ///< crystal K-edges (segment breaks)

    bool empty() const { return energies_keV.empty(); }
    size_t index(size_t e, size_t c, size_t p) const {
        const size_t np = phis_deg.empty() ? 1 : phis_deg.size();
        return (e * cos_thetas.size() + c) * np + p;
    }

    /// Build the per-(ct,phi) segmented energy interpolants. Must be called
    /// after filling / deserializing and before eval().
    void finalize();

    /// ln eta at (E, ct, phi); clamps and reports it via `clamped`.
    double eval_ln(double energy_keV, double cos_theta, double phi_deg,
                   bool& clamped) const;

    /// Interpolated fractional MC sigma at the query point (bilinear in
    /// (lnE, ct) at the nearest phi node).
    double node_frac_sigma(double energy_keV, double cos_theta,
                           double phi_deg) const;

private:
    struct SegCurve {                    // one (ct, phi) node's energy curve
        std::vector<Pchip> segs;         // per K-edge segment, over (lnE, ln eta)
        std::vector<double> seg_lo, seg_hi;  // lnE range per segment
        double eval(double lnE, bool& clamped) const;
    };
    std::vector<SegCurve> curves_;       // [c][p] flattened c-major
};

// ---------------------------------------------------------------------------
// Near-field model (spec Eq. 6)
// ---------------------------------------------------------------------------

/// Tabulated near-field multiplier: ln N on a per-shape-energy tensor grid
/// (cos_theta x distance-from-crystal-face, cm), with per-node MC sigma.
/// Interpolation: PCHIP over ln(d) per (E, cos_theta) node (built in
/// finalize()), PCHIP across cos_theta, linear in ln E between shape
/// energies. Clamped outside the grid (never extrapolated); the outermost
/// distance node is an ln N = 0 anchor so N fades smoothly to 1 at the
/// breakpoint gate. N = 1 outside the measured breakpoint d_break(E, theta).
/// Table is measured at phi = 0 (axial); boxes reuse it for all phi
/// (documented limitation -- the only box fixture is far-field profile).
struct NearFieldModel {
    std::vector<double> energies_keV;   ///< shape energies, ascending
    std::vector<double> cos_thetas;     ///< ascending, e.g. 0.02 .. 1.0
    std::vector<double> dists_cm;       ///< ascending, from crystal-face origin
    std::vector<double> ln_n;           ///< [e][c][d] energy-major flattened
    std::vector<double> frac_sigma;     ///< per-node MC fractional sigma
    /// Breakpoint distances (cm, from crystal-face origin) on
    /// (energies_keV x break_cos_thetas), row-major [e][c]:
    std::vector<double> break_cos_thetas;
    std::vector<double> break_d_cm;

    bool empty() const { return energies_keV.empty(); }
    size_t index(size_t e, size_t c, size_t d) const {
        return (e * cos_thetas.size() + c) * dists_cm.size() + d;
    }

    /// Build the per-(E, cos_theta) ln(d) interpolants. Must be called after
    /// filling / deserializing and before ln_boost() / node_frac_sigma().
    void finalize();

    /// ln N at (E, cos_theta, d); clamps to the grid on every axis.
    double ln_boost(double energy_keV, double cos_theta, double d_cm) const;

    /// Interpolated fractional MC sigma at the query point (trilinear on
    /// sigma^2 in (lnE, cos_theta, ln d)).
    double node_frac_sigma(double energy_keV, double cos_theta,
                           double d_cm) const;

    double breakpoint_d_cm(double energy_keV, double cos_theta) const;

private:
    std::vector<Pchip> d_curves_;   ///< [e][c] flattened e-major, over (ln d, ln N)
};

// ---------------------------------------------------------------------------
// eps_tot tier (spec Eq. 4)
// ---------------------------------------------------------------------------

enum class TotEffTier : uint8_t {
    KernelExact,   ///< bare crystal: eps_tot = K_{mu - mu_RS}
    BCurve,        ///< canned scintillator: eps_tot = b(E) * K_{mu - mu_RS}
    EtaTotTable    ///< HPGe-class: eps_tot = k(E) * eta_tot(E,theta) * K
};

struct TotEffPayload {
    TotEffTier tier = TotEffTier::KernelExact;
    std::vector<double> b_energies_keV;   ///< BCurve: ~8 log-spaced nodes
    std::vector<double> ln_b;
    EtaTable eta_tot;                     ///< EtaTotTable tier only

    void finalize();
    double ln_b_at(double energy_keV) const;   ///< clamped PCHIP over (lnE, ln b)

private:
    Pchip b_curve_;
};

// ---------------------------------------------------------------------------
// Grounding (spec Eqs. 7a-8b)
// ---------------------------------------------------------------------------

/// One measured calibration point, kept verbatim for provenance and
/// re-grounding. `model_eff` is the simulated efficiency at the SAME
/// geometry, so k_i = measured / model (Eq. 7a).
struct GroundingPoint {
    double energy_keV = 0.0;
    double measured_eff = 0.0;        ///< absolute FEP efficiency
    double model_eff = 0.0;           ///< eps_sim at the point's geometry
    double frac_stat_sigma = 0.0;     ///< independent (peak area, BR)
    double frac_cert_sigma = 0.0;     ///< 100% correlated within source_key
    std::string source_key;           ///< calibration-source identity
    double distance_cm = 0.0;         ///< point's own geometry
    double cos_theta = 1.0;
    double phi_deg = 0.0;
};

/// Empirical grounding-transfer inflation sigma_transfer(d, theta, E):
/// ~0 on-axis far-field, growing off-axis / close-in / toward low E
/// (S7-measured; constants from the spec sec 4 grounding table -- Level-1
/// values; a Level-2 nuisance fit would shrink the near term to ~1%).
struct SigmaTransferModel {
    double far_onaxis = model_sigma::transfer_far_onaxis;        ///< far-field on-axis floor
    double offaxis_mid = model_sigma::transfer_offaxis_mid;      ///< x sin^2(theta), mid/high E
    double offaxis_low_e = model_sigma::transfer_offaxis_low_e;  ///< extra x sin^2(theta) at low E
    double low_e_ref_keV = model_sigma::transfer_low_e_ref_keV;  ///< where the low-E term is fully on
    double mid_e_ref_keV = model_sigma::transfer_mid_e_ref_keV;  ///< where the low-E term is off
    double near_contact = model_sigma::transfer_near_contact;    ///< at contact (d ~ a), no Level-2
    double near_gate_a = model_sigma::transfer_near_gate_a;      ///< near term active below this many a
    /// Half-saturation in sin^2(theta) of the off-axis term.  <= 0 selects the legacy
    /// unbounded sin^2(theta) form, so a response deserialized from a file written before
    /// this field existed behaves exactly as it did.
    double offaxis_s2_half = model_sigma::transfer_offaxis_s2_half;

    /// The three mechanisms separately - the on-axis floor, the off-axis (angle-flat eta)
    /// residual and the near-field residual; eval() is their quadrature sum.  A covariance treats
    /// each as its own fully-correlated common mode (rank-one block): their per-energy magnitudes
    /// differ (the off-axis term ramps up below mid_e_ref_keV), and one combined block would
    /// over-correlate energies whose magnitudes differ.
    struct Components { double far_onaxis = 0.0, offaxis = 0.0, near = 0.0; };
    Components components(double d_over_a, double cos_theta, double energy_keV) const;

    /// d in units of the transverse half-extent a.
    double eval(double d_over_a, double cos_theta, double energy_keV) const;
};

/// k(E) grounding fit: ln k on hat-basis knots in ln E, with the full GLS
/// covariance (Eq. 7c-7d). Linear in the coefficients, so
/// Cov[ln k(E), ln k(E')] = B(E) C B(E')^T (Eq. 8a) is exact.
struct GroundingBlock {
    bool curve_derived = false;     ///< sampled from a fitted legacy curve
                                    ///< (lower quality) instead of raw peaks
    std::vector<GroundingPoint> points;
    std::vector<double> knot_ln_energies;  ///< hat-basis knots (ln keV)
    std::vector<double> ln_k;              ///< coefficients at knots
    std::vector<double> cov;               ///< row-major NxN covariance
    SigmaTransferModel transfer;

    bool empty() const { return knot_ln_energies.empty(); }
    /// ln k at E (clamped to the knot range; sets `clamped` outside it).
    double eval_ln_k(double energy_keV, bool& clamped) const;
    /// Var[ln k(E)] from the fit covariance.
    double var_ln_k(double energy_keV) const;
    /// Cov[ln k(E1), ln k(E2)].
    double cov_ln_k(double e1_keV, double e2_keV) const;
};

// ---------------------------------------------------------------------------
// Uncertainty floors + provenance
// ---------------------------------------------------------------------------

/// Coverage-tuned per-{quantity x regime} model floors (fractional 1-sigma;
/// spec sec 4/5). Defaults are the campaign's conservative envelope
/// (model_sigma); the generator inflates the FEP pair on a minor model-form
/// failure of its closed loop.
struct SigmaFloors {
    double fep_far = model_sigma::fep_far_floor;
    double fep_near = model_sigma::fep_near_floor;
    double tot_far = model_sigma::tot_far_floor;
    double tot_near = model_sigma::tot_near_floor;
    double near_regime_a = model_sigma::near_regime_a;   ///< near regime: d < this many a
};

enum class ResponseProfile : uint8_t {
    FarField,   ///< eta(E,theta) only; no near model; near queries flag
    General,    ///< + near-field (Omega, theta_eff) model  [default]
    Contact     ///< + denser near scan / near eps_tot table
};

const char* to_string(ResponseProfile p);

/// How a response was produced. Recorded rather than inferred: a quick-MC
/// transfer and a measured-curve transfer both leave an angle-flat eta table
/// and a model_transfer envelope, so the payload alone cannot tell them apart,
/// and a host that wants to show a stored detector's current method (rather
/// than a default) would be guessing.
enum class ProductionMethod : int {
    FullMc = 0,           ///< ResponseGenerator::generate(), full energy x angle(x distance) scan
    QuickMcTransfer = 1,  ///< ResponseGenerator::generate() with transfer_mode
    CurveTransfer = 2     ///< make_transfer_response(): no Monte Carlo at all
};

const char* to_string(ProductionMethod m);

struct ResponseProvenance {
    ProductionMethod method = ProductionMethod::FullMc;
    std::string ceelo_version;        ///< library version string
    std::string created_utc;          ///< ISO-8601, informational
    ResponseProfile profile = ResponseProfile::General;
    double node_fep_precision = 0.003;///< per-node MC precision target
    /// Half-width (keV) of the full-energy-peak window this response's FEP was
    /// scored with (physics/FepWindow.h).  A model that credits in-window
    /// Compton must use the SAME window, or it is calibrated against the wrong
    /// truth - which is why this travels with the response instead of being
    /// assumed.
    double fep_window_keV = kDefaultFepWindowKeV;
    uint64_t generation_seed = 0;     ///< base seed (per-node seeds derive)
    int kernel_n_rays = 2048;         ///< quadrature rays used in evaluation
    double valid_e_min_keV = 0.0, valid_e_max_keV = 0.0;
    double min_distance_cm = 0.0;     ///< d-validity floor (far-field profile)
    std::string detector_name;        ///< informational label
};

// ---------------------------------------------------------------------------
// Accuracy certificate (persisted metadata)
// ---------------------------------------------------------------------------

/// An honest, persisted record of how well the assembled response reproduces a
/// fresh MC probe bank -- the output of ResponseGenerator::certify(). It is
/// METADATA ABOUT the response, not part of the response content: it is
/// serialized in the XML (so a stored DRF carries its own scorecard) but is
/// EXCLUDED from content_hash(), so the same response with or without a
/// certificate hashes identically. Absent by default (empty()).
struct AccuracyCertificate {
    bool converged = false;       ///< D-b: closed loop met tolerance;
                                  ///< D-a: true when the probe pass ran
    int iterations = 0;           ///< refinement iterations (D-a: 0)
    double cpu_seconds = 0.0;     ///< certificate probe-bank MC cost
    uint64_t probe_seed_base = 0; ///< generation base seed the probes derive from
    /// |model/mc - 1| percentiles over converged probes (mc_sig/mc <= 0.05):
    double fep_median = 0.0, fep_p95 = 0.0, fep_max = 0.0;
    double tot_median = 0.0, tot_p95 = 0.0;

    /// One scored probe. `tag` classifies the probe family (0 = random; D-b
    /// adds structured tags). `pass` is the noise-aware tolerance verdict.
    struct Row {
        double E_keV = 0.0, d_cm = 0.0, cos_theta = 1.0, phi_deg = 0.0;
        /// Full-energy peak: the MC truth and the model, each with its sigma.
        double mc = 0.0, mc_sig = 0.0, model = 0.0, model_sig = 0.0;
        /// Total efficiency, the same four.  Kept per row rather than summarized
        /// away because the tot_* regime floors have to be DERIVED from this
        /// distribution, and two percentiles cannot be deconvolved against the
        /// MC noise that produced them.
        double mc_tot = 0.0, mc_tot_sig = 0.0, model_tot = 0.0, model_tot_sig = 0.0;
        uint8_t tag = 0;
        bool pass = false;
    };
    std::vector<Row> rows;

    bool empty() const { return rows.empty(); }
};

// ---------------------------------------------------------------------------
// The response object
// ---------------------------------------------------------------------------

/// Storable + evaluable detector response. Non-copyable (owns the Material
/// instances its Geometry references); share via shared_ptr.
///
/// Thread-safety: all eval methods are const and touch only immutable state
/// after construction/finalize(); safe for concurrent queries.
class DetectorResponse {
public:
    DetectorResponse() = default;
    DetectorResponse(const DetectorResponse&) = delete;
    DetectorResponse& operator=(const DetectorResponse&) = delete;

    // --- assembly (generator / deserialization API) ---
    GeometryDescriptor descriptor;
    std::vector<MuTable> mu_tables;      ///< one per referenced material
    EtaTable eta_fep;
    NearFieldModel near_field;
    TotEffPayload tot_eff;
    GroundingBlock grounding;
    SigmaFloors floors;
    ResponseProvenance provenance;

    /// Optional geometry-transfer sigma envelope, applied UNCONDITIONALLY when
    /// present (unlike GroundingBlock::transfer, which needs grounding). Set by
    /// the EFFTRAN transfer producers (make_transfer_response / transfer_mode)
    /// so an ungrounded angle-flat response inflates sigma off-axis/near, where
    /// the un-modeled eta(E,theta) residual lives. Absent by default -> stored
    /// responses are byte-identical.
    std::optional<SigmaTransferModel> model_transfer;

    /// Dead-layer / endcap Compton scatter-in recapture fraction for the TOTAL
    /// kernel (see ResponseKernel.h). Folds forward-Compton scatter-in from the
    /// passive layers into K_noRayleigh, correcting the near-field total
    /// under-prediction for thick-dead-layer (HPGe-class) detectors. Applied to
    /// both the anchor and target kernels so the transfer ratio stays consistent.
    /// 0 (default) is bit-identical to the removal-only kernel; the transfer
    /// producers set the calibrated value. FEP is unaffected.
    double scatter_in_recapture = 0.0;

    /// Persisted accuracy scorecard (ResponseGenerator::certify). Additive
    /// metadata: serialized in the XML but EXCLUDED from content_hash, so a
    /// response is identical with or without it. Empty by default.
    AccuracyCertificate certificate;

    /// Shielded-eps_tot build-up model (Stage E3 A1 seam).  NON-serialized,
    /// settable, default null.  When set AND a caller passes a ShieldContext to
    /// an eps_total_* overload, eps_total is multiplied by max(1, buildup_model)
    /// and the build-up sigma floor is folded in.  Null / no-ShieldContext
    /// callers get byte-identical behavior.  CeeLo leaves it null; InterSpec
    /// installs a GadrasShieldScatter-backed ratio model.
    BuildupModel buildup_model;

    /// Build the Geometry + lookup structures from the descriptor/tables.
    /// Call after filling the public fields (from_xml_string does it).
    void finalize();
    bool finalized() const { return !owned_materials_.empty() || geometry_built_; }

    // --- geometry access ---
    const Geometry& geometry() const { return geometry_; }
    double transverse_half_extent() const { return descriptor.transverse_half_extent(); }

    /// Convert a (theta, phi, distance-from-reference-point) query to a source
    /// position in the crystal frame (z = 0 at the CRYSTAL SOLID's face, source
    /// in front at negative z).
    ///
    /// `dist_cm` is interpreted in descriptor.reference_point's convention;
    /// this is the ONLY place that applies endcap_front_offset_cm(). Every
    /// caller that builds a source position from a user-supplied distance must
    /// go through here (or reproduce it exactly, as
    /// CeeLoUtils::makeTransferResponse does for the anchor position) --
    /// otherwise the response is queried at a different place than it was
    /// generated for, which biases results silently.
    Eigen::Vector3d query_position(double theta_rad, double phi_rad,
                                   double dist_cm) const;

    /// The REFERENCE POINT itself, in the crystal frame: the origin for a
    /// CrystalFace response, (0, 0, -endcap_front_offset_cm()) for EndcapFront.
    /// Equivalently query_position(0, 0, 0).
    ///
    /// A caller that needs a DIRECTION between a queried source position and
    /// "the detector" must use this, not the origin: the two differ by the
    /// endcap offset, which is a several-degree parallax for a source at
    /// contact. Aiming at the origin instead tilts anything built on that
    /// direction (e.g. an aperture fan mapped into the caller's own frame)
    /// away from the geometry the caller thinks it is describing.
    Eigen::Vector3d reference_point_position() const {
        return query_position(0.0, 0.0, 0.0);
    }

    // --- point-source evaluation ---
    /// dist_cm measured from descriptor.reference_point along the (theta,
    /// phi) direction; theta from the detector axis (0 = on-axis front).
    EffResult eps_fep(double energy_keV, double theta_rad, double phi_rad,
                      double dist_cm) const;
    EffResult eps_total(double energy_keV, double theta_rad, double phi_rad,
                        double dist_cm) const;

    /// Position variants (crystal-face frame). The quadrature-taking
    /// overloads let fit loops with fixed positions reuse the (expensive)
    /// quadrature across energies.
    ApertureQuadrature make_quadrature(const Eigen::Vector3d& src_cm) const;
    EffResult eps_fep_at(double energy_keV, const Eigen::Vector3d& src_cm) const;
    /// `sc` (optional, default null) applies the build-up seam: when non-null
    /// and buildup_model is set, eps_total is scaled by max(1, buildup_model)
    /// and its sigma inflated.  Passing null is byte-identical to no seam.
    EffResult eps_total_at(double energy_keV, const Eigen::Vector3d& src_cm,
                           const ShieldContext* sc = nullptr) const;
    EffResult eps_fep_at(double energy_keV, const Eigen::Vector3d& src_cm,
                         const ApertureQuadrature& q) const;
    EffResult eps_total_at(double energy_keV, const Eigen::Vector3d& src_cm,
                           const ApertureQuadrature& q,
                           const ShieldContext* sc = nullptr) const;

    /// Per-element extended-source evaluation (spec Eq. 5): eps_fep at an
    /// element position with a per-ray source-transmission factor folded
    /// into the kernel. `t_src(dir)` must return the source-geometry
    /// survival factor along `dir` FROM the element -- computed with the
    /// FEP survival removal mu (fep_survival_removal_mu; prefer the
    /// material-aware kn_in_window_fraction(E, win, mat) f_win when the
    /// source/shield material is known) for eps_fep, or mu_total for
    /// eps_total.
    EffResult eps_fep_element(double energy_keV, const Eigen::Vector3d& src_cm,
                              const ApertureQuadrature& q,
                              const std::function<double(const Eigen::Vector3d&)>& t_src) const;
    EffResult eps_total_element(double energy_keV, const Eigen::Vector3d& src_cm,
                                const ApertureQuadrature& q,
                                const std::function<double(const Eigen::Vector3d&)>& t_src,
                                const ShieldContext* sc = nullptr) const;

    // --- kernel with the STORED mu tables ---
    /// K(E) over a quadrature using the generation-time MuTables (NOT the
    /// live cross-section data): the FEP kernel with MuChoice::Total, the
    /// eps_tot kernel with MuChoice::NoRayleigh.
    double kernel_K(double energy_keV, const ApertureQuadrature& q, MuChoice mu,
                    const std::function<double(const Eigen::Vector3d&)>* t_src = nullptr) const;

    /// The kernel, decomposed so a host can supply a per-ray source transmission of its OWN type -
    /// e.g. an autodiff scalar, which the std::function above cannot carry.
    ///
    ///     K = sum_i w_out[i] * t_src(dirs_out[i])
    ///
    /// and with no t_src, sum(w_out) == kernel_K(...).  The weights depend only on (energy, ray,
    /// stored mu tables): no fit parameter enters them, so a host may hold them fixed while
    /// differentiating through its own t_src.
    ///
    /// Use the fep_/total_ pair rather than picking a MuChoice: the eps_total kernel is
    /// TIER-DEPENDENT (the EtaTotTable tier uses MuChoice::Total, not NoRayleigh) and applies
    /// scatter_in_recapture, and getting that wrong is silent.
    void fep_ray_weights(double energy_keV, const ApertureQuadrature& q,
                         std::vector<double>& w_out,
                         std::vector<Eigen::Vector3d>& dirs_out) const;
    void total_ray_weights(double energy_keV, const ApertureQuadrature& q,
                           std::vector<double>& w_out,
                           std::vector<Eigen::Vector3d>& dirs_out) const;

    /// Per-ray FEP interaction probability from the STORED mu tables, PARALLEL to `q.rays` (0 for
    /// a ray without an active chord) - the counterpart of fep_ray_weights for an etendue line set
    /// (io/DetectorEtendue.h), where the host keeps per-line bookkeeping and needs one entry per
    /// line rather than the compacted weights.  Identity, pinned by tests/test_detector_etendue:
    ///     sum_i q.rays[i].omega_w * p_out[i] == sum(w_out) of fep_ray_weights(...)
    void fep_line_probabilities(double energy_keV, const ApertureQuadrature& q,
                                std::vector<double>& p_out) const;

    /// Everything in eps_fep EXCEPT the kernel: exp(ln_eta + ln_N + ln_k), with the flag and the
    /// fractional sigma the full query would have reported.  So
    ///     eps_fep(E, pos) == fep_prefactor(E, pos, q).value * kernel_K(E, q, MuChoice::Total)
    /// and a host assembling K itself keeps the near-field/off-axis/grounding physics - and the
    /// validity flag - rather than silently dropping them.
    EffResult fep_prefactor(double energy_keV, const Eigen::Vector3d& src_cm,
                            const ApertureQuadrature& q) const;

    /// The eps_total counterpart of #fep_prefactor (tier-dependent; 1.0 for the bare-crystal tier).
    EffResult total_prefactor(double energy_keV, const Eigen::Vector3d& src_cm,
                              const ApertureQuadrature& q) const;

    /// Passive-layer transmission envelope over a quadrature with the stored
    /// mu tables (used for the collimator hole-fraction gate).
    double kernel_transmitted(double energy_keV, const ApertureQuadrature& q) const;

    // --- multi-energy covariance for fits (spec Eq. 8b) ---
    /// Row-major NxN fractional covariance of eps_fep between `energies_keV` at ONE query
    /// geometry (a crystal-face-frame position, as eps_fep_at takes):
    ///
    ///     C_ij = sum_t m_t[i] * m_t[j] + Cov[ln k](E_i, E_j) + delta_ij * node2_i
    ///
    /// m_t are the model-envelope terms of the per-query sigma budget (regime floor,
    /// behind-plane, collimator shadow, near-field-unmodeled, and the far / off-axis / near
    /// components of model_transfer and grounding.transfer): each says "the model may be off by
    /// x% here" for a reason shared by every energy, so each is a fully-correlated common mode
    /// that must NOT average down over a fit's peaks.  node2 is the MC node variance (independent
    /// per energy) and Cov[ln k] the anchor / grounding fit covariance - the data-derived part.
    ///
    /// INVARIANT (pinned by tests): C_ii == (eps_fep_at(E_i, src_cm, q).sigma / value)^2 - the
    /// same budget through the same code (fep_budget), so the two cannot drift.  PSD by
    /// construction.  When `model_part` is given it receives the envelope-only matrix
    /// sum_t m_t[i] m_t[j] (same layout), so a host can separate the ad hoc envelopes from the
    /// uncertainty its own data supports.
    ///
    /// The quadrature only feeds the collimator shadow gate; the overload without one builds it
    /// (as eps_fep_at does).  The distance form mirrors eps_fep(E, theta, phi, dist) and goes
    /// through query_position().  eps_total has no covariance API: hosts consume its value only;
    /// its budget is the same struct, so one could be added the same way.
    std::vector<double> frac_covariance(const std::vector<double>& energies_keV,
                                        const Eigen::Vector3d& src_cm,
                                        const ApertureQuadrature& q,
                                        std::vector<double>* model_part = nullptr) const;
    std::vector<double> frac_covariance(const std::vector<double>& energies_keV,
                                        const Eigen::Vector3d& src_cm,
                                        std::vector<double>* model_part = nullptr) const;
    std::vector<double> frac_covariance(const std::vector<double>& energies_keV,
                                        double theta_rad, double phi_rad, double dist_cm,
                                        std::vector<double>* model_part = nullptr) const;

    // --- XML (one codec for generator + InterSpec) ---
    /// Root element <CeeLoResponse version="1">; InterSpec convention
    /// "version" attribute. Doubles printed with max_digits10 so a
    /// save/load/save cycle is string-identical.
    static const int sm_xmlSerializationVersion;  // = 1
    std::string to_xml_string() const;
    static std::shared_ptr<DetectorResponse> from_xml_string(const std::string& xml);

    /// FNV-1a hash of the canonical XML payload (stable content identity;
    /// grounding and grid changes change it). INVARIANT to the accuracy
    /// certificate (which is metadata about the response, not content), so a
    /// response hashes the same with or without one.
    uint64_t content_hash() const;

private:
    /// Serialize to XML; `include_certificate` gates the <Certificate> element
    /// so content_hash() can hash the certificate-free payload (invariance).
    std::string serialize_xml(bool include_certificate) const;

    struct SigmaBudget;  // one query's fractional sigma budget: data-derived vs model envelopes
    struct EvalCommon;   // internal per-query bundle
    /// Geometry, flags and the geometry-only envelope terms (behind-plane, model_transfer,
    /// collimator shadow) shared by the FEP and total paths.
    EvalCommon common_eval(double energy_keV, const Eigen::Vector3d& src_cm,
                           const ApertureQuadrature& q) const;
    /// THE one FEP sigma budget: the near-field gate, grounding (k and its transfer envelope),
    /// the eta node sigma and the regime floor, on top of what common_eval filled in.  Raises
    /// flags on `ec`; returns the near-field boost and grounding ln k for the value.
    /// fep_prefactor (hence eps_fep) and frac_covariance are its only callers - which is what
    /// makes the covariance diagonal equal the per-query sigma by construction.
    void fep_budget(double energy_keV, EvalCommon& ec, double& ln_N, double& ln_k) const;
    /// Shared ray loop behind kernel_K and the *_ray_weights accessors, so the decomposition and
    /// the thing it decomposes cannot drift apart.
    void kernel_ray_weights_impl(double energy_keV, const ApertureQuadrature& q, MuChoice mu,
                                 double recap, std::vector<double>& w_out,
                                 std::vector<Eigen::Vector3d>& dirs_out) const;

    EffResult eps_fep_impl(double energy_keV, const Eigen::Vector3d& src_cm,
                           const ApertureQuadrature& q,
                           const std::function<double(const Eigen::Vector3d&)>* t_src) const;
    EffResult eps_total_impl(double energy_keV, const Eigen::Vector3d& src_cm,
                             const ApertureQuadrature& q,
                             const std::function<double(const Eigen::Vector3d&)>* t_src,
                             const ShieldContext* sc = nullptr) const;

    std::vector<std::unique_ptr<Material>> owned_materials_;
    Geometry geometry_;
    bool geometry_built_ = false;
    // Material* (as used in quadrature rays) -> mu-table index.
    std::vector<std::pair<const Material*, size_t>> mat_to_mu_;
};

} // namespace ceelo

#endif // CEELO_IO_DETECTOR_RESPONSE_H
