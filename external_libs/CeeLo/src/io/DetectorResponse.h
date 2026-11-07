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
/// k(E)*eta_tot(E,theta)*K -- the tier is chosen at generation time.  Close in,
/// any tier is multiplied by a measured near-field N_tot when the response has
/// one (TotEffPayload::near_field; General/Contact profiles), as eps_fep is by N.
///
/// Every query returns {value, sigma, flag}; sigma combines interpolated
/// node-MC variance, coverage-tuned per-regime model floors, grounding
/// covariance and the empirical grounding-transfer inflation (spec sec 5).
///
/// XML: to_xml_string()/from_xml_string() -- ONE codec shared by the
/// generator and by InterSpec (which embeds the fragment in its DRF XML).

#include "geometry/ConvexSeparation.h"
#include "geometry/Geometry.h"
#include "geometry/GeometryDescriptor.h"
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
/// THE CORPUS.  Values marked MEASURED come from scoring generated responses against fresh,
/// never-fitted MC on a stratified (E, d, theta) grid: 36 detectors - 4 CeeLo presets, 10 ANGLE
/// imports, 22 GADRAS Detector.dat imports - spanning NaI / HPGe / LaBr3 / CZT and cylinder / box,
/// generation and probes both at node_fep_precision = 0.003.  Harness: the `--envelope-*` flags on
/// InterSpec's test_CeeLoDrfIntegration (`envelope_corpus_measure`, `envelope_transfer_from_mc`);
/// raw per-probe CSVs and analysis under scratch/20260913_envelope_study/.  Re-run it when the
/// engine changes.
///
/// THE STATISTIC.  These are consumed by a chi-square (InterSpec's
/// compute_efficiency_whitening), so each measured value is the floor that makes the RMS pull
/// about MC truth equal 1 over the relevant probes, with the response's own declared node sigma in
/// the denominator.  That is NOT a 68% interval - the residual is heavy-tailed, so at these values
/// about 75% (FEP) and 84% (total) of probes fall inside one sigma, and the two calibrations
/// differ by ~1.6x.  State which one you targeted if you re-derive.
///
/// NOT INCLUDED, deliberately: CeeLo-vs-GEANT4 model error, bounded at <= 0.48% (FEP) and
/// <= 0.23% (total) on the bare/lightly-shielded G4 configs {1,2,3,5,6,25,26} from the committed
/// references in tests/data/{geant4,ceelo}_reference (scratch/20260913_envelope_study/
/// mc_vs_geant4.py).  It is a bound rather than a measurement - the two codes' own counting noise
/// is ~73% of the observed spread - it is code-vs-code rather than against data, and its
/// applicability to an arbitrary user detector is not established.  A careful comparison against
/// measured data is what would justify folding it in.
namespace model_sigma {
    /// theta > 90 degrees: no nodes behind the face plane; the value is a clamped guess.
    ///
    /// MEASURED, and deliberately NOT reduced.  Direct MC behind the face plane (three
    /// detectors, 100-180 degrees, 60 keV to 1.3 MeV - `envelope_refuse_grade_terms`, since
    /// `probe_bank` samples cos_theta over [cos_theta_min, 1] and never goes there) puts the
    /// clamped guess within 1-4% of truth, against this 30%.
    ///
    /// That measures the model against the solid as DESCRIBED, and the descriptor is knowingly
    /// incomplete exactly there: CeeLo models no attenuator behind the crystal (see LayerSpec),
    /// which is where a real detector keeps its PMT, cryostat or electronics.  So the 1-4%
    /// says the extrapolation is sound for a bare back, not that a query behind a real
    /// detector is good to 4%.  The envelope stays wide because what it is covering is the
    /// missing back structure, which this measurement cannot see.
    ///
    /// The query also raises ResponseFlag::NeedsMc, so a host can decline rather than use it.
    ///
    /// Since back structure can be described (a layer's z_end past the envelope inside it - see
    /// LayerSpec), this applies only where that measurement does not reach: a source behind the
    /// crystal's BACK plane on a detector whose descriptor describes nothing there.  Beside the
    /// detector, and behind one with a described back, the terms below apply instead.
    constexpr double behind_plane = 0.30;

    /// Behind the face plane (cos_theta < 0), on a response whose angular tables stop at the
    /// grazing node and are not angle-flat: the eta and near-field corrections are the grazing
    /// node's, extrapolated.  ResponseFlag::OutOfRangeClamped.
    ///
    /// Measured (2026-11, InterSpec's test_BackIncidence BackEnvelopeMeasure) against CeeLo MC
    /// beside and behind two HPGe and two NaI detectors (2332 points, 40-1332 keV, contact to far
    /// field; NeedsMc points excluded): 0.05 leaves every region's RMS pull at or below 0.77, the
    /// same coverage 0.08 gave.
    constexpr double behind_extrapolation = 0.05;

    /// Behind the crystal's back plane on a detector whose back structure IS described: the
    /// description is approximate (a couple of layers where a real detector has a cryostat, a
    /// PMT, electronics), and this is what that costs for a REASONABLE description.  No constant
    /// covers a materially wrong one: in the same measurements, leaving out 2 cm of steel behind
    /// an HPGe changed the efficiency directly behind it 30-fold at 122 keV, while a NaI PMT
    /// described as 3 mm of steel rather than its real layers was within ~1% at median but
    /// 60-80% off at 40-60 keV.  Hosts should say so when a source is behind the back plane.
    constexpr double back_structure_approx = 0.10;

    /// A response with no NearFieldModel queried inside the near-field gate: the kernel-only
    /// near-field error, i.e. the near-field boost a ray-traced kernel cannot know.
    ///
    /// MEASURED as the RMS of |exp(lnN) - 1| over the stored NearFieldModel grid at its worst
    /// distance (d/a ~ 1), pooled over three fully-characterized detectors.  Equal to
    /// `transfer_near_contact` on purpose: it is the same limitation, and fep_budget applies
    /// whichever one carries it, never both (see there for the angular structure).
    ///
    /// Gated at `near_regime_a`, not at provenance.min_distance_cm: the boost is still ~1.8-2.1%
    /// out at 2.4-3a, which a 2a gate left uncovered.
    constexpr double near_unmodeled = 0.06;

    /// Collimator shadow, by transmitted hole fraction s: below `shadow_refuse_s` the query is
    /// refuse-grade (sigma ~100%); up to `shadow_ramp_s` the sigma ramps linearly from
    /// `shadow_ramp_max` down to 0.
    ///
    /// NOT DERIVED - none of the four.  No corpus detector carries a collimator.  They are cheap
    /// to measure once points are placed by hand: s = kernel_transmitted/omega_frac_active is
    /// computable without any MC, so probes can be placed at chosen s.
    ///
    /// Known defect independent of the values: the ramp is DISCONTINUOUS at the refuse boundary,
    /// giving 0.417 just above shadow_refuse_s and 1.0 just below it.
    constexpr double shadow_refuse = 1.0;
    constexpr double shadow_refuse_s = 0.05;
    constexpr double shadow_ramp_s = 0.30;
    constexpr double shadow_ramp_max = 0.5;

    /// Model-form floor on the (ratio-only) build-up correction of a shielded eps_total, added in
    /// quadrature when a ShieldContext is supplied.
    ///
    /// NOT DERIVED.  CeeLo ships the build-up seam empty and the host installs the model, so
    /// deriving this means comparing a corrected shielded eps_total against MC through real
    /// shields with source geometry - a different experiment, and its own piece of work.
    constexpr double buildup_floor = 0.10;

    /// Far-field peak-efficiency floor - MEASURED over 36 detectors.  One constant serves every
    /// crystal class: solving per class gives 0.475% excluding CdTe and 0.356% for CdTe alone, so
    /// the class spread does not justify a split here (it does for the total).
    constexpr double fep_far_floor = 0.005;

    /// Near-field peak-efficiency floor - MEASURED, and the answer is that there is no near
    /// excess: solving inside `near_regime_a` gives 0.21% against 0.43% outside it, over the
    /// same detectors at the same generation profile.  The parameterization is if anything
    /// MORE accurate close in, because the near-field table is measured there.
    ///
    /// Set equal to `fep_far_floor` rather than to the smaller measured value - a floor that
    /// DROPS close in would be a strange promise, and the far value covers both regimes.  The
    /// peak near/far split is therefore inert by construction; the TOTAL one is not (3.0x).
    constexpr double fep_near_floor = 0.005;

    /// Far-field total-efficiency floor - MEASURED, EXCLUDING CdTe-class crystals.
    constexpr double tot_far_floor = 0.006;

    /// Far-field total-efficiency floor for CdTe/CZT-class crystals - MEASURED separately because
    /// one constant cannot serve: CdTe solves 3.2x higher than everything else.  Selected on
    /// crystal COMPOSITION (`crystal_is_cdte_class`), not on the material's name, since names are
    /// free-form user text.
    ///
    /// EMPIRICAL, and the root cause is not understood.  The excess tracks neither detector size
    /// (the corpus's smallest crystal, a bare 0.5 cm3 CZT, is its most accurate) nor `TotEffTier`
    /// (EtaTotTable spans the whole range); both were tested and rejected.  The leading remaining
    /// idea is that an uncollided kernel plus an angle-flat b(E) cannot carry scatter-in from the
    /// housing - see scratch/20260915_nonbare_det_uncert_investigate_prompt.md.  DELETE this
    /// constant once the model is fixed; it should not outlive the defect.
    constexpr double tot_far_floor_cdte = 0.020;

    /// Near-field total-efficiency floor - MEASURED.  Unlike the peak efficiency, the total does
    /// degrade close in: 1.92% inside `near_regime_a` against 0.64% outside, a factor 3.0.
    ///
    /// One value serves every class measured near-field, but that is weak evidence about CdTe: the
    /// near subset's only two CdTe detectors are the two showing no far-field excess either, and
    /// the three that motivate `tot_far_floor_cdte` were never measured near-field.
    constexpr double tot_near_floor = 0.020;

    /// Where the near regime begins, in transverse half-extents.
    ///
    /// NOT DETERMINED by the data.  Re-solving the total floors for gates from 1.5 to 5 leaves the
    /// near/far ratio at 2.0-2.1 and the achieved coverage within half a point, because the total
    /// error is a smooth ramp with distance (~1.7% at contact to ~0.6% at 4-6a), not a step.  Any
    /// gate in that range performs about equally.  A ramp would fit the physics better than a step
    /// and would need a new field.
    constexpr double near_regime_a = 4.0;

    /// Multiplier the closed loop applies to the FEP floors on a minor model-form failure.
    ///
    /// NOT DERIVED, and dead in practice: only `GenerationOptions::closed_loop` applies it, and
    /// nothing in InterSpec enables that.
    constexpr double generator_floor_inflation = 1.25;

    // --- SigmaTransferModel defaults ---------------------------------------------------------
    //
    // All four were measured on eps_FEP only (`envelope_transfer_from_mc` scores the FEP rows),
    // but common_eval applies them to eps_total as well.  There is no total-efficiency measurement
    // behind any of them.
    //
    // The measurement anchors a transfer on a MEASURED MC curve at one distance and scores it
    // against MC at every OTHER distance, angle and energy, so every scored point is held out.
    // (Anchoring on a model and re-querying the anchor measures nothing - it is a tautology.)

    /// Transfer floor on axis in the far field.  MEASURED at 0.669% RMS held out, and 0.815% at
    /// the shortest far distance probed (d/a = 5), i.e. above this value by 1.3-1.6x.  Left at
    /// 0.005 rather than raised: it is a floor under a term the anchor covariance also carries,
    /// and raising it would double-count a measured curve's own uncertainty.
    constexpr double transfer_far_onaxis = 0.005;

    /// Off-axis amplitude, mid and high energy - MEASURED.  This is the amplitude the saturating
    /// form approaches (see `transfer_offaxis_s2_half`), not a per-sin^2 slope.
    constexpr double transfer_offaxis_mid = 0.037;

    /// Extra off-axis amplitude at low energy, entering as (mid + low_e * w^2) with w the ln ramp
    /// between the two reference energies below.  MEASURED: splitting the corpus at
    /// `transfer_mid_e_ref_keV` gives amplitudes of 3.68% above and 5.23% below, a factor 1.42.
    constexpr double transfer_offaxis_low_e = 0.016;

    /// The low-energy ramp's endpoints.  NOT DERIVED - `transfer_mid_e_ref_keV` was used as the
    /// SPLIT POINT when the amplitude either side of it was measured, which says nothing about
    /// whether 45 and 150 are the right endpoints.  Testing that needs the amplitude resolved
    /// across energy rather than pooled into two bins.
    constexpr double transfer_low_e_ref_keV = 45.0;
    constexpr double transfer_mid_e_ref_keV = 150.0;

    /// Half-saturation of the off-axis term in sin^2(theta) - MEASURED at sin^2(23 deg); fitting
    /// above and below `transfer_mid_e_ref_keV` separately gives sin^2(20 deg) and sin^2(27 deg),
    /// so one value serves.
    ///
    /// The residual SATURATES with angle - past ~30 degrees you view the crystal from the side and
    /// the distribution of path lengths through it stops changing much - where sin^2(theta) grows
    /// without bound.  `A * s2/(s2 + s2_half)` beat sin^2 on HELD-OUT angles in every split tried
    /// (fit {15,45,75} predict {30,60}: RMS 0.49% vs 1.15%, and the reverse 0.75% vs 1.53%), which
    /// is the bar a change of functional form has to clear.
    ///
    /// <= 0 selects the legacy unbounded sin^2 form, which is what a stored response written
    /// before this field existed gets on read.
    ///
    /// HOW MUCH OF THE OFF-AXIS ERROR IS REDUCIBLE, since adding angular resolution is the obvious
    /// idea and buys less than it looks like: a shared shape takes the 3.08% RMS residual to
    /// 2.65%, adding an aspect-ratio amplitude to 2.50%, and a PERFECT per-(detector, angle)
    /// correction only to 2.28% - because just 45% of the variance is a fixed bias per (detector,
    /// angle) and 55% varies with ENERGY at fixed angle.  Resolving energy off axis is what
    /// collapses it, and that is a full characterization, not more anchor angles.
    ///
    /// Predictors tried and rejected, so they are not retried: crystal aspect ratio, the kernel's
    /// solid-angle-weighted mean chord, and that chord weighted by interaction probability
    /// 1 - exp(-mu(E)L).  On a common subset all three correlate at |r| ~ 0.35-0.38 and leave
    /// ~1.9%; the interaction weighting buys nothing over plain geometry despite carrying real
    /// energy dependence.  The entry chord says where a photon first interacts, while full-energy
    /// containment depends on escape FROM that point - a different geometric quantity.  See
    /// `envelope_chord_predictor`.
    constexpr double transfer_offaxis_s2_half = 0.153;

    /// Transfer near-field term at contact, ramping to zero at `transfer_near_gate_a`.
    ///
    /// MEASURED as the RMS over the full (E, cos theta) grid at the worst distance (d/a ~ 1).  The
    /// angular range matters: the error is U-shaped in cos theta - for one HPGe at d/a = 1 it runs
    /// 9.9% at grazing, a minimum of 1.6% near cos theta ~ 0.77, and 3.9% on axis - so an
    /// on-axis-only measurement gives 4.0% and under-covers grazing queries by about two.  A single
    /// number cannot carry that shape, and this is the honest one for a query of unknown angle.
    ///
    /// KNOWN SHORTCOMING: a query that is both near-grazing and close-in is still under-covered
    /// (~10% real against ~6% declared).  Carrying it needs an angular term this struct does not
    /// have, and the shape is not sin^2.
    constexpr double transfer_near_contact = 0.06;

    /// Where the transfer's near term switches off, in transverse half-extents.  With the
    /// amplitude above the ramp tracks the measured decline out to ~4.5a, where the residual is
    /// ~0.8% and the term is nearly off.
    constexpr double transfer_near_gate_a = 5.0;
}  // namespace model_sigma

// ---------------------------------------------------------------------------
// Geometry descriptor: geometry/GeometryDescriptor.h (always built). Its XML
// codec (to_xml_string / from_xml_string) is defined with this file's.
// ---------------------------------------------------------------------------

/// True when the crystal is CdTe-class (CdTe, CZT): Cd + Te carry more than half the
/// crystal's mass.  Keyed on COMPOSITION rather than on the material's name, because names
/// are free-form user text ("CZT", "CdZnTe", "Cd0.9Zn0.1Te", ...) and a name test would
/// quietly stop matching.
///
/// ONE definition, used by both response-building paths (ResponseGenerator::generate and
/// make_transfer_response).  It selects model_sigma::tot_far_floor_cdte, which is an
/// empirical patch over a modelling gap that is not understood, so this is expected to be
/// DELETED along with that constant rather than extended.
bool crystal_is_cdte_class(const GeometryDescriptor& gd);


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
///
/// "Both flanks are nodes" is the producer's job, and `finalize()` no longer
/// trusts it: a segment left with a single node is MERGED into a neighbour
/// (see finalize()), because a one-node segment would otherwise interpolate
/// as a constant and fabricate a discontinuity at the edge - for a Ge K-edge
/// falling between a file's 10 and 12 keV nodes that was a ~95x jump that
/// existed nowhere in the source data.
class EtaTable {
public:
    std::vector<double> energies_keV;   ///< ascending; edge flanks as nodes
    std::vector<double> cos_thetas;     ///< ascending, e.g. 0.02 .. 1.0
    std::vector<double> phis_deg;       ///< empty (axial) or quadrant nodes
    std::vector<double> ln_eta;         ///< [e][c][p] energy-major flattened
    std::vector<double> frac_sigma;     ///< per-node MC fractional sigma
    /// Crystal K-edges (segment breaks).  `finalize()` sorts, de-duplicates and
    ///  drops out-of-span entries in place, so a producer or a deserialized file
    ///  need not supply them ordered.
    std::vector<double> edges_keV;

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

    /// True when every cos-theta column holds the same values - an angle-flat table (curve and
    /// flat transfers), whose clamp in cos-theta loses nothing.
    bool is_angle_flat() const;

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
    std::vector<double> dists_cm;       ///< ascending, from crystal-face origin (the side &
                                        ///< back table: a + clearance, see back_near_field)
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
    EtaTotTable,   ///< HPGe-class: eps_tot = k(E) * eta_tot(E,theta) * K
    /// eps_tot is NOT CHARACTERIZED: the response carries FEP only, and the
    /// eps_total_* queries return 0 flagged NeedsMc rather than a number.
    ///
    /// This exists because KernelExact is a POSITIVE claim, not a default:
    /// ResponseGenerator only selects it after checking the bare kernel against
    /// MC to 1% (see pick_tot_tier). A producer with no total-efficiency data at
    /// all -- an FEP-only import, or a curve transfer whose source DRF had no
    /// total curve -- was previously left at the KernelExact default, so it
    /// silently served a bare-crystal kernel as if it were a verified total.
    /// For a real HPGe that is not a small error: the kernel omits the passive
    /// housing and every peak-to-total effect, and at low energy it falls BELOW
    /// the response's own eps_fep, which is physically impossible. A host gating
    /// cascade-summing on "does this response have a total?" got a confident yes
    /// and a wrong correction. Being un-representable is the honest answer, so
    /// it is a tier rather than a flag no caller has to read.
    NotCharacterized
};

struct TotEffPayload {
    TotEffTier tier = TotEffTier::KernelExact;
    std::vector<double> b_energies_keV;   ///< BCurve: ~8 log-spaced nodes
    std::vector<double> ln_b;
    EtaTable eta_tot;                     ///< EtaTotTable tier only

    /// Optional measured near-field multiplier for eps_total, applied (any characterized tier)
    /// below its own breakpoint exactly as DetectorResponse::near_field is for eps_fep: ln N_tot
    /// is ln(MC total / the tier's far-field model) at each node.  Without it every close-in total
    /// is the far-field value carried by the kernel alone - 1-4% off at grazing angles for a
    /// 5x5 cm HPGe (measured against direct MC).  The General and Contact profiles fill it from
    /// the near-field nodes' totals, which cost nothing extra since those nodes run full
    /// transport anyway; a far-field-profile or transfer response has none.  Empty = absent: not
    /// serialized, so such responses keep their bytes and content_hash.
    NearFieldModel near_field;

    /// The total's twin of DetectorResponse::back_near_field: ln(MC total / the far-field model) at
    /// the same nodes, on the same lattice and coordinate.  Empty = absent.
    NearFieldModel back_near_field;

    void finalize();
    double ln_b_at(double energy_keV) const;   ///< clamped PCHIP over (lnE, ln b)

    /// False only for #TotEffTier::NotCharacterized - i.e. whether an
    /// eps_total query returns a modeled value at all.
    bool characterized() const { return tier != TotEffTier::NotCharacterized; }

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
    struct Components { double far_onaxis = 0.0, offaxis = 0.0, near_field = 0.0; };
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

/// Per-node MC precision policy (the D0 speedup dial). Uniform reproduces the
/// historical single-target behaviour bit-for-bit; RelaxMild installs the
/// energy-graded map measured in the D0 policy memo (~2-2.5x cheaper high-E
/// nodes at a bounded p95 cost).
enum class PrecisionProfile : uint8_t {
    Uniform,    ///< node_fep_precision everywhere (default; unchanged)
    RelaxMild   ///< graded target precision by (stage, energy) -- see relax_mild_precision_map
};

const char* to_string(PrecisionProfile p);

/// Side and back incidence scanned at generation (GenerationOptions::back_coverage).
enum class BackCoverage : uint8_t { None, Coarse, Fine };

const char* to_string(BackCoverage b);

struct ResponseProvenance {
    ProductionMethod method = ProductionMethod::FullMc;
    std::string ceelo_version;        ///< library version string
    std::string created_utc;          ///< ISO-8601, informational
    ResponseProfile profile = ResponseProfile::General;
    double node_fep_precision = 0.003;///< per-node MC precision target
    /// The remaining generation choices, so a host can offer the same run again.  Set by
    /// ResponseGenerator::generate() (back_coverage on the full-MC grid, n_anchor_angles in
    /// transfer mode); empty for a response made otherwise, or saved before they were recorded.
    std::optional<PrecisionProfile> precision_profile;
    std::optional<BackCoverage> back_coverage;
    std::optional<int> n_anchor_angles;
    /// Depth (cm, behind the crystal face) of the detector body the side & back near-field tables
    /// assume - GenerationOptions::back_body_depth_cm, or 0 where it was dropped or not asked for.
    /// Used in evaluation (DetectorResponse::back_near_body), so written and hashed when > 0.
    double back_body_depth_cm = 0.0;
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

    /// Side & back near field: ln N behind the face plane (cos theta < 0), measured at generation
    /// (GenerationOptions::back_coverage with the General or Contact profile).  Same layout and
    /// interpolation as `near_field`, but its distance axis is d = a + the CLEARANCE past
    /// back_near_body() - the nearest point of the detector's outside, its back extended over the
    /// assumed body - rather than the distance from the crystal-face origin, which says little
    /// about how near a source beside a long detector is.  A source inside that body (outside the
    /// described detector, but where a real one keeps its cryostat or electronics) is left to the
    /// kernel, as on a response without the table.  phi = 0 only, as `near_field`.  Empty =
    /// absent: not serialized, so a response without it keeps its bytes and content_hash.
    NearFieldModel back_near_field;

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

    /// The outside the back near-field tables measure clearance from: the described outer
    /// envelope, its back face moved to provenance.back_body_depth_cm when that is further back.
    GeometryDescriptor::Envelope back_near_body() const;

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

    /// FNV-1a hash of the stored content - exactly what to_xml_string() saves,
    /// hashed from the members (stable content identity; grounding and grid
    /// changes change it, and a save/load round trip does not). INVARIANT to the
    /// accuracy certificate (which is metadata about the response, not
    /// content), so a response hashes the same with or without one.
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
    /// ln N_tot from tot_eff.near_field below its breakpoint (0 without one) - behind the face
    /// plane, from tot_eff.back_near_field - adding its node sigma to `ec`'s budget: the total's
    /// counterpart of fep_budget's near-field term.
    double total_near_ln(double energy_keV, EvalCommon& ec) const;
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
