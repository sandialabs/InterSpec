#ifndef CeeLoUtils_h
#define CeeLoUtils_h
/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
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

#include "InterSpec_config.h"

#include <memory>
#include <string>
#include <vector>

#include "io/DetectorResponse.h"
#include "io/EfficiencyTransfer.h"

struct Material;
struct GadrasDetectorDat;
struct AngleMaterial;
struct AngleOutxContents;
class DetectorPeakResponse;

namespace SandiaDecay{ class SandiaDecayDataBase; }

/** Small shared helpers for mapping InterSpec objects onto the CeeLo
 Monte-Carlo library's types.  Used by the detector-geometry MC UI
 (DetectorGeometryInput / MakeMcResponseForDrf), the fixed-geometry response
 builder, and the cascade-summing truth-generation tests.
 */
namespace CeeLoUtils
{
  /** The fractional 1-sigma assumed for an efficiency anchor point when the DRF carries no
   uncertainty information at all, and the floor applied when it does: a fitted curve's tiny
   statistical sigma must not claim the anchor is exact.  The only such defaults on the
   InterSpec side, kept here so they cannot multiply.

   sm_default_anchor_frac_sigma - MEASURED.  The canonical "states no uncertainty" case is a
   GADRAS Detector.dat + Efficiency.csv, and `gadras_efficiency_cross_validation` measures how
   far that stated curve sits from Monte-Carlo truth for the geometry the file describes - 160
   points over 20 detectors, 60 keV to 2 MeV, giving median |MC/stated - 1| of 2.3% (LaBr3),
   2.6% (NaI), 3.0% (CZT) and 5.8% (HPGe).  Converting those medians to an equivalent one-sigma
   (x1.48 for a half-normal) spans about 3.4% to 8.6%, and 0.05 sits inside that - generous for
   LaBr3 and NaI, about 1.7x thin for HPGe.  Note this deliberately includes GEOMETRY-description
   error, because that is part of how wrong such a DRF is, which is what an anchor sigma covers.
   Note also it is a median-to-sigma conversion, a different statistic from the RMS-pull
   calibration behind `ceelo::model_sigma`.

   sm_min_anchor_frac_sigma - NOT DERIVED, and the question is not its size.  It exists to stop a
   fitted curve's tiny statistical sigma claiming an anchor is exact; since
   `curveAnchorWithCovarianceForDrf` carries the curve's full covariance rather than a per-point
   sigma, the guard may have less to do.  Deciding that needs a look at what covariances
   MakeDrfFit actually produces, not a corpus measurement.
   */
  constexpr double sm_default_anchor_frac_sigma = 0.05;
  constexpr double sm_min_anchor_frac_sigma = 0.01;

  /** Converts an InterSpec MaterialDB material to a self-contained CeeLo
   material spec (per-element mass fractions; nuclide fractions folded into
   their element).  Throws for elements CeeLo has no data for (Z > 92).
   */
  ceelo::MaterialSpec to_ceelo_material( const Material &mat );

  /** Builds a CeeLo material straight from an ANGLE inline material definition
   (a density plus `<elements>`, `<compound>` or `<compounds>`), with no
   `MaterialDB` involved - ANGLE files routinely define their own materials, and
   those definitions are more authoritative than any name lookup.

   `<compound>` atom counts are converted to mass fractions using SandiaDecay's
   natural atomic masses, and `<compounds>` are merged mass-weighted.

   Throws std::runtime_error if `mat` has no usable composition or density, or
   names an element outside Z in [1,92] (which the Monte Carlo has no
   cross-sections for). */
  ceelo::MaterialSpec toCeeloMaterial( const AngleMaterial &mat,
                                       const SandiaDecay::SandiaDecayDataBase *db );

  /** The single-position absolute-efficiency curve a DRF pins to data, ready
   to anchor an EFFTRAN-style efficiency transfer (see
   external_libs/CeeLo/src/io/EfficiencyTransfer.h).
   */
  struct TransferAnchor
  {
    /** Absolute full-energy-peak efficiency at the reference position,
     in-vacuum, strictly ascending in energy. */
    ceelo::AnchorCurve curve;

    /** On-axis source distance of the anchor, in cm, measured from the
     geometry descriptor's reference point (user convention - typically the
     endcap front, NOT the crystal face). */
    double ref_distance_cm = 0.0;

    /** True when the anchor was sampled from the DRF's fitted efficiency
     curve rather than taken from raw measured points. */
    bool curve_derived = false;
  };//struct TransferAnchor

  /** Builds the efficiency-transfer anchor for a DRF.

   Prefers the DRF's raw measured points when they exist: points taken at
   different distances are transferred to one reference distance (the most
   common measurement distance, or `override_ref_distance_cm` when > 0) with
   the geometry kernel ratio K(E,d_ref)/K(E,d_i) - see #GeometryKernel; the
   per-point uncertainty is sqrt(stat^2 + cert^2) (the source-certificate
   correlation is conservatively flattened onto the diagonal).  Otherwise the fitted intrinsic curve is
   sampled at 16 log-spaced energies - plus flanking samples just either side
   of each crystal K-edge, since the transfer's eta interpolant is segmented at
   the edges but gets nodes only at anchor energies - and converted to absolute
   efficiency at `override_ref_distance_cm` (or, when <= 0, the distance the
   curve is pinned at: absoluteEfficiencyDistance() for far-field-absolute
   DRFs, else max( 50cm, 10 * transverse-half-extent )).

   Throws std::runtime_error for fixed-geometry DRFs, or when fewer than two
   usable anchor points can be constructed.
   */
  TransferAnchor transferAnchorForDrf(
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const ceelo::GeometryDescriptor &geom,
                      const double override_ref_distance_cm );

  /** Samples the DRF's total-efficiency curve (when it has one) at the FEP
   anchor's energies and reference distance, for use as the transfer's
   total-efficiency anchor - transferred measured totals (endcap/dead-layer
   effects included) beat the bare-crystal fallback tier for cascade-summing.
   Returns an empty curve when the DRF has no total-efficiency information.
   */
  ceelo::AnchorCurve totalTransferAnchorForDrf(
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const TransferAnchor &anchor );

  /** Packages an anchored efficiency curve as a CeeLo efficiency-transfer
   DetectorResponse: the anchor is held fixed and every (energy, position)
   query is transferred through the ray-traced kernel ratio - so evaluations at
   other distances (and mildly off-axis) are geometry-correct, with the
   attached SigmaTransferModel inflating the reported uncertainty off-axis and
   near-field.  `tot_curve` (from #totalTransferAnchorForDrf; may be empty)
   becomes the total-efficiency anchor when it has at least two points, else
   the bare-crystal kernel-exact total tier is used.

   Touches nothing but its arguments (in particular, no DetectorPeakResponse),
   so it may run on a worker thread with value-captured inputs.  No Monte
   Carlo is run - this is deterministic and sub-second.
   */
  /** Gives `drf` off-axis and near-field support, for free, when it already has everything the
   support needs: a stated geometry (#DetectorPeakResponse::geometry) and a measured efficiency to
   anchor on.  Builds the measured-curve transfer (#makeTransferResponse - no Monte Carlo, sub-
   second) and attaches it with #DetectorPeakResponse::setCeeloResponse.

   Without this, a detector imported from a `Detector.dat` plus its `Efficiency.csv` answers every
   query with the flat-disk solid angle and cannot answer off-axis at all, despite the import having
   handed us both halves of a proper transfer.

   A no-op (returning false) when the DRF already carries a response, states no geometry, is
   fixed-geometry, or cannot produce a usable anchor.  Never throws.
   */
  bool attachCurveTransferResponse( DetectorPeakResponse &drf );

  std::shared_ptr<ceelo::DetectorResponse> makeTransferResponse(
                      const ceelo::GeometryDescriptor &geom,
                      const TransferAnchor &anchor,
                      const ceelo::AnchorCurve &tot_curve,
                      const std::string &detector_name );

  /** MC-free full-energy-peak kernel K(E, position) of a detector geometry - the attenuation-weighted
   effective solid angle the EFFTRAN transfer is built on (see external_libs/CeeLo/src/io/EfficiencyTransfer.h).

   Used to put measurements taken at different source distances onto one footing: the absolute
   efficiency at distance d is eta(E)*K(E,d), so dividing a measured point by
   #intrinsicFactor(E, d) gives the same far-field intrinsic efficiency that
   DetectorPeakResponse::intrinsicEfficiencyEval reports for a geometry-modeled DRF, whatever
   distance the point was taken at.  Sub-millisecond per evaluation; the ray quadrature is built once
   per distinct distance.  Not thread-safe (build one per thread); the geometry is copied in.
   */
  class GeometryKernel
  {
  public:
    explicit GeometryKernel( const ceelo::GeometryDescriptor &geom, int n_rays = 2048 );
    ~GeometryKernel();

    /** K(E) for an on-axis point source `face_dist_cm` from the detector face. */
    double kernel( double energy_keV, double face_dist_cm );

    /** The on-axis distance (cm, from the face) the far-field intrinsic efficiency is defined at:
     max(1000*a, 100 cm), with `a` the transverse half-extent - as intrinsicEfficiencyEval uses. */
    double farFieldDistanceCm() const;

    /** g(E,d) = K(E,d) * Omega_disk(d_far,a) / K(E,d_far): absolute efficiency at `face_dist_cm`
     divided by this is the far-field intrinsic efficiency.  Reduces to the flat-disk solid angle
     fraction far from the detector. */
    double intrinsicFactor( double energy_keV, double face_dist_cm );

    /** d ln(g)/d(d) at `face_dist_cm`, per cm (central difference) - about -2/d far from the
     detector; multiply by a distance uncertainty (cm) for the fractional efficiency uncertainty
     it induces. */
    double intrinsicFactorDistanceSlope( double energy_keV, double face_dist_cm );

  private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
  };//class GeometryKernel

  /** Convenience over #GeometryKernel::intrinsicFactor for parallel vectors of energies (keV) and
   face distances (cm). */
  std::vector<double> farFieldIntrinsicFactors( const ceelo::GeometryDescriptor &geom,
                                                const std::vector<double> &energies_keV,
                                                const std::vector<double> &face_dist_cm,
                                                int n_rays = 2048 );

  /** Convenience over #GeometryKernel::intrinsicFactorDistanceSlope (per cm). */
  std::vector<double> farFieldIntrinsicFactorDistanceSlopes( const ceelo::GeometryDescriptor &geom,
                                                const std::vector<double> &energies_keV,
                                                const std::vector<double> &face_dist_cm,
                                                int n_rays = 2048 );

  /** The transfer anchor a fitted DRF should carry: the DRF's efficiency CURVE, densely sampled
   (24 log-spaced energies plus flanks either side of each crystal K-edge), converted to absolute
   efficiency at `ref_distance_cm` (or, when <= 0, at #GeometryKernel::farFieldDistanceCm) with the
   geometry kernel, together with the curve's full fractional covariance (a fitted equation's
   coefficient covariance propagated to the sample energies, or its node covariance) as
   `AnchorCurve::frac_cov` - so correlations between energies survive into the response - and the
   raw measured points, when the DRF has them, as provenance.

   Unlike #transferAnchorForDrf this never anchors on the raw points themselves: a smooth fitted
   curve is the better estimate of the efficiency between measured energies, and its covariance
   is the right uncertainty for it.  Reads the DRF's curve directly, so it is unaffected by any
   response the DRF may already carry.  `curve_derived` is true in the returned anchor.

   Throws std::runtime_error for an invalid or fixed-geometry DRF.
   */
  TransferAnchor curveAnchorWithCovarianceForDrf(
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const ceelo::GeometryDescriptor &geom,
                      const double ref_distance_cm );

  /** InterSpec measures every source-to-detector distance from the DETECTOR FACE - the endcap
   front, the surface a user can put a ruler against.  That is the only convention the user ever
   sees.  CeeLo's `GeometryDescriptor::reference_point` is a convention internal to the library that
   InterSpec never consults: every query position is formed HERE, face-referenced, and handed to the
   position-taking CeeLo API (`eps_fep_at`, `eps_total_at`, `make_quadrature`, the position forms of
   `frac_covariance`).  Nothing in InterSpec may call the distance-taking `eps_fep` / `eps_total` /
   `frac_covariance(energies, theta, phi, dist)` / `query_position` / `reference_point_position`,
   which would silently re-interpret the distance in whatever convention
   the descriptor happens to carry (a CrystalFace descriptor would put the source one endcap offset
   too close).  The flat-disk model lives in the same convention: `DetectorPeakResponse::efficiency`
   adds the detector setback to the face distance to reach the crystal.

   Returns the source position in CeeLo's crystal-face frame (cm; z = 0 at the crystal face, the
   source in front at negative z) for a source `dist_from_face_cm` from the detector face along the
   (theta, phi) direction, theta measured from the detector axis.
   */
  Eigen::Vector3d sourcePositionFromFace( const ceelo::GeometryDescriptor &gd,
                                          double theta_rad, double phi_rad,
                                          double dist_from_face_cm );

  /** THE far-field on-axis query InterSpec's intrinsic (solid-angle-free) views of a response use:
   DetectorPeakResponse::intrinsicEfficiencyEval, its no-geometry efficiencyFracCovariance, and
   setLegacyEfficiencyFromResponse all evaluate at this one position, so the stored curve, the
   live query and its covariance cannot come from different places.  The distance from the face is
   max(1000*a, 100 cm), `a` the transverse half-extent - well beyond every near-field term.
   */
  double farFieldDistanceCm( const ceelo::GeometryDescriptor &gd );
  Eigen::Vector3d farFieldSourcePosition( const ceelo::GeometryDescriptor &gd );

  /** The detector face (endcap front) in the crystal-face frame, (0, 0, -endcap_front_offset_cm):
   the point InterSpec's distances are measured to, and the point a source-side ray aims at.
   */
  Eigen::Vector3d detectorFacePosition( const ceelo::GeometryDescriptor &gd );

  /** Converts a distance CeeLo quotes from the crystal-face origin (`Provenance::min_distance_cm`,
   `EvalCommon::d_cm`, ...) into the on-axis distance from the detector face that InterSpec shows the
   user; clamped at zero.
   */
  double faceDistanceFromCrystalOrigin( const ceelo::GeometryDescriptor &gd, double d_crystal_cm );

  /** Maps the physical detector model parsed from an ANGLE file (see
   InterSpec/AngleOutxImport.h) onto a CeeLo geometry descriptor - used to seed
   the Make MC Response tool ("generic detector" import mode).

   Materials are resolved in three steps: an inline ANGLE definition is used
   directly (see #toCeeloMaterial - no `MaterialDB` needed); otherwise the name
   is mapped through ANGLE's 19 predefined materials onto a `MaterialDB` entry,
   with a chemical-formula fallback; anything still unresolved becomes a
   near-transparent spacer that preserves the layer's physical extent but not
   its attenuation, and a human-readable note is appended to `warnings`.  A thin
   implanted contact is folded into the dead layer; `contactPin` and any piece
   behind the crystal are dropped with a warning.  The crystal's front-edge bulletization and
   the core's rounded tip are carried through (`bullet_radius_cm` /
   `BoreHoleConfig::rounded_tip`), and dropped with a warning only when they
   violate a CeeLo geometry precondition (see #relaxGeometryFeatures).  The
   vacuum gap becomes a near-transparent spacer layer, so the crystal's recess
   behind the endcap is preserved.  Lengths in `contents` are PhysicalUnits and
   converted to cm here.  `reference_point` is set to EndcapFront (ANGLE's
   convention).

   Throws std::runtime_error when `contents.hasGeometry` is false or the
   crystal dimensions are non-positive.
   */
  ceelo::GeometryDescriptor buildAngleGeometry( const AngleOutxContents &contents,
                                                std::vector<std::string> &warnings );

  /** Maps a `ceelo::GeometryProblem` onto the MakeMcResponseForDrf.xml message
   id ("dgi-err-...") that explains it to the user.

   Kept here rather than in CeeLo so the library stays free of InterSpec and Wt
   string resources; `ceelo::to_string()` is the developer-facing counterpart.
   */
  const char *geometryProblemMsgId( const ceelo::GeometryProblem problem );

  /** Drops the optional CeeLo crystal-shape features `gd` cannot legally carry
   - the front-edge fillet and/or the rounded bore tip - so `build_geometry()`
   and `configure_calculator()` cannot trip a (debug-only) precondition and
   trace garbage in a release build.  Each drop appends a note to `warnings`.

   Features are dropped, never clamped: a clamped fillet is a silently wrong
   number, whereas a dropped one is the sharp-edged approximation InterSpec
   shipped before the fillet was modeled at all.  The crystal, bore and dead
   layer are never modified - those are the file's primary content, and a
   descriptor still reporting problems() afterwards is a hard error for the
   caller.

   For IMPORTED geometry only.  User-entered geometry is validated in
   DetectorGeometryInput::toDescriptor(), which rejects rather than silently
   rewrites what was typed.
   */
  void relaxGeometryFeatures( ceelo::GeometryDescriptor &gd,
                              std::vector<std::string> &warnings );

  /** Builds a far-field-absolute seed DRF from an ANGLE file's measured
   reference efficiency points, pinned at `contents.referenceDistanceCm`.  The
   points are also attached as raw measured points, so #transferAnchorForDrf
   can anchor an EFFTRAN-style transfer directly on them.

   Throws std::runtime_error when fewer than two positive reference points are
   present, or the crystal diameter / reference distance are non-positive.
   */
  std::shared_ptr<DetectorPeakResponse> buildAngleSeedDrf(
                      const AngleOutxContents &contents );


  //=========================== GADRAS Detector.dat ===========================
  // A GADRAS `Detector.dat` describes the same thing an ANGLE file does - a
  //  physical detector - so it maps onto a CeeLo geometry the same way (see
  //  #buildAngleGeometry).  The differences are all in what GADRAS does NOT
  //  record: no bore geometry, attenuators given only as an effective atomic
  //  number plus an areal density, and (usually) no crystal material at all.

  /** Builds a `ceelo::MaterialSpec` from a GADRAS stoichiometric formula, e.g.
   "Bi4Ge3O12", "CsI", "Cd0.96Zn0.04Te", "C14[H2]12".

   The subscripts are ATOM COUNTS (fractional counts allowed; a missing count
   means 1), converted to mass fractions through the natural atomic masses.
   `[H2]`-style isotope brackets are folded onto their element.

   Deliberately NOT `MaterialDB::materialFromChemicalFormula`, whose
   `Material::parseChemicalFormula` treats the subscripts as MASS FRACTIONS and
   requires an explicit count on every element - it would read "Bi4Ge3O12" as
   masses 4:3:12 and would throw outright on "CsI".

   Throws std::runtime_error on an unparseable formula, an unknown element, a
   non-positive density, or an element the Monte Carlo has no data for (Z > 92).
   */
  ceelo::MaterialSpec materialFromGadrasFormula( const std::string &formula,
                                                 const double density_g_per_cm3,
                                                 const std::string &name );

  /** The crystal material for a GADRAS material-table name (e.g. "NaI", "BGO").

   Prefers the CeeLo built-ins - whose names match the GADRAS table spellings
   exactly - so the geometry form's crystal combo matches by `.name`; otherwise
   builds from the table's formula and density via #materialFromGadrasFormula.

   Throws std::runtime_error when the name is empty or not in the table.
   */
  ceelo::MaterialSpec gadrasCrystalMaterial( const std::string &gadras_material_name );

  /** A stand-in material for a GADRAS attenuator, which is specified only as an
   effective atomic number and an areal density (g/cm2) - never a composition or
   a thickness.

   CeeLo has no generic areal-density attenuator (`Geometry::add_attenuator`
   takes a real material plus thicknesses), so one is synthesized: the two
   elements bracketing `atomic_number`, mass-weighted so the mixture's mass
   attenuation coefficient interpolates between them, at the density
   `areal_density_g_cm2 / thickness_cm` that makes the layer's areal density
   come out exactly right.

   Mass-weighting is what makes this correct: a mixture's mu/rho is the
   mass-weighted sum of its components' mu/rho, which is the same quantity
   InterSpec's own `MassAttenuation::massAttenuationCoefficientFracAN`
   interpolates across atomic number.

   Throws std::runtime_error unless `areal_density_g_cm2` and `thickness_cm` are
   both positive and `atomic_number` is within [1,92].
   */
  ceelo::MaterialSpec genericAttenuatorMaterial( const double atomic_number,
                                                 const double areal_density_g_cm2,
                                                 const double thickness_cm );

  /** Maps a parsed GADRAS `Detector.dat` onto a CeeLo geometry descriptor - the
   counterpart of #buildAngleGeometry, and the basis of the "generic detector"
   import mode.

   The crystal material comes from `GadrasDetectorDat::materialName()`; when the
   file names none, a stand-in is used and noted in `warnings` rather than
   failing the import - the geometry form lets the user set it.

   The crystal shape and dimensions come from `GadrasDetectorDat::inferShape()`,
   the dead layer from param 51 (applied to the front AND sides - the file gives
   one thickness with no direction), and the inner/outer attenuators become
   synthesized generic layers (#genericAttenuatorMaterial) whose thicknesses
   divide the GADRAS setback in proportion to their areal densities, so that
   `endcap_front_offset_cm()` reproduces `setbackCm()` exactly and
   `reference_point` can be EndcapFront - GADRAS's own convention.

   Everything GADRAS cannot express is dropped with a note in `warnings`, never
   silently defaulted: the coaxial bore (absent from the format), the side
   shield when its per-face coverage is partial (CeeLo layers are concentric,
   with no azimuthal coverage), and the back shield (CeeLo models no attenuator
   behind the crystal).

   Throws std::runtime_error when the crystal dimensions are non-positive (e.g.
   the shipped "HPGe 40%", whose length is 0) or the crystal material cannot be
   resolved.
   */
  ceelo::GeometryDescriptor buildGadrasGeometry( const GadrasDetectorDat &dat,
                                                 std::vector<std::string> &warnings );


  //====================== MC response -> legacy DRF curve =====================

  /** Fills @p drf's ordinary (non-CeeLo) intrinsic efficiency curve by sampling
   @p response - the "backbone" efficiency points a Monte-Carlo characterization
   produces.

   Needed because a DRF built from geometry alone has no measured curve of its
   own: without this it carries the placeholder `setIntrinsicEfficiencyFormula("1.0")`
   and is not valid, serializable or exportable on its own, even though every
   #DetectorPeakResponse::EffEval query would have answered correctly through the
   attached response.

   Sampling reproduces `intrinsicEfficiencyEval`'s CeeLo branch exactly - on
   axis, at the response's own far field, divided by the same disk solid angle -
   so the curve and the CeeLo dispatch cannot disagree.  Points are log-spaced
   over the response's validated energy range, plus a pair either side of each
   crystal K-edge so the interpolated curve never bridges one.  The per-point MC
   sigma is carried into `EnergyEffPoint::efficiencyUncert`.

   The detector diameter, geometry type and energy range are taken from the
   response; @p drf's setback and DRF-source are preserved across the call (they
   are members #DetectorPeakResponse::setEfficiencyPoints otherwise resets).

   Throws std::runtime_error for a null response, or if fewer than two positive
   points can be sampled.
   */
  void setLegacyEfficiencyFromResponse( DetectorPeakResponse &drf,
                      const std::shared_ptr<const ceelo::DetectorResponse> &response,
                      const size_t num_points = 48 );

  /** The canonical text for a generic attenuator - one specified only by an
   effective atomic number and an areal density, as GADRAS does.  Used as the
   material name in a #genericAttenuatorMaterial and shown in the geometry
   form's layer rows, where it reads as an editable value rather than a material
   name nothing can resolve.
   */
  std::string genericAttenuatorName( const double atomic_number,
                                     const double areal_density_g_cm2 );

  /** Reads #genericAttenuatorName back.  Returns false (leaving the outputs
   untouched) when @p text is not in that form, so a caller can fall through to
   resolving it as an ordinary material name.
   */
  bool parseGenericAttenuatorName( const std::string &text, double &atomic_number,
                                   double &areal_density_g_cm2 );

  /** Whether a layer-material text means "nothing there": empty, or one of the
   spellings of a vacuum gap ("void", "vacuum", "galactic", "galactic vacuum",
   "none"), compared trimmed and case-insensitively.  The geometry form treats
   such a layer as a gap that keeps its thickness.
   */
  bool isVacuumMaterialName( const std::string &name );

  /** The near-transparent spacer that stands in for a vacuum gap (rho ~ 1e-25
   g/cm3, single H): preserves a layer's physical extent - and so the crystal
   recess - while attenuating ~nothing.  MaterialDB-independent; the same
   spacer the ANGLE import uses for cryostat gaps.  Named "vacuum", which
   #isVacuumMaterialName recognizes, so it round-trips through the geometry form.
   */
  ceelo::MaterialSpec vacuumMaterialSpec();
}//namespace CeeLoUtils

#endif //CeeLoUtils_h
