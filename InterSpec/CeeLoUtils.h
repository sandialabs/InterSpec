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
struct AngleOutxContents;
class DetectorPeakResponse;

/** Small shared helpers for mapping InterSpec objects onto the CeeLo
 Monte-Carlo library's types.  Used by the detector-geometry MC UI
 (DetectorGeometryInput / MakeMcResponseForDrf), the fixed-geometry response
 builder, and the cascade-summing truth-generation tests.
 */
namespace CeeLoUtils
{
  /** Converts an InterSpec MaterialDB material to a self-contained CeeLo
   material spec (per-element mass fractions; nuclide fractions folded into
   their element).  Throws for elements CeeLo has no data for (Z > 92).
   */
  ceelo::MaterialSpec to_ceelo_material( const Material &mat );

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

   Prefers the DRF's raw measured points when they exist and were all taken at
   a single distance (within 1%); the per-point uncertainty is then
   sqrt(stat^2 + cert^2) (the source-certificate correlation is conservatively
   flattened onto the diagonal).  Otherwise the fitted intrinsic curve is
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
  std::shared_ptr<ceelo::DetectorResponse> makeTransferResponse(
                      const ceelo::GeometryDescriptor &geom,
                      const TransferAnchor &anchor,
                      const ceelo::AnchorCurve &tot_curve,
                      const std::string &detector_name );

  /** Maps the physical detector model parsed from an ANGLE file (see
   InterSpec/AngleOutxImport.h) onto a CeeLo geometry descriptor - used to seed
   the Make MC Response tool ("generic detector" import mode).

   Materials are resolved via the session `MaterialDB` singleton, with
   chemical-formula fallbacks for ANGLE trade names (Mylar/PET, Brass, Carbon
   fiber); anything still unresolved is skipped and a human-readable note
   appended to `warnings`.  The germanium contact is folded into the side dead
   layer; `contactPin` is dropped.  The crystal's front-edge bulletization and
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
}//namespace CeeLoUtils

#endif //CeeLoUtils_h
