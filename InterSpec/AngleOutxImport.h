#ifndef AngleOutxImport_h
#define AngleOutxImport_h
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

#include <string>
#include <vector>
#include <memory>
#include <utility>

class DetectorPeakResponse;

/** Plain-data result of parsing an ANGLE `.outx`/`.xml` efficiency file.

 ANGLE files carry far more than the geometry-transferred efficiency curve
 (`<results>`) InterSpec historically imported: a full physical detector model
 (`<detector>`) and a measured reference efficiency curve
 (`<referenceEfficiencyCurve>`).  This struct holds the raw contents so callers
 can either (Mode B) use `fixedGeomDrf` as today, or (Mode A) map the geometry
 onto a `ceelo::GeometryDescriptor` and anchor an efficiency transfer / MC on
 the reference curve.

 All lengths are stored in InterSpec `PhysicalUnits` units (the file's mm/cm
 values multiplied by `PhysicalUnits::mm`/`::cm`); downstream code divides by
 `PhysicalUnits::cm` when building a CeeLo descriptor.  Materials are kept by
 name only, so the parser needs no `MaterialDB` (name resolution happens in the
 UI layer).  Point/layer/region counts are all left generic - nothing is
 assumed about how many of each a given file contains.
 */
struct AngleOutxContents
{
  /** The fixed-geometry DRF built from `<results>` - the historical import
   product (Mode B).  Always populated on success. */
  std::shared_ptr<DetectorPeakResponse> fixedGeomDrf;

  // --- detector geometry (PhysicalUnits length units), materials by name ---
  std::string detName;
  std::string detDescription;
  std::string detType;

  double crystalRadius = 0.0;
  double crystalLength = 0.0;
  double bulletizingRadius = 0.0;

  bool hasCore = false;
  double coreRadius = 0.0;
  double coreDepth = 0.0;

  double deadLayerFront = 0.0;
  double deadLayerSide = 0.0;

  std::string contactMaterial;
  double contactSide = 0.0;

  /** Vacuum gap in front of / beside the crystal - a spacing, not an
   attenuator (attenuation ~ 0). */
  double vacuumFront = 0.0;
  double vacuumSide = 0.0;

  /** One concentric endcap/window/housing layer.  `front`/`side` are the
   thicknesses (PhysicalUnits units); `material` is the ANGLE material name. */
  struct Layer
  {
    std::string material;
    double front = 0.0;
    double side = 0.0;
  };//struct Layer

  /** Endcap, window, and housing pieces, ordered innermost -> outermost. */
  std::vector<Layer> layers;

  // --- reference anchor ---
  /** Measured reference efficiency points: (energy in keV, absolute
   efficiency), ascending in energy. */
  std::vector<std::pair<float,float>> referencePoints;

  /** Reference on-axis source distance, in cm (derived from the reference
   curve's `<holder height>`); approximate, meant to be user-editable. */
  double referenceDistanceCm = 0.0;

  // --- sample source geometry (the top-level `<geometry>`/`<source>` block) ---
  /** These describe the SAMPLE the file's `<results>` efficiency is for (not the
   reference curve).  They are not used to build the DRF, but let callers (e.g.
   cross-validation tests, or a future "seed the Activity/Shielding source
   geometry" link) reconstruct the measured source.  Zero/empty when absent. */

  /** On-axis standoff from the endcap front (i.e. the outer detector face) to
   the sample's NEAR (bottom) face - where the sample begins, not its center.
   This is the empty gap between detector and sample, taken from the sample
   holder's `<height>` (the sample rests on top of the holder).  For an extended
   source, add half of `sourceHeight` to reach the source center; for a
   point/thin source the two coincide.  In cm. */
  double sampleDistanceCm = 0.0;

  /** Sample source radius / height, in cm (0 for a point/unspecified source). */
  double sourceRadius = 0.0;
  double sourceHeight = 0.0;

  /** Sample source matrix material: ANGLE trade name and (if given) density. */
  std::string sourceMaterialName;
  double sourceDensity = 0.0;

  /** One element of the sample source material composition. `massFraction` is
   a fraction in [0,1] (ANGLE mass-percent divided by 100). */
  struct SourceElement
  {
    std::string symbol;
    double massFraction = 0.0;
  };//struct SourceElement

  std::vector<SourceElement> sourceComposition;

  // --- fitted region equations (developer-checks only) ---
  /** One fitted log-log efficiency region.  `equation` is the raw text of the
   XML comment ANGLE writes after each `<region>` (e.g.
   "log10 e = c0 + c1 * log10 Eg + c2 * (log10 Eg)^2 + ...").  Extracted only
   for a self-consistency dev-check; not interpreted otherwise. */
  struct FitRegion
  {
    double startKeV = 0.0;
    double endKeV = 0.0;
    std::string equation;
  };//struct FitRegion

  std::vector<FitRegion> fitRegions;

  /** Whether a usable detector geometry / reference curve was parsed; gate
   Mode A availability on these. */
  bool hasGeometry = false;
  bool hasReference = false;
};//struct AngleOutxContents

#endif //AngleOutxImport_h
