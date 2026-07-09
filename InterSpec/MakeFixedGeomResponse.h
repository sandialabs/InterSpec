#ifndef MakeFixedGeomResponse_h
#define MakeFixedGeomResponse_h
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

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <functional>

#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/GammaInteractionCalc.h"

class DetectorPeakResponse;

/** Computes a FIXED-GEOMETRY detector response for one specific
 source/shielding scene via CeeLo Monte-Carlo (EfficiencyCalculator over the
 assembled scene - NOT the angular ResponseGenerator characterization), and
 packages it as a DetectorPeakResponse:
   - full-energy-peak + total efficiency-per-emitted-photon curves for the
     whole scene (source distribution, source matrix, shields, air),
   - per-node MC uncertainties as the DRF efficiency uncertainty, and
   - the scene itself embedded as the DRF's fixed-geometry setup blob
     (DetectorPeakResponse::fixedGeometrySetupXml), which the Act/Shield tool
     then displays read-only.

 Only available when nothing geometric is being fit and the current DRF has a
 CeeLo MC response attached (that response's GeometryDescriptor provides the
 detector model the scene is built around).
 */
class MakeFixedGeomResponse
{
public:
  /** The scene description embedded in the produced DRF.  Owned by the
   Act/Shield layer (DetectorPeakResponse treats the XML as opaque).
   */
  struct Setup
  {
    GammaInteractionCalc::GeometryType geometry
                              = GammaInteractionCalc::GeometryType::Spherical;

    /** Detector face to the CENTER of the geometry (PhysicalUnits length). */
    double distance = 0.0;

    std::vector<ShieldingSourceFitCalc::ShieldingInfo> shieldings;

    std::string toXmlString() const;

    /** Throws std::runtime_error on malformed input. */
    void fromXmlString( const std::string &xml );
  };//struct Setup

  /** Whether the scene can be expressed in a CeeLo scene: no generic (AN/AD)
   layers, the source (self-atten/trace host, if any) is the innermost layer,
   and the geometry type maps onto a CeeLo source shape.  If not, `reason`
   (when non-null) is set to a human-readable explanation.
   */
  static bool sceneRepresentable( const Setup &setup, std::string *reason );

  /** Runs the Monte-Carlo (synchronously - call from a worker thread) and
   builds the fixed-geometry DRF.

   @param base_drf The current detector; must have a CeeLo response attached
          (its descriptor is the detector model).  Name/FWHM carry over; the
          far-field curves, MC response, and measured points are replaced by
          the scene-specific fixed-geometry curves.
   @param setup The scene.
   @param extra_energies_keV Extra energy nodes (e.g. the fit peaks' gamma
          lines) unioned with a log grid over the valid range.
   @param fep_precision Per-node MC fractional precision target (e.g. 0.005).
   @param progress Called with fraction complete (0..1); may be empty.
   @param cancel When set true, throws std::runtime_error("cancelled").

   Throws std::runtime_error on invalid input, cancellation, or MC failure.
   */
  static std::shared_ptr<DetectorPeakResponse> computeFixedGeomDrf(
                    const std::shared_ptr<const DetectorPeakResponse> &base_drf,
                    const Setup &setup,
                    const std::vector<double> &extra_energies_keV,
                    const double fep_precision,
                    const std::function<void(double)> &progress,
                    const std::shared_ptr<std::atomic<bool>> &cancel );
};//class MakeFixedGeomResponse

#endif //MakeFixedGeomResponse_h
