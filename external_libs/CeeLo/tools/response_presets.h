#ifndef CEELO_TOOLS_RESPONSE_PRESETS_H
#define CEELO_TOOLS_RESPONSE_PRESETS_H
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

// Shared detector-geometry presets for golden fixtures, report drivers, and
// benchmark tools (single source of truth; formerly private to
// make_golden_response.cpp).

#include "io/DetectorResponse.h"
#include "materials/Material.h"

#include <string>
#include <stdexcept>

namespace ceelo {

inline GeometryDescriptor preset_descriptor(const std::string& name) {
    GeometryDescriptor gd;
    if (name == "nai3x3") {
        // Campaign Z2 sentinel: NaI 3"x3" + 0.5 mm Al can.
        gd.shape = DetectorShape::Cylinder;
        gd.dimensions_cm = {3.81, 7.62};
        gd.materials = {MaterialSpec::from(make_NaI()),
                        MaterialSpec::from(make_Aluminum())};
        gd.crystal_material_index = 0;
        LayerSpec can;
        can.material_index = 1;
        can.front_thickness_cm = 0.05;
        can.side_thickness_cm = 0.05;
        can.z_end_cm = 7.62;
        gd.layers.push_back(can);
    } else if (name == "hpge_coax") {
        // Campaign Z5 sentinel: large coax HPGe, bore + Li dead layer + Al endcap.
        gd.shape = DetectorShape::Cylinder;
        gd.dimensions_cm = {4.0, 8.0};
        gd.materials = {MaterialSpec::from(make_HPGe()),
                        MaterialSpec::from(make_Aluminum())};
        gd.crystal_material_index = 0;
        gd.bore = BoreHoleConfig{0.5, 6.5};
        gd.dead_layer = DeadLayerConfig{0.07, 0.07, 0.0};
        LayerSpec endcap;
        endcap.material_index = 1;
        endcap.front_thickness_cm = 0.15;
        endcap.side_thickness_cm = 0.10;
        endcap.z_end_cm = 8.0;
        gd.layers.push_back(endcap);
        gd.reference_point = ReferencePoint::EndcapFront;
    } else if (name == "czt_box") {
        // Campaign Z6: bare CZT 1x1x0.5 cm box (phi axis).
        gd.shape = DetectorShape::Box;
        gd.dimensions_cm = {0.5, 0.5, 0.5};
        gd.materials = {MaterialSpec::from(make_CZT())};
        gd.crystal_material_index = 0;
        gd.symmetry = ResponseSymmetry::Quadrant;
    } else if (name == "detective_x") {
        // ORTEC Detective-X-like coax HPGe, from PUBLIC specs (approximate:
        // ~65 mm dia x 50 mm P-type coax, Al endcap, ~0.7 mm Li dead layer).
        gd.shape = DetectorShape::Cylinder;
        gd.dimensions_cm = {3.25, 5.0};
        gd.materials = {MaterialSpec::from(make_HPGe()),
                        MaterialSpec::from(make_Aluminum())};
        gd.crystal_material_index = 0;
        gd.bore = BoreHoleConfig{0.45, 3.8};
        gd.dead_layer = DeadLayerConfig{0.07, 0.07, 0.0};
        LayerSpec endcap;
        endcap.material_index = 1;
        endcap.front_thickness_cm = 0.15;
        endcap.side_thickness_cm = 0.10;
        endcap.z_end_cm = 5.0;
        gd.layers.push_back(endcap);
        gd.reference_point = ReferencePoint::EndcapFront;
    } else {
        throw std::runtime_error("unknown preset: " + name);
    }
    return gd;
}

}  // namespace ceelo

#endif  // CEELO_TOOLS_RESPONSE_PRESETS_H
