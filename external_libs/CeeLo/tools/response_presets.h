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
        gd.set_dimensions(CylinderDims{3.81, 7.62});
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
        gd.set_dimensions(CylinderDims{4.0, 8.0});
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
        gd.set_dimensions(BoxDims{0.5, 0.5, 0.5});
        gd.materials = {MaterialSpec::from(make_CZT())};
        gd.crystal_material_index = 0;
        gd.symmetry = ResponseSymmetry::Quadrant;
    } else if (name == "detective_x") {
        // ORTEC Detective-X-like coax HPGe, from PUBLIC specs (approximate:
        // ~65 mm dia x 50 mm P-type coax, Al endcap, ~0.7 mm Li dead layer).
        gd.set_dimensions(CylinderDims{3.25, 5.0});
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
    } else if (name == "micro_detective") {
        // ORTEC Micro-Detective as used by JRC for the open U/Pu data: HPGe 50 mm dia x 30 mm (the JRC paper's
        //  length; with it a fixed CeeLo efficiency reproduces the JRC CBNM enrichments), 8 mm bullet, rounded bore
        //  8.7 x 21.3 mm (guessed), 0.7 mm dead layers, 0.5 mm PET over the face, mount cup 1.5 mm Al + 0.38 mm
        //  stainless, 1 mm Al endcap whose outer face is 1.0 cm in front of the crystal (a vacuum shell fills the gap).
        const double len = 3.0;
        gd.set_dimensions(CylinderDims{2.5, len});
        gd.bullet_radius_cm = 0.8;
        gd.materials = {MaterialSpec::from(make_HPGe()),
                        MaterialSpec::from(Material("PET", 1.38, {{1, 0.041959}, {6, 0.625017}, {8, 0.333025}})),
                        MaterialSpec::from(make_Aluminum()),
                        MaterialSpec::from(make_StainlessSteel304()),
                        MaterialSpec::from(Material("Vacuum", 1.0E-9, {{7, 1.0}}))};
        gd.crystal_material_index = 0;
        gd.bore = BoreHoleConfig{0.435, 2.13, true};
        gd.dead_layer = DeadLayerConfig{0.07, 0.07, 0.0};
        const auto layer = [](int mat, double front, double side, double z0, double z1) {
            LayerSpec l;
            l.material_index = mat;
            l.front_thickness_cm = front;
            l.side_thickness_cm = side;
            l.z_start_cm = z0;
            l.z_end_cm = z1;
            return l;
        };
        gd.layers.push_back(layer(1, 0.05, 0.0, 0.0, len));             // PET over the face
        gd.layers.push_back(layer(2, 0.0, 0.15, 0.0, 4.5));             // mount cup, Al
        gd.layers.push_back(layer(3, 0.0, 0.038, 0.0, 4.5));            // mount cup, stainless
        gd.layers.push_back(layer(4, 1.0 - 0.05 - 0.1, 0.4, 0.0, 4.5)); // vacuum gap to the endcap
        gd.layers.push_back(layer(2, 0.1, 0.1, 0.0, 5.0));              // endcap
        gd.reference_point = ReferencePoint::EndcapFront;
    } else {
        throw std::runtime_error("unknown preset: " + name);
    }
    return gd;
}

}  // namespace ceelo

#endif  // CEELO_TOOLS_RESPONSE_PRESETS_H
