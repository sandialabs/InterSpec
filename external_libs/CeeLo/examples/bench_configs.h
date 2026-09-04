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

/// @file bench_configs.h
/// @brief Shared builders for the GEANT4-validated benchmark configurations.
///
/// Single source of truth for configs 1, 2, 3, 5, 6, 7 (point sources) and
/// 8, 11, 12 (source-effect configs) so that benchmark and validation tools
/// run exactly the geometries documented in DESIGN.md.

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace bench {

/// Owns the materials referenced by a configured EfficiencyCalculator.
/// EfficiencyCalculator stores raw `const Material*` pointers, so the
/// materials must outlive the calculator.
struct ConfigSetup {
    std::vector<std::unique_ptr<ceelo::Material>> materials;
    std::string description;
};

inline std::vector<double> default_energies_for_cfg(int cfg) {
    switch (cfg) {
        case 1:
        case 2:
        case 3:
            return {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1173, 1332, 2000, 3000};
        case 5:
            return {30, 50, 80, 100, 200, 300, 500, 662, 800, 1000, 1500};
        case 6:
            return {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000};
        case 7:
            return {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000};
        case 8:
            return {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614};
        case 11:
            return {100, 200, 300, 500, 662, 1000, 1173, 1332, 2000, 3000};
        case 12:
            return {59, 100, 200, 300, 500, 662, 1000, 1173, 1332, 2614};
        // Out-of-sample stress configs (13-17, 19): used to validate the
        // biasing auto-enable policy outside the GEANT4-anchored set. Not
        // gated by compare_validation.py (no G4 references).
        case 13:
            return {59, 662};
        case 14:
            return {662, 1332, 3000};
        case 15:
            return {186, 662, 1460};
        case 16:
            return {122, 662};
        case 17:
            return {662, 1332};
        case 19:
            return {300, 662};
        // HPGe sharp/bulletized pair: weighted low, where the front-edge
        // geometry matters most.
        case 25:
        case 26:
            return {45, 59.5, 88, 122, 344, 662, 1332};
        // Deep-Fe-shield low-energy pair (2026-09-04): where the analytic
        // Rayleigh-deflection loss lives, and where CeeLo's own scattered-FEP
        // stream had never been checked against GEANT4 (cfg 11 starts at 200 keV).
        case 27:
        case 28:
            return {60, 88, 122};
        default:
            return {};
    }
}

/// Configure `calc` for benchmark config `cfg`. Materials are stored in
/// `setup.materials` (must outlive `calc`). Returns false for unknown configs.
inline bool make_config(int cfg, ceelo::EfficiencyCalculator& calc,
                        ConfigSetup& setup) {
    using namespace ceelo;
    auto add_mat = [&setup](Material m) -> const Material* {
        setup.materials.push_back(std::make_unique<Material>(std::move(m)));
        return setup.materials.back().get();
    };

    switch (cfg) {
        case 1: {
            const Material* nai = add_mat(make_NaI());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
            setup.description = "3\"x3\" NaI bare, point source on-axis 10cm";
            return true;
        }
        case 2: {
            const Material* nai = add_mat(make_NaI());
            const Material* al = add_mat(make_Aluminum());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.add_attenuator(al, 0.1, 0.1, 0.0, 7.62);
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
            setup.description = "3\"x3\" NaI + 1mm Al, point source on-axis 10cm";
            return true;
        }
        case 3: {
            const Material* labr3 = add_mat(make_LaBr3());
            const Material* al = add_mat(make_Aluminum());
            calc.set_detector(DetectorShape::Cylinder, labr3, {2.54, 5.08});
            calc.add_attenuator(al, 0.05, 0.05, 0.0, 5.08);
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
            setup.description = "2\"x2\" LaBr3 + 0.5mm Al, point source on-axis 5cm";
            return true;
        }
        case 5: {
            const Material* czt = add_mat(make_CZT());
            calc.set_detector(DetectorShape::Box, czt, {0.5, 0.5, 0.5});
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
            setup.description = "1x1x0.5cm CZT bare, point source on-axis 5cm";
            return true;
        }
        case 6: {
            const Material* nai = add_mat(make_NaI());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            const double d = 15.0;
            const double theta = 45.0 * M_PI / 180.0;
            calc.set_point_source(
                Eigen::Vector3d(d * std::sin(theta), 0.0, -d * std::cos(theta)));
            setup.description = "3\"x3\" NaI bare, point source 45deg off-axis 15cm";
            return true;
        }
        case 7: {
            const Material* nai = add_mat(make_NaI());
            const Material* al = add_mat(make_Aluminum());
            const Material* pb = add_mat(make_Lead());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.add_attenuator(al, 0.1, 0.1, 0.0, 7.62);
            calc.add_attenuator(pb, 0.2, 0.2, 0.0, 7.62);
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -15.0));
            setup.description = "3\"x3\" NaI + 1mm Al + 2mm Pb, point source on-axis 15cm";
            return true;
        }
        case 8: {
            const Material* nai = add_mat(make_NaI());
            const Material* al = add_mat(make_Aluminum());
            const Material* water = add_mat(make_Water());
            const Material* pe = add_mat(make_Polyethylene());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.add_attenuator(al, 0.05, 0.05, 0.0, 7.62);
            calc.set_marinelli_beaker(
                /*well_inner_radius=*/4.3,
                /*well_depth=*/6.0,
                /*outer_radius=*/7.5,
                /*fill_height=*/4.0,
                /*endcap_to_beaker=*/0.5,
                water, pe, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI + 0.5mm Al, Marinelli beaker (water, 2mm PE)";
            return true;
        }
        case 11: {
            const Material* nai = add_mat(make_NaI());
            const Material* fe = add_mat(make_Iron());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
            calc.add_source_shield(fe, 0.5);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI bare, point source 10cm + 0.5cm Fe shield";
            return true;
        }
        case 12: {
            const Material* nai = add_mat(make_NaI());
            const Material* ss = add_mat(make_StainlessSteel304());
            const Material* cellulose = add_mat(make_Cellulose());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -25.0),
                Eigen::Vector3d(5.0, 7.5, 10.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(cellulose);
            calc.add_source_shield(ss, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 10x15x20cm SS304 box + cellulose, 15cm";
            return true;
        }

        // ------ Out-of-sample stress configs for biasing-policy validation ------
        // No GEANT4 references; used for FOM A/B and Gate-A consistency only.

        case 13: {
            // Thin low-Z shield: scattered share of total ~0; the policy
            // should approach the pure direct/cone limit.
            const Material* nai = add_mat(make_NaI());
            const Material* al = add_mat(make_Aluminum());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
            calc.add_source_shield(al, 0.1);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI bare, point 10cm + 1mm Al shield (stress)";
            return true;
        }
        case 14: {
            // Very thick high-Z shield: total efficiency dominated by
            // shield-scattered photons; weight/variance stress case.
            const Material* nai = add_mat(make_NaI());
            const Material* pb = add_mat(make_Lead());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));
            calc.add_source_shield(pb, 2.0);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI bare, point 10cm + 2cm Pb shield (stress)";
            return true;
        }
        case 15: {
            // Large low-density source volume, no shield: per-position cones,
            // source self-attenuation only.
            const Material* nai = add_mat(make_NaI());
            const Material* water = add_mat(make_Water());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -22.0),
                                        /*radius=*/10.0, /*half_length=*/10.0);
            calc.set_source_material(water);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 20cm-dia x 20cm water cylinder at 12cm (stress)";
            return true;
        }
        case 16: {
            // Close geometry (~1 cm gap): near-degenerate cones, recoil
            // electrons from source water can reach the crystal.
            const Material* nai = add_mat(make_NaI());
            const Material* water = add_mat(make_Water());
            const Material* pe = add_mat(make_Polyethylene());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -1.6),
                                        /*radius=*/2.5, /*half_length=*/0.5);
            calc.set_source_material(water);
            calc.add_source_shield(pe, 0.1);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 5cm-dia x 1cm water puck + 1mm PE, 1cm gap (stress)";
            return true;
        }
        case 17: {
            // Tiny detector + big source: very small detector solid angle,
            // forcing synergy.
            const Material* czt = add_mat(make_CZT());
            const Material* soil = add_mat(make_Soil());
            const Material* ss = add_mat(make_StainlessSteel304());
            calc.set_detector(DetectorShape::Box, czt, {0.5, 0.5, 0.5});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -15.0),
                Eigen::Vector3d(5.0, 5.0, 5.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(soil);
            calc.add_source_shield(ss, 0.1);
            calc.enable_source_electron_transport(true);
            setup.description = "1x1x0.5cm CZT, 10cm soil cube + 1mm SS304, 10cm (stress)";
            return true;
        }
        case 19: {
            // Off-axis shielded point source: cone-axis correctness check
            // (axis = normalize(-src_pos), not -z).
            const Material* nai = add_mat(make_NaI());
            const Material* fe = add_mat(make_Iron());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            const double d = 15.0;
            const double theta = 45.0 * M_PI / 180.0;
            calc.set_point_source(
                Eigen::Vector3d(d * std::sin(theta), 0.0, -d * std::cos(theta)));
            calc.add_source_shield(fe, 0.5);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI bare, point 45deg off-axis 15cm + 0.5cm Fe (stress)";
            return true;
        }
        // ------ High-Z skin-escape validation (cfg-12 geometry, wall = Pb / W) ------
        // Identical to cfg 12 except the box wall material, to validate the
        // skin-escape albedo gate's high-Z (capped) behaviour against GEANT4.
        case 20: {
            const Material* nai = add_mat(make_NaI());
            const Material* pb = add_mat(make_Lead());
            const Material* cellulose = add_mat(make_Cellulose());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -25.0),
                Eigen::Vector3d(5.0, 7.5, 10.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(cellulose);
            calc.add_source_shield(pb, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 10x15x20cm Pb-wall box + cellulose, 15cm";
            return true;
        }
        case 21: {
            const Material* nai = add_mat(make_NaI());
            const Material* w = add_mat(make_Tungsten());
            const Material* cellulose = add_mat(make_Cellulose());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -25.0),
                Eigen::Vector3d(5.0, 7.5, 10.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(cellulose);
            calc.add_source_shield(w, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 10x15x20cm W-wall box + cellulose, 15cm";
            return true;
        }
        case 22: {
            // Intermediate-Z (Sn, Z=50) wall: maps the albedo-cap onset between
            // Fe (cfg 12, Z=26, un-capped) and W/Pb (Z=74/82, capped).
            const Material* nai = add_mat(make_NaI());
            const Material* sn = add_mat(make_Tin());
            const Material* cellulose = add_mat(make_Cellulose());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -25.0),
                Eigen::Vector3d(5.0, 7.5, 10.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(cellulose);
            calc.add_source_shield(sn, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 10x15x20cm Sn-wall box + cellulose, 15cm";
            return true;
        }

        case 23: {
            // Light-Z (Al, Z=13) wall: the empirical gate FLOORS Al to fully
            // analog (eta0(Al)=0.153 < 0.16 floor), ignoring Al's real ~15%
            // backscatter; the principled tails include it. G4 adjudicates the
            // light-Z end where the Fe-anchored gate has no calibration.
            const Material* nai = add_mat(make_NaI());
            const Material* al = add_mat(make_Aluminum());
            const Material* cellulose = add_mat(make_Cellulose());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -25.0),
                Eigen::Vector3d(5.0, 7.5, 10.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(cellulose);
            calc.add_source_shield(al, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 10x15x20cm Al-wall box + cellulose, 15cm";
            return true;
        }

        case 24: {
            // Very-low-Z (polyethylene, Z_eff~5.3) wall: both methods ~analog;
            // confirms the light-Z limit (walk accurate, ~no backscatter to add)
            // matches G4.
            const Material* nai = add_mat(make_NaI());
            const Material* pe = add_mat(make_Polyethylene());
            const Material* cellulose = add_mat(make_Cellulose());
            calc.set_detector(DetectorShape::Cylinder, nai, {3.81, 7.62});
            calc.set_rectangular_source(
                Eigen::Vector3d(0.0, 0.0, -25.0),
                Eigen::Vector3d(5.0, 7.5, 10.0),
                Eigen::Matrix3d::Identity());
            calc.set_source_material(cellulose);
            calc.add_source_shield(pe, 0.2);
            calc.enable_source_electron_transport(true);
            setup.description = "3\"x3\" NaI, 10x15x20cm PE-wall box + cellulose, 15cm";
            return true;
        }

        // Configs 25/26 are a matched pair: the same coaxial HPGe with a sharp
        // and a bulletized front edge.  Efficiencies are compared against
        // GEANT4 individually, but the quantity the pair exists to validate is
        // the *difference* between them -- the geometry change, isolated from
        // the physics residuals both configs share.  Dimensions follow the
        // ANGLE GEM35-70 crystal.  No dead layer: write_gdml() does not export
        // one, so including it would make the two codes model different solids.
        case 25: {
            const Material* ge = add_mat(make_HPGe());
            calc.set_detector(DetectorShape::Cylinder, ge, {2.915, 6.89});
            calc.set_bore_hole(0.495, 5.54);
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
            setup.description = "GEM35-70 HPGe coax, SHARP front edge, point source 5cm";
            return true;
        }

        case 26: {
            const Material* ge = add_mat(make_HPGe());
            calc.set_detector(DetectorShape::Cylinder, ge, {2.915, 6.89});
            calc.set_bullet_radius(0.8);
            calc.set_bore_hole(0.495, 5.54, /*rounded_tip=*/true);
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -5.0));
            setup.description =
                "GEM35-70 HPGe coax, BULLETIZED front edge + rounded bore tip, point source 5cm";
            return true;
        }

        // Configs 27/28: the cfg-26 crystal behind a 0.5 cm Fe spherical shell
        // around a point source, at contact-like (2 cm) and 10 cm.  At 60 keV
        // the shell is 4.4 non-Rayleigh mfp and 0.36 Rayleigh mfp, so ~9% of the
        // peak is Rayleigh-deflected photons that are then absorbed - the term
        // InterSpec's analytic model was missing.  This pair measures how well
        // CeeLo places the scattered-into-peak stream there (u/s split printed
        // by both codes), since its cfg 8/12 stream reads 2-4% low at 59 keV.
        case 27:
        case 28: {
            const Material* ge = add_mat(make_HPGe());
            const Material* fe = add_mat(make_Iron());
            calc.set_detector(DetectorShape::Cylinder, ge, {2.915, 6.89});
            calc.set_bullet_radius(0.8);
            calc.set_bore_hole(0.495, 5.54, /*rounded_tip=*/true);
            const double d_cm = (cfg == 27) ? 2.0 : 10.0;
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -d_cm));
            calc.add_source_shield(fe, 0.5);
            setup.description = std::string("GEM35-70 HPGe coax bulletized, point source + 0.5 cm Fe shell, ")
                                + ((cfg == 27) ? "2 cm" : "10 cm");
            return true;
        }

        default:
            return false;
    }
}

} // namespace bench
