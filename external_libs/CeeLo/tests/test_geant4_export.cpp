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

/// @file test_geant4_export.cpp
/// @brief The exported GDML must describe the solid the tracer traces.
///
/// The GEANT4 references are generated from this export, so a discrepancy
/// between the exported solid and the traced solid does not show up as a test
/// failure -- it shows up as a physics disagreement that looks like a transport
/// bug.  These tests close that loop directly: for a crystal, integrate the
/// scoring chord the ray tracer returns over the front face, and compare with
/// the volume of the polycone that was written out.
///
/// The bulletized/bored profile is where this matters.  A merge that truncated
/// the fillet arc once cost 3.6% of the crystal volume for a bore deep enough
/// to reach the front fillet -- silently, and only for configurations no
/// benchmark used.

#define BOOST_TEST_MODULE Geant4ExportTests
#include <boost/test/unit_test.hpp>

#include "efficiency/EfficiencyCalculator.h"
#include "export/Geant4Export.h"
#include "geometry/Geometry.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

struct ZPlane { double rmin, rmax, z; };

/// Pull the crystal polycone's z-planes back out of an exported GDML file.
/// Returns empty when the crystal was written as a plain <tube> instead.
std::vector<ZPlane> read_crystal_polycone(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::vector<ZPlane> planes;
    bool in_crystal = false;

    auto attr = [](const std::string& s, const std::string& key) {
        const size_t p = s.find(key + "=\"");
        if (p == std::string::npos) return 0.0;
        return std::atof(s.c_str() + p + key.size() + 2);
    };

    while (std::getline(f, line)) {
        if (line.find("<polycone name=\"CrystalSolid\"") != std::string::npos) {
            in_crystal = true;
            continue;
        }
        if (!in_crystal) continue;
        if (line.find("</polycone>") != std::string::npos) break;
        if (line.find("<zplane") != std::string::npos) {
            planes.push_back({attr(line, "rmin"), attr(line, "rmax"), attr(line, "z")});
        }
    }
    return planes;
}

/// Volume of the conical frusta the z-planes describe -- exact, since GDML
/// interpolates linearly between planes.
double polycone_volume(const std::vector<ZPlane>& p) {
    double v = 0.0;
    for (size_t i = 0; i + 1 < p.size(); ++i) {
        const double dz = p[i + 1].z - p[i].z;
        if (dz == 0.0) continue;
        const double outer = p[i].rmax * p[i].rmax + p[i].rmax * p[i + 1].rmax
                           + p[i + 1].rmax * p[i + 1].rmax;
        const double inner = p[i].rmin * p[i].rmin + p[i].rmin * p[i + 1].rmin
                           + p[i + 1].rmin * p[i + 1].rmin;
        v += M_PI / 3.0 * dz * (outer - inner);
    }
    return v;
}

/// Volume the transport code actually sees: integrate its own scoring chord
/// over the front face.  Axis-parallel rays, so this is exact up to the radial
/// quadrature.
double traced_volume(const Geometry& g, double R) {
    const int N = 40000;
    const Eigen::Vector3d dir(0.0, 0.0, 1.0);
    std::vector<PathSegment> segs;
    double v = 0.0;
    for (int i = 0; i < N; ++i) {
        const double rho = R * (i + 0.5) / N;
        g.trace_ray(Eigen::Vector3d(rho, 0.0, -5.0), dir, segs);
        double chord = 0.0;
        for (const auto& s : segs) if (s.is_scoring) chord += s.length();
        v += chord * 2.0 * M_PI * rho * (R / N);
    }
    return v;
}

std::string tmp_gdml(const char* tag) {
    std::ostringstream oss;
    oss << "ceelo_export_test_" << tag << ".gdml";
    return oss.str();
}

// GEM35-70.
constexpr double R = 2.915, L = 6.89, RB = 0.8, RBORE = 0.495, DEPTH = 5.54;

} // namespace


BOOST_AUTO_TEST_SUITE(Geant4Export)

BOOST_AUTO_TEST_CASE(exported_polycone_matches_traced_solid) {
    // Every bore/fillet combination the API admits, including the deep bore
    // whose closed end reaches back into the front fillet -- the case that
    // once truncated the arc -- and a bore wider than the fillet ring radius.
    struct Case { const char* tag; double r_b, bore_r, depth; bool tip; };
    const Case cases[] = {
        {"sharp_flat",   0.0, RBORE, DEPTH, false},
        {"bullet_flat",  RB,  RBORE, DEPTH, false},
        {"bullet_round", RB,  RBORE, DEPTH, true },
        {"bullet_nobore",RB,  0.0,   0.0,   false},
        {"deep_flat",    RB,  RBORE, 6.50,  false},
        {"deep_round",   RB,  RBORE, 6.50,  true },
        {"wide_bore",    RB,  2.00,  DEPTH, false},
        {"shallow",      RB,  RBORE, 0.60,  true },
    };

    Material ge = make_HPGe();
    for (const auto& c : cases) {
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &ge, {R, L});
        if (c.r_b > 0.0) calc.set_bullet_radius(c.r_b);
        if (c.bore_r > 0.0) calc.set_bore_hole(c.bore_r, c.depth, c.tip);

        const std::string path = tmp_gdml(c.tag);
        calc.export_geant4_gdml(path, /*vacuum_world=*/true);

        const auto planes = read_crystal_polycone(path);
        BOOST_REQUIRE_MESSAGE(!planes.empty(),
                              "no CrystalSolid polycone exported for " << c.tag);

        // A polycone is only well formed if z is monotone and no plane has a
        // bore wider than the crystal there.
        for (size_t i = 0; i < planes.size(); ++i) {
            BOOST_REQUIRE_MESSAGE(planes[i].rmin <= planes[i].rmax + 1e-12,
                                  c.tag << ": rmin > rmax at plane " << i);
            if (i) BOOST_REQUIRE_MESSAGE(planes[i].z >= planes[i - 1].z - 1e-12,
                                         c.tag << ": z not monotone at plane " << i);
        }

        const double v_export = polycone_volume(planes);
        const double v_traced = traced_volume(calc.geometry(), R);
        BOOST_CHECK_MESSAGE(std::abs(v_export - v_traced) / v_traced < 2e-4,
                            c.tag << ": exported " << v_export
                                  << " cm^3 vs traced " << v_traced << " cm^3");
        std::remove(path.c_str());
    }
}

BOOST_AUTO_TEST_CASE(exported_fillet_volume_is_right_not_just_close) {
    // Normalising to the whole crystal hides the fillet: the removed corner is
    // only 1.3% of the volume, so a badly wrong arc still lands within a
    // fraction of a percent of the total.  Compare the *difference* against
    // the Pappus closed form instead.
    Material ge = make_HPGe();

    auto exported_volume = [&](double r_b) {
        EfficiencyCalculator calc;
        calc.set_detector(DetectorShape::Cylinder, &ge, {R, L});
        if (r_b > 0.0) calc.set_bullet_radius(r_b);
        // No bore: isolate the fillet.
        const std::string path = tmp_gdml(r_b > 0.0 ? "fv_bullet" : "fv_sharp");
        calc.export_geant4_gdml(path, /*vacuum_world=*/true);
        const auto planes = read_crystal_polycone(path);
        double v;
        if (planes.empty()) {
            v = M_PI * R * R * L;   // sharp crystals still export as a <tube>
        } else {
            v = polycone_volume(planes);
        }
        std::remove(path.c_str());
        return v;
    };

    const double rho_c = R - RB;
    const double v_removed_exact = 2.0 * M_PI * rho_c * RB * RB * (1.0 - M_PI / 4.0)
                                 + M_PI * RB * RB * RB / 3.0;
    const double v_removed_export = exported_volume(0.0) - exported_volume(RB);

    // 0.1 % *of the fillet*, which is ~1.3e-3 % of the crystal -- tight enough
    // that a wrong arc radius or a truncated sweep cannot pass.
    BOOST_CHECK_CLOSE(v_removed_export, v_removed_exact, 0.1);
}

BOOST_AUTO_TEST_CASE(plain_cylinder_still_exports_as_a_tube) {
    // Pre-existing configs must keep byte-identical geometry: only crystals
    // with a bore and/or a fillet become polycones.
    Material ge = make_HPGe();
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &ge, {R, L});

    const std::string path = tmp_gdml("plain");
    calc.export_geant4_gdml(path, /*vacuum_world=*/true);

    std::ifstream f(path);
    const std::string txt((std::istreambuf_iterator<char>(f)),
                          std::istreambuf_iterator<char>());
    BOOST_CHECK(txt.find("<tube name=\"CrystalOuterTube\"") != std::string::npos);
    BOOST_CHECK(txt.find("<polycone") == std::string::npos);
    std::remove(path.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
