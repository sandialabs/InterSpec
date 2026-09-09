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
#include "test_fep_window.h"
#include "export/Geant4Export.h"
#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
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

/// Whole exported file as a string.
std::string slurp(const std::string& path) {
    std::ifstream f(path);
    return std::string((std::istreambuf_iterator<char>(f)),
                       std::istreambuf_iterator<char>());
}

double attr_of(const std::string& line, const std::string& key) {
    const size_t p = line.find(key + "=\"");
    if (p == std::string::npos) return -1.0;
    return std::atof(line.c_str() + p + key.size() + 2);
}

/// One solid of the exported source chain, in file order (innermost first).
struct SrcSolid {
    std::string name;
    std::string kind;   ///< "sphere" | "tube" | "box"
    double d0 = 0.0;    ///< sphere rmax / tube rmax / box full-x
    double d1 = 0.0;    ///< sphere rmin / tube full-z  / box full-y
    double d2 = 0.0;    ///< box full-z
};

std::vector<SrcSolid> read_src_solids(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::vector<SrcSolid> out;
    while (std::getline(f, line)) {
        const size_t n = line.find("name=\"Src");
        if (n == std::string::npos) continue;
        if (line.find("<sphere") == std::string::npos
            && line.find("<tube") == std::string::npos
            && line.find("<box") == std::string::npos) continue;
        const size_t q0 = line.find('"', n) + 1;
        const size_t q1 = line.find('"', q0);
        SrcSolid s;
        s.name = line.substr(q0, q1 - q0);
        if (line.find("<sphere") != std::string::npos) {
            s.kind = "sphere"; s.d0 = attr_of(line, "rmax"); s.d1 = attr_of(line, "rmin");
        } else if (line.find("<tube") != std::string::npos) {
            s.kind = "tube";   s.d0 = attr_of(line, "rmax"); s.d1 = attr_of(line, "z");
        } else {
            s.kind = "box";    s.d0 = attr_of(line, "x");
            s.d1 = attr_of(line, "y"); s.d2 = attr_of(line, "z");
        }
        out.push_back(s);
    }
    return out;
}

/// Material of each Src* logical volume, and the daughter it carries.
struct SrcVol { std::string lv, mat, daughter_lv; };

std::vector<SrcVol> read_src_volumes(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::vector<SrcVol> out;
    bool in_src = false;
    while (std::getline(f, line)) {
        if (line.find("<volume name=\"Src") != std::string::npos) {
            const size_t q0 = line.find('"') + 1, q1 = line.find('"', q0);
            out.push_back({line.substr(q0, q1 - q0), "", ""});
            in_src = true;
            continue;
        }
        if (!in_src) continue;
        if (line.find("</volume>") != std::string::npos) { in_src = false; continue; }
        if (line.find("<materialref") != std::string::npos) {
            const size_t q0 = line.find("ref=\"") + 5, q1 = line.find('"', q0);
            out.back().mat = line.substr(q0, q1 - q0);
        } else if (line.find("<volumeref") != std::string::npos) {
            const size_t q0 = line.find("ref=\"") + 5, q1 = line.find('"', q0);
            out.back().daughter_lv = line.substr(q0, q1 - q0);
        }
    }
    return out;
}

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
        calc.set_fep_window_keV(kTestFepWindowKeV);
        calc.set_detector(&ge, CylinderDims{R, L});
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
        calc.set_fep_window_keV(kTestFepWindowKeV);
        calc.set_detector(&ge, CylinderDims{R, L});
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
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&ge, CylinderDims{R, L});

    const std::string path = tmp_gdml("plain");
    calc.export_geant4_gdml(path, /*vacuum_world=*/true);

    std::ifstream f(path);
    const std::string txt((std::istreambuf_iterator<char>(f)),
                          std::istreambuf_iterator<char>());
    BOOST_CHECK(txt.find("<tube name=\"CrystalOuterTube\"") != std::string::npos);
    BOOST_CHECK(txt.find("<polycone") == std::string::npos);
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// Source cores.  These close the export<->tracer loop for the concentric source
// stack the same way the polycone tests do for the crystal: the GDML is parsed
// back and checked against what the ray tracer actually traverses.  The GEANT4
// references are generated FROM this export, so a mismatch here would not show
// up as a test failure -- it would show up as a physics disagreement that looks
// like a transport bug.
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cored_sphere_export_matches_traced_segments) {
    // A soil shell [2,3] cm with a 1 cm iron core and a 1 cm lead core inside
    // it, leaving a 0-1 cm cavity, plus an outer Fe shield.
    Material nai = make_NaI(), soil = make_Soil();
    Material fe = make_Iron(), pb = make_Lead();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    calc.set_spherical_source(Eigen::Vector3d(0, 0, -10.0), 3.0,
                              Eigen::Matrix3d::Identity(), 2.0);
    calc.set_source_material(&soil);
    calc.add_source_core(&fe, 1.0);   // fills [1,2]
    calc.add_source_core(&pb, 1.0);   // fills [0,1]... leaves nothing
    calc.add_source_shield(&fe, 0.5);

    const std::string path = tmp_gdml("cored_sphere");
    calc.export_geant4_gdml(path);   // must NOT throw any more

    const std::vector<SrcSolid> solids = read_src_solids(path);
    BOOST_REQUIRE_EQUAL(solids.size(), 4u);   // Pb core, Fe core, soil, Fe shield

    // File order is innermost-first, which is what GDML needs (a volume must be
    // defined before the mother that references it).
    BOOST_CHECK_EQUAL(solids[0].name, "SrcCoreSolid1");
    BOOST_CHECK_EQUAL(solids[1].name, "SrcCoreSolid0");
    BOOST_CHECK_EQUAL(solids[2].name, "SrcMaterialSolid");
    BOOST_CHECK_EQUAL(solids[3].name, "SrcShieldSolid0");

    // Every solid is FULL (a daughter displaces its mother), so only the
    // innermost keeps the 1e-4 cm centre hole G4 navigation wants.
    BOOST_CHECK_CLOSE(solids[0].d1, 1e-4, 1e-6);
    for (size_t i = 1; i < solids.size(); ++i)
        BOOST_CHECK_SMALL(solids[i].d1, 1e-12);

    BOOST_CHECK_CLOSE(solids[0].d0, 1.0, 1e-9);
    BOOST_CHECK_CLOSE(solids[1].d0, 2.0, 1e-9);
    BOOST_CHECK_CLOSE(solids[2].d0, 3.0, 1e-9);
    BOOST_CHECK_CLOSE(solids[3].d0, 3.5, 1e-9);

    // The closure: a ray from the centre outward must cross exactly those
    // boundaries, in that order, with those lengths.
    const SourceGeometry& sg = calc.source_geometry();
    std::vector<SourceGeometry::SourcePathSegment> segs;
    sg.trace_source_segments(Eigen::Vector3d(0, 0, -10.0), Eigen::Vector3d(0, 0, 1),
                             662.0, segs);
    BOOST_REQUIRE_EQUAL(segs.size(), 4u);
    double prev = 0.0;
    for (size_t i = 0; i < 4; ++i) {
        BOOST_CHECK_CLOSE(segs[i].length, solids[i].d0 - prev, 1e-6);
        prev = solids[i].d0;
    }
    BOOST_CHECK_EQUAL(segs[0].material, &pb);
    BOOST_CHECK_EQUAL(segs[1].material, &fe);
    BOOST_CHECK_EQUAL(segs[2].material, &soil);
    BOOST_CHECK_EQUAL(segs[3].material, &fe);

    // Nesting: each volume carries the next one in, and only the outermost is
    // placed in the world.  This is what makes /gps/pos/confine SrcMaterialPV
    // sample the shell alone -- IsSourceConfined() locates the DEEPEST volume.
    const std::vector<SrcVol> vols = read_src_volumes(path);
    BOOST_REQUIRE_EQUAL(vols.size(), 4u);
    BOOST_CHECK_EQUAL(vols[0].daughter_lv, "");              // innermost
    BOOST_CHECK_EQUAL(vols[1].daughter_lv, "SrcCoreLV1");
    BOOST_CHECK_EQUAL(vols[2].daughter_lv, "SrcCoreLV0");    // cores inside source
    BOOST_CHECK_EQUAL(vols[3].daughter_lv, "SrcMaterialLV");

    const std::string txt = slurp(path);
    BOOST_CHECK(txt.find("<subtraction") == std::string::npos);  // no booleans
    BOOST_CHECK(txt.find("SrcShieldPV0") != std::string::npos);  // world holds outermost

    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(partly_filled_core_leaves_a_vacuum_cavity) {
    // Cores that do not reach the centre must leave a real void -- the GDML twin
    // of partly_filled_core_leaves_a_void_but_keeps_the_distances.
    Material nai = make_NaI(), soil = make_Soil(), fe = make_Iron();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    calc.set_spherical_source(Eigen::Vector3d(0, 0, -10.0), 3.0,
                              Eigen::Matrix3d::Identity(), 2.0);
    calc.set_source_material(&soil);
    calc.add_source_core(&fe, 0.5);   // fills [1.5, 2]; [0,1.5] stays empty

    const std::string path = tmp_gdml("partial_core");
    calc.export_geant4_gdml(path);

    const std::vector<SrcSolid> solids = read_src_solids(path);
    BOOST_REQUIRE_EQUAL(solids.size(), 3u);
    BOOST_CHECK_EQUAL(solids[0].name, "SrcVoidSolid");
    BOOST_CHECK_CLOSE(solids[0].d0, 1.5, 1e-9);
    BOOST_CHECK_EQUAL(solids[1].name, "SrcCoreSolid0");
    BOOST_CHECK_CLOSE(solids[1].d0, 2.0, 1e-9);

    // The cavity must be vacuum, not the world material: the MC charges it no
    // attenuation whatever the world is made of.
    const std::vector<SrcVol> vols = read_src_volumes(path);
    BOOST_REQUIRE(!vols.empty());
    BOOST_CHECK_EQUAL(vols[0].lv, "SrcVoidLV");
    BOOST_CHECK_EQUAL(vols[0].mat, "Vacuum");

    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(hollow_source_without_cores_exports_a_vacuum_cavity) {
    // No cores at all: the cavity is still a void daughter rather than a
    // subtraction solid, so the shell is expressed the same way either way.
    Material nai = make_NaI(), soil = make_Soil();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    calc.set_spherical_source(Eigen::Vector3d(0, 0, -10.0), 3.0,
                              Eigen::Matrix3d::Identity(), 2.0);
    calc.set_source_material(&soil);

    const std::string path = tmp_gdml("hollow_nocore");
    calc.export_geant4_gdml(path);

    const std::vector<SrcSolid> solids = read_src_solids(path);
    BOOST_REQUIRE_EQUAL(solids.size(), 2u);
    BOOST_CHECK_EQUAL(solids[0].name, "SrcVoidSolid");
    BOOST_CHECK_CLOSE(solids[0].d0, 2.0, 1e-9);
    BOOST_CHECK_EQUAL(solids[1].name, "SrcMaterialSolid");
    BOOST_CHECK_CLOSE(solids[1].d0, 3.0, 1e-9);
    BOOST_CHECK(slurp(path).find("<subtraction") == std::string::npos);

    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(cores_are_additive_in_the_export) {
    // One 2 cm core and four 0.5 cm cores of the same material describe the same
    // scene; the exported material distribution must agree.  This is the GDML
    // twin of cores_are_additive, and the G-E additivity control in GEANT4.
    Material nai = make_NaI(), soil = make_Soil(), fe = make_Iron();

    auto build = [&](EfficiencyCalculator& c, int n_cores) {
        c.set_fep_window_keV(kTestFepWindowKeV);
        c.set_detector(&nai, CylinderDims{3.81, 7.62});
        c.set_spherical_source(Eigen::Vector3d(0, 0, -10.0), 3.0,
                               Eigen::Matrix3d::Identity(), 2.0);
        c.set_source_material(&soil);
        for (int i = 0; i < n_cores; ++i)
            c.add_source_core(&fe, 2.0 / n_cores);
    };

    EfficiencyCalculator one, many;
    build(one, 1);
    build(many, 4);

    const std::string p1 = tmp_gdml("core_add1"), p4 = tmp_gdml("core_add4");
    one.export_geant4_gdml(p1);
    many.export_geant4_gdml(p4);

    const std::vector<SrcSolid> s1 = read_src_solids(p1);
    const std::vector<SrcSolid> s4 = read_src_solids(p4);
    BOOST_REQUIRE_EQUAL(s1.size(), 2u);   // one core + source
    BOOST_REQUIRE_EQUAL(s4.size(), 5u);   // four cores + source

    // Same outer boundary for the iron, same source shell, and every
    // intermediate boundary is an iron/iron interface that changes nothing.
    BOOST_CHECK_CLOSE(s1[0].d0, 2.0, 1e-9);
    BOOST_CHECK_CLOSE(s4[3].d0, 2.0, 1e-9);
    BOOST_CHECK_CLOSE(s1[1].d0, s4[4].d0, 1e-9);

    // And the tracer agrees: a subdivided core merges back into one run.
    std::vector<SourceGeometry::SourcePathSegment> g1, g4;
    one.source_geometry().trace_source_segments(
        Eigen::Vector3d(0, 0, -10.0), Eigen::Vector3d(0, 0, 1), 662.0, g1);
    many.source_geometry().trace_source_segments(
        Eigen::Vector3d(0, 0, -10.0), Eigen::Vector3d(0, 0, 1), 662.0, g4);
    BOOST_REQUIRE_EQUAL(g1.size(), g4.size());
    for (size_t i = 0; i < g1.size(); ++i)
        BOOST_CHECK_CLOSE(g1[i].length, g4[i].length, 1e-9);

    std::remove(p1.c_str());
    std::remove(p4.c_str());
}

BOOST_AUTO_TEST_CASE(nested_cylinder_exports_a_closed_cavity_not_a_pipe) {
    // The export twin of nested_cylinder_is_not_a_pipe.  Before nesting, the
    // cylinder source solid was written as a plain through-bore <tube> that
    // ignored cyl_inner_half_length() entirely -- a nested stack came out a pipe.
    Material nai = make_NaI(), soil = make_Soil(), fe = make_Iron();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    // R=3, half-length 3; cavity r=2, half-length 2 => CLOSED, not a through-bore.
    calc.set_cylindrical_source(Eigen::Vector3d(0, 0, -10.0), 3.0, 3.0,
                                Eigen::Matrix3d::Identity(), 2.0, 2.0);
    calc.set_source_material(&soil);
    calc.add_source_core(&fe, 2.0, 2.0);

    const std::string path = tmp_gdml("nested_cyl");
    calc.export_geant4_gdml(path);

    const std::vector<SrcSolid> solids = read_src_solids(path);
    BOOST_REQUIRE_EQUAL(solids.size(), 2u);
    BOOST_CHECK_EQUAL(solids[0].kind, "tube");
    BOOST_CHECK_CLOSE(solids[0].d0, 2.0, 1e-9);   // core radius
    BOOST_CHECK_CLOSE(solids[0].d1, 4.0, 1e-9);   // FULL z = 2*2, NOT 2*3
    BOOST_CHECK_CLOSE(solids[1].d0, 3.0, 1e-9);
    BOOST_CHECK_CLOSE(solids[1].d1, 6.0, 1e-9);

    // The tracer sees the same thing, and it is what distinguishes the two: a
    // ray leaving the axis at |z| = 2.5 (past the cavity end) crosses the full
    // 3 cm radius of soil.  Modelled as a through-bore it would cross only
    // 3 - 2 = 1 cm, because the bore would run the whole length.
    const SourceGeometry& sg = calc.source_geometry();
    std::vector<SourceGeometry::SourcePathSegment> segs;
    sg.trace_source_segments(Eigen::Vector3d(0, 0, -10.0 + 2.5),
                             Eigen::Vector3d(1, 0, 0), 662.0, segs);
    double soil_len = 0.0;
    for (const auto& sgm : segs) if (sgm.material == &soil) soil_len += sgm.length;
    BOOST_CHECK_CLOSE(soil_len, 3.0, 1e-6);

    EfficiencyCalculator pipe;
    pipe.set_fep_window_keV(kTestFepWindowKeV);
    pipe.set_detector(&nai, CylinderDims{3.81, 7.62});
    pipe.set_cylindrical_source(Eigen::Vector3d(0, 0, -10.0), 3.0, 3.0,
                                Eigen::Matrix3d::Identity(), 2.0);  // through-bore
    pipe.set_source_material(&soil);
    std::vector<SourceGeometry::SourcePathSegment> psegs;
    pipe.source_geometry().trace_source_segments(
        Eigen::Vector3d(0, 0, -10.0 + 2.5), Eigen::Vector3d(1, 0, 0), 662.0, psegs);
    double pipe_len = 0.0;
    for (const auto& sgm : psegs) if (sgm.material == &soil) pipe_len += sgm.length;
    BOOST_CHECK_CLOSE(pipe_len, 1.0, 1e-6);

    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(cored_box_exports_nested_boxes) {
    Material nai = make_NaI(), soil = make_Soil(), fe = make_Iron();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    calc.set_rectangular_source(Eigen::Vector3d(0, 0, -20.0),
                                Eigen::Vector3d(5.0, 6.0, 7.0),
                                Eigen::Matrix3d::Identity(),
                                Eigen::Vector3d(3.0, 4.0, 5.0));
    calc.set_source_material(&soil);
    calc.add_source_core(&fe, 3.0, 4.0, 5.0);   // fills the cavity exactly

    const std::string path = tmp_gdml("cored_box");
    calc.export_geant4_gdml(path);

    const std::vector<SrcSolid> solids = read_src_solids(path);
    BOOST_REQUIRE_EQUAL(solids.size(), 2u);
    BOOST_CHECK_EQUAL(solids[0].kind, "box");
    BOOST_CHECK_CLOSE(solids[0].d0, 6.0,  1e-9);   // 2 * 3
    BOOST_CHECK_CLOSE(solids[0].d1, 8.0,  1e-9);
    BOOST_CHECK_CLOSE(solids[0].d2, 10.0, 1e-9);
    BOOST_CHECK_CLOSE(solids[1].d0, 10.0, 1e-9);   // 2 * 5
    BOOST_CHECK(slurp(path).find("<subtraction") == std::string::npos);

    std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(source_without_material_is_vacuum_not_shield_material) {
    // An extended source with no set_source_material() is a real volume that the
    // MC charges NOTHING for.  Under nesting it must still be emitted, as vacuum,
    // or its shield -- now a full solid rather than a hollow shell -- would
    // silently fill the source region with shield material.  Regression guard for
    // exactly that: the bug is invisible in the solid dimensions and shows up only
    // as a physics disagreement.
    Material nai = make_NaI(), pb = make_Lead();

    EfficiencyCalculator calc;
    calc.set_fep_window_keV(kTestFepWindowKeV);
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0, 0, -15.0), 3.0, 3.0);
    calc.add_source_shield(&pb, 0.3, 0.1);   // note: no set_source_material()

    const std::string path = tmp_gdml("nomat");
    calc.export_geant4_gdml(path);

    const std::vector<SrcSolid> solids = read_src_solids(path);
    BOOST_REQUIRE_EQUAL(solids.size(), 2u);
    BOOST_CHECK_EQUAL(solids[0].name, "SrcMaterialSolid");
    BOOST_CHECK_CLOSE(solids[0].d0, 3.0, 1e-9);   // the source volume is present
    BOOST_CHECK_CLOSE(solids[1].d0, 3.3, 1e-9);

    const std::vector<SrcVol> vols = read_src_volumes(path);
    BOOST_REQUIRE_EQUAL(vols.size(), 2u);
    BOOST_CHECK_EQUAL(vols[0].lv,  "SrcMaterialLV");
    BOOST_CHECK_EQUAL(vols[0].mat, "Vacuum");     // NOT Pb
    BOOST_CHECK_EQUAL(vols[1].mat, "Pb");
    BOOST_CHECK_EQUAL(vols[1].daughter_lv, "SrcMaterialLV");

    std::remove(path.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
