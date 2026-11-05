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

/// @file env_scatter_export.cpp
/// @brief Emit GEANT4 inputs (GDML + GPS macro) matching a CeeLo
///        EnvironmentGeometry, for validating the environment-scatter kernel.
///
/// The geometry written here is deliberately the SAME abstraction CeeLo scores:
/// an air world, optional material slabs (ground / room surfaces), and a
/// **scoring sphere made of air** at the detector position. The sphere is
/// physically transparent — it is air inside air — so it is a pure scoring
/// boundary and does not perturb the transport, exactly matching CeeLo's
/// assumption that the scoring body is a void in the field.
///
/// The sphere is named `active_crystal` because the harness resolves its scoring
/// volume by that name; with `--sphere-flux R` the harness then tallies photon
/// TRACK LENGTH inside it, binned in (E, mu) — phi = L/V, which is precisely the
/// observable `NextEventEstimator` produces. Scoring the same quantity on both
/// sides is what makes the comparison apples-to-apples; a disk-crossing scorer
/// would be a DIFFERENT observable and any mismatch would be indistinguishable
/// from a physics discrepancy.
///
/// No detector is placed. That matches the CeeLo kernel (which reports
/// incident flux, not detected counts) and isolates the environment physics.
///
/// Usage:
///   env_scatter_export --mode field --distance 500 --h-src 100 --h-det 100 \
///                      --energy 662 --sphere-radius 5 --events 200000000 \
///                      --prefix v1_field_662

#include "materials/Material.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

struct Slab {
    std::string name;
    double cx, cy, cz;          // centre, cm
    double hx, hy, hz;          // half-dimensions, cm
    std::string material;
};

} // namespace

int main(int argc, char** argv)
{
    std::string mode = "field", material = "soil", prefix = "envscatter";
    double distance = 500.0, h_src = 100.0, h_det = 100.0;
    double energy = 662.0, sphere_radius = 4.0;
    double world_half = 0.0;          // cm; 0 => auto
    double slab_thickness = 60.0;     // cm
    // room extents (cm): height, half-width, back, front
    double room_h = 250.0, room_w = 150.0, room_b = 100.0, room_f = 100.0;
    long long events = 100000000;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto nxt = [&]() { return (i + 1 < argc) ? argv[++i] : ""; };
        if (a == "--mode") mode = nxt();
        else if (a == "--material") material = nxt();
        else if (a == "--distance") distance = std::atof(nxt());
        else if (a == "--h-src") h_src = std::atof(nxt());
        else if (a == "--h-det") h_det = std::atof(nxt());
        else if (a == "--energy") energy = std::atof(nxt());
        else if (a == "--sphere-radius") sphere_radius = std::atof(nxt());
        else if (a == "--world-half") world_half = std::atof(nxt());
        else if (a == "--slab-thickness") slab_thickness = std::atof(nxt());
        else if (a == "--room") {
            room_h = std::atof(nxt()); room_w = std::atof(nxt());
            room_b = std::atof(nxt()); room_f = std::atof(nxt());
        }
        else if (a == "--events") events = std::atoll(nxt());
        else if (a == "--prefix") prefix = nxt();
        else { std::cerr << "unknown option " << a << "\n"; return 1; }
    }

    const Material air = make_Air();
    const Material soil = make_Soil();
    const Material concrete = make_Concrete();
    const bool is_room = (mode == "room");
    const bool is_air  = (mode == "air");
    const Material& slab_mat = (material == "concrete") ? concrete : soil;
    // Absolute height has no physical meaning in homogeneous air.  Canonical
    // air campaign nodes use 0/0, so emit that identity even if the caller left
    // the field-oriented height defaults unchanged.
    const double fit_h_src = is_air ? 0.0 : h_src;
    const double fit_h_det = is_air ? 0.0 : h_det;

    // Detector at the origin, source on the -z axis. The ground plane is normal
    // to +y, placed so the detector sits h_det above it. A height DIFFERENCE is
    // expressed by tilting nothing -- instead the source is raised, which is
    // equivalent for a flat ground and keeps every solid axis-aligned (G4
    // navigation is happier, and it matches CeeLo's field() plane, which solves
    // for exactly these two perpendicular distances).
    // CeeLo's field() puts the source at |r| = distance and TILTS the ground
    // plane to meet both heights. Here the ground stays axis-aligned and the
    // source is raised instead -- so src_z must be shortened to keep the
    // source-detector SEPARATION equal to `distance`, not just its z component.
    // Getting this wrong offsets the whole comparison by (d_true/d)^2: with
    // h_src=30, h_det=100, d=500 the naive src_z=-500 gives |r|=504.9 cm and a
    // spurious -1.9% in every flux.
    const double src_y = fit_h_src - fit_h_det;
    if (std::abs(src_y) >= distance) {
        std::cerr << "|h_src - h_det| must be < distance\n";
        return 1;
    }
    const double src_z = -std::sqrt(distance * distance - src_y * src_y);

    std::vector<Slab> slabs;
    // The world must contain the contributing halo.
    //
    // GEANT4 kills any track that leaves the world, while CeeLo's field() is an
    // infinite half-space, so a world that is too small charges CeeLo for flux
    // the reference structurally cannot produce -- a one-directional bias that
    // is invisible in a total-vs-total comparison because it looks exactly like
    // a physics excess. A 15 m half-size at the campaign's 5 m standoff (for
    // example max(3*distance, 1500 cm)) loses 2-3% of the scattered flux.
    // Measured by emulating this boundary inside CeeLo
    // (EnvScatterConfig::world_half_size_cm),
    // soil field at 5 m:
    //
    //     world half   122 keV   662 keV   2614 keV
    //         15 m      -3.1%     -2.2%      -2.1%
    //         30 m      -1.4%     -0.2%      +0.1%
    //         60 m      -0.2%     -0.2%      -0.8%
    //
    // The binding scale is geometric rather than attenuation: contributions fall
    // off about as r^-3 per unit radius, so the tail is thin but long, and the
    // air mean free path (54 m at 122 keV, 107 m at 662) is not what sets it.
    // 60 m puts the truncation below the statistical precision of a feasible
    // reference run, so that is the floor.
    //
    // A larger world costs GEANT4 tracking time. When that matters, keep the
    // world small and set CeeLo's `world_half_size_cm` to the SAME value so the
    // comparison is geometry-matched -- but then say so, because such a
    // reference is not the physical infinite-halfspace case.
    double wh = (world_half > 0.0) ? world_half
                                   : std::max(12.0 * distance, 6000.0);
    if (!is_air) {
        // Ground: top face at y = -h_det. Deep enough for saturated albedo
        // (>= ~5 mean free paths; 60 cm of soil at 662 keV is ~7).
        slabs.push_back({"ground", 0.0, -h_det - 0.5 * slab_thickness, 0.0,
                         wh, 0.5 * slab_thickness, wh,
                         "MatSlab"});
    }
    if (is_room) {
        const double y_ceiling = -h_det + room_h;
        slabs.push_back({"ceiling", 0.0, y_ceiling + 0.5 * slab_thickness, 0.0,
                         room_w + slab_thickness, 0.5 * slab_thickness,
                         (distance + room_b + room_f) * 0.5 + slab_thickness,
                         "MatSlab"});
        const double zc = -0.5 * distance;             // room centre along z
        const double zh = 0.5 * distance + 0.5 * (room_b + room_f);
        slabs.push_back({"wall_px", room_w + 0.5 * slab_thickness,
                         -h_det + 0.5 * room_h, zc,
                         0.5 * slab_thickness, 0.5 * room_h, zh, "MatSlab"});
        slabs.push_back({"wall_nx", -room_w - 0.5 * slab_thickness,
                         -h_det + 0.5 * room_h, zc,
                         0.5 * slab_thickness, 0.5 * room_h, zh, "MatSlab"});
        slabs.push_back({"wall_back", 0.0, -h_det + 0.5 * room_h,
                         room_b + 0.5 * slab_thickness,
                         room_w, 0.5 * room_h, 0.5 * slab_thickness, "MatSlab"});
        slabs.push_back({"wall_front", 0.0, -h_det + 0.5 * room_h,
                         -distance - room_f - 0.5 * slab_thickness,
                         room_w, 0.5 * room_h, 0.5 * slab_thickness, "MatSlab"});
        wh = std::max(wh, (distance + room_b + room_f) + 4.0 * slab_thickness);
    }

    // ---- GDML ----
    const std::string gdml = prefix + ".gdml";
    std::ofstream os(gdml);
    if (!os) { std::cerr << "cannot write " << gdml << "\n"; return 1; }
    os << std::setprecision(17);
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
       << "<gdml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
       << " xsi:noNamespaceSchemaLocation=\"http://service-spi.web.cern.ch/"
          "service-spi/app/releases/GDML/schema/gdml.xsd\">\n"
       << "  <!-- CEEFIT-GEOMETRY-v1 geometry_tier=" << mode
       << " material=" << (is_air ? "none" : material)
       << " energy_keV=" << energy
       << " distance_cm=" << distance
       << " source_height_cm=" << fit_h_src
       << " detector_height_cm=" << fit_h_det
       << " scoring_radius_cm=" << sphere_radius
       << " world_half_extent_cm=" << wh << " -->\n"
       << "  <define/>\n  <materials>\n";
    // Compose from NIST ELEMENTAL materials (G4_O, G4_Si, ...) by mass
    // fraction. GDML allows a material to be built from other materials, so this
    // avoids hand-rolling isotopic compositions while still reproducing CeeLo's
    // make_Soil()/make_Concrete() exactly: same mass fractions, same density.
    // G4_AIR and G4_CONCRETE are taken straight from NIST.
    {
        static const char* kSym[31] = {
            "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
            "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
            "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"};
        auto emit = [&](const char* name, const Material& m) {
            os << "    <material name=\"" << name << "\" state=\"solid\">\n"
               << "      <D value=\"" << m.density() << "\" unit=\"g/cm3\"/>\n";
            for (const MaterialComponent& c : m.composition())
                os << "      <fraction n=\"" << c.mass_fraction
                   << "\" ref=\"G4_" << kSym[c.Z] << "\"/>\n";
            os << "    </material>\n";
        };
        emit("MatSoil", soil);
        emit("MatConcrete", concrete);
    }
    os << "  </materials>\n  <solids>\n"
       << "    <box name=\"world_s\" x=\"" << 2 * wh << "\" y=\"" << 2 * wh
       << "\" z=\"" << 2 * wh << "\" lunit=\"cm\"/>\n"
       << "    <sphere name=\"score_s\" rmin=\"0\" rmax=\"" << sphere_radius
       << "\" startphi=\"0\" deltaphi=\"360\" starttheta=\"0\""
       << " deltatheta=\"180\" aunit=\"deg\" lunit=\"cm\"/>\n";
    for (const Slab& s : slabs)
        os << "    <box name=\"" << s.name << "_s\" x=\"" << 2 * s.hx
           << "\" y=\"" << 2 * s.hy << "\" z=\"" << 2 * s.hz
           << "\" lunit=\"cm\"/>\n";
    os << "  </solids>\n  <structure>\n";
    // Scoring sphere: AIR, so it is physically transparent and is a pure
    // scoring boundary. Named active_crystal because the harness looks that up.
    os << "    <volume name=\"active_crystal\">\n"
       << "      <materialref ref=\"G4_AIR\"/>\n"
       << "      <solidref ref=\"score_s\"/>\n    </volume>\n";
    for (const Slab& s : slabs)
        os << "    <volume name=\"" << s.name << "\">\n"
           << "      <materialref ref=\""
           << ((material == "concrete") ? "MatConcrete" : "MatSoil")
           << "\"/>\n      <solidref ref=\"" << s.name << "_s\"/>\n"
           << "    </volume>\n";
    os << "    <volume name=\"world\">\n"
       << "      <materialref ref=\"G4_AIR\"/>\n"
       << "      <solidref ref=\"world_s\"/>\n"
       << "      <physvol><volumeref ref=\"active_crystal\"/>"
       << "<position name=\"score_p\" x=\"0\" y=\"0\" z=\"0\" unit=\"cm\"/>"
       << "</physvol>\n";
    for (const Slab& s : slabs)
        os << "      <physvol><volumeref ref=\"" << s.name << "\"/>"
           << "<position name=\"" << s.name << "_p\" x=\"" << s.cx << "\" y=\""
           << s.cy << "\" z=\"" << s.cz << "\" unit=\"cm\"/></physvol>\n";
    os << "    </volume>\n  </structure>\n"
       << "  <setup name=\"Default\" version=\"1.0\">"
       << "<world ref=\"world\"/></setup>\n</gdml>\n";
    os.close();

    // ---- GPS macro: ISOTROPIC point source. Cone biasing is invalid here --
    // environment scatter arrives from all directions, which is the entire
    // reason the kernel exists.
    const std::string mac = prefix + ".mac";
    std::ofstream ms(mac);
    ms << std::setprecision(17)
       << "/run/initialize\n"
       << "/process/em/verbose 0\n"
       << "/run/verbose 0\n"
       << "/event/verbose 0\n"
       << "/tracking/verbose 0\n"
       << "/gps/particle gamma\n"
       << "/gps/energy " << energy << " keV\n"
       << "/gps/pos/type Point\n"
       << "/gps/pos/centre 0 " << src_y << " " << src_z << " cm\n"
       << "/gps/ang/type iso\n"
       << "/run/beamOn " << events << "\n";
    ms.close();

    std::cout << "Wrote " << gdml << " and " << mac << "\n"
              << "  mode=" << mode << (is_air ? "" : (" material=" + material))
              << " d=" << distance << " cm  h_src=" << fit_h_src
              << " h_det=" << fit_h_det << " cm\n"
              << "  source at (0, " << src_y << ", " << src_z << ") cm; "
              << "scoring sphere R=" << sphere_radius << " cm at origin\n"
              << "  world half-size " << wh << " cm; slabs: " << slabs.size()
              << "\n  Run:  ceelo_g4val " << gdml << " " << mac
              << " out.csv --sphere-flux " << sphere_radius << " --seed SEED\n";
    return 0;
}
