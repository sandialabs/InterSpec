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

#include "export/Geant4Export.h"
#include "geometry/Geometry.h"
#include "geometry/SourceGeometry.h"
#include "materials/Material.h"

#include <Eigen/Core>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <set>
#include <vector>
#include <cmath>
#include <iomanip>

namespace ceelo {

namespace {

// ---------------------------------------------------------------------------
// Element data table: Z → (symbol, IUPAC name, atomic weight g/mol)
// ---------------------------------------------------------------------------
struct ElementInfo {
    uint8_t     Z;
    const char* symbol;
    const char* name;
    double      atomic_weight_g_per_mol;
};

static const ElementInfo kElementData[] = {
    {  1, "H",  "Hydrogen",    1.008    },
    {  2, "He", "Helium",      4.003    },
    {  3, "Li", "Lithium",     6.941    },
    {  4, "Be", "Beryllium",   9.012    },
    {  5, "B",  "Boron",      10.811    },
    {  6, "C",  "Carbon",     12.011    },
    {  7, "N",  "Nitrogen",   14.007    },
    {  8, "O",  "Oxygen",     15.999    },
    {  9, "F",  "Fluorine",   18.998    },
    { 10, "Ne", "Neon",       20.180    },
    { 11, "Na", "Sodium",     22.990    },
    { 12, "Mg", "Magnesium",  24.305    },
    { 13, "Al", "Aluminum",   26.982    },
    { 14, "Si", "Silicon",    28.086    },
    { 15, "P",  "Phosphorus", 30.974    },
    { 16, "S",  "Sulfur",     32.065    },
    { 17, "Cl", "Chlorine",   35.450    },
    { 18, "Ar", "Argon",      39.948    },
    { 19, "K",  "Potassium",  39.098    },
    { 20, "Ca", "Calcium",    40.078    },
    { 21, "Sc", "Scandium",   44.956    },
    { 22, "Ti", "Titanium",   47.867    },
    { 23, "V",  "Vanadium",   50.942    },
    { 24, "Cr", "Chromium",   51.996    },
    { 25, "Mn", "Manganese",  54.938    },
    { 26, "Fe", "Iron",       55.845    },
    { 27, "Co", "Cobalt",     58.933    },
    { 28, "Ni", "Nickel",     58.693    },
    { 29, "Cu", "Copper",     63.546    },
    { 30, "Zn", "Zinc",       65.380    },
    { 31, "Ga", "Gallium",    69.723    },
    { 32, "Ge", "Germanium",  72.630    },
    { 33, "As", "Arsenic",    74.922    },
    { 34, "Se", "Selenium",   78.971    },
    { 35, "Br", "Bromine",    79.904    },
    { 36, "Kr", "Krypton",    83.798    },
    { 37, "Rb", "Rubidium",   85.468    },
    { 38, "Sr", "Strontium",  87.620    },
    { 39, "Y",  "Yttrium",    88.906    },
    { 40, "Zr", "Zirconium",  91.224    },
    { 41, "Nb", "Niobium",    92.906    },
    { 42, "Mo", "Molybdenum", 95.960    },
    { 43, "Tc", "Technetium", 98.000    },
    { 44, "Ru", "Ruthenium", 101.070    },
    { 45, "Rh", "Rhodium",   102.906    },
    { 46, "Pd", "Palladium", 106.420    },
    { 47, "Ag", "Silver",    107.868    },
    { 48, "Cd", "Cadmium",   112.414    },
    { 49, "In", "Indium",    114.818    },
    { 50, "Sn", "Tin",       118.710    },
    { 51, "Sb", "Antimony",  121.760    },
    { 52, "Te", "Tellurium", 127.600    },
    { 53, "I",  "Iodine",    126.904    },
    { 54, "Xe", "Xenon",     131.293    },
    { 55, "Cs", "Cesium",    132.905    },
    { 56, "Ba", "Barium",    137.327    },
    { 57, "La", "Lanthanum", 138.905    },
    { 58, "Ce", "Cerium",    140.116    },
    { 59, "Pr", "Praseodymium",140.908  },
    { 60, "Nd", "Neodymium", 144.242    },
    { 62, "Sm", "Samarium",  150.360    },
    { 63, "Eu", "Europium",  151.964    },
    { 64, "Gd", "Gadolinium",157.250    },
    { 65, "Tb", "Terbium",   158.925    },
    { 66, "Dy", "Dysprosium",162.500    },
    { 67, "Ho", "Holmium",   164.930    },
    { 68, "Er", "Erbium",    167.259    },
    { 69, "Tm", "Thulium",   168.934    },
    { 70, "Yb", "Ytterbium", 173.045    },
    { 71, "Lu", "Lutetium",  174.967    },
    { 72, "Hf", "Hafnium",   178.490    },
    { 73, "Ta", "Tantalum",  180.948    },
    { 74, "W",  "Tungsten",  183.840    },
    { 75, "Re", "Rhenium",   186.207    },
    { 76, "Os", "Osmium",    190.230    },
    { 77, "Ir", "Iridium",   192.217    },
    { 78, "Pt", "Platinum",  195.084    },
    { 79, "Au", "Gold",      196.967    },
    { 80, "Hg", "Mercury",   200.592    },
    { 81, "Tl", "Thallium",  204.383    },
    { 82, "Pb", "Lead",      207.200    },
    { 83, "Bi", "Bismuth",   208.980    },
    { 84, "Po", "Polonium",  209.000    },
    { 85, "At", "Astatine",  210.000    },
    { 86, "Rn", "Radon",     222.000    },
    { 87, "Fr", "Francium",  223.000    },
    { 88, "Ra", "Radium",    226.000    },
    { 89, "Ac", "Actinium",  227.000    },
    { 90, "Th", "Thorium",   232.038    },
    { 91, "Pa", "Protactinium", 231.036 },
    { 92, "U",  "Uranium",   238.029    },
};

const ElementInfo* find_element(uint8_t Z) {
    for (const auto& el : kElementData) {
        if (el.Z == Z) return &el;
    }
    return nullptr;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Sanitize a material name for use as an XML id (replace spaces with underscores).
std::string sanitize_name(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out += (c == ' ' || c == '(' || c == ')') ? '_' : c;
    }
    return out;
}

// Collect all unique Z values from a list of materials.
std::set<uint8_t> collect_z_values(const std::vector<const Material*>& mats) {
    std::set<uint8_t> zs;
    // Always include Air components: N(7), O(8), Ar(18)
    zs.insert(7);
    zs.insert(8);
    zs.insert(18);
    for (const Material* m : mats) {
        if (!m) continue;
        for (const auto& comp : m->composition()) {
            zs.insert(comp.Z);
        }
    }
    return zs;
}

// Format a double with fixed precision, suppressing trailing zeros.
std::string fmt(double v, int prec = 6) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(prec) << v;
    return oss.str();
}

// ---------------------------------------------------------------------------
// GDML generation
// ---------------------------------------------------------------------------

void write_gdml_elements(std::ostream& out, const std::set<uint8_t>& zs) {
    for (uint8_t Z : zs) {
        const ElementInfo* el = find_element(Z);
        if (!el) continue;
        out << "    <element name=\"" << el->symbol << "\" formula=\"" << el->symbol
            << "\" Z=\"" << static_cast<int>(Z) << "\">\n"
            << "      <atom type=\"A\" unit=\"g/mol\" value=\"" << fmt(el->atomic_weight_g_per_mol, 4) << "\"/>\n"
            << "    </element>\n";
    }
}

void write_gdml_material(std::ostream& out, const Material* mat) {
    std::string safe = sanitize_name(mat->name());
    out << "    <material name=\"" << safe << "\" state=\"solid\">\n"
        << "      <D type=\"density\" unit=\"g/cm3\" value=\"" << fmt(mat->density(), 5) << "\"/>\n";
    for (const auto& comp : mat->composition()) {
        const ElementInfo* el = find_element(comp.Z);
        if (!el) continue;
        out << "      <fraction n=\"" << fmt(comp.mass_fraction, 6) << "\" ref=\"" << el->symbol << "\"/>\n";
    }
    out << "    </material>\n";
}

void write_gdml_air(std::ostream& out) {
    // Standard dry air composition by mass: N2(75.52%), O2(23.20%), Ar(1.28%).
    // Fractions sum to 1.0 to avoid GEANT4 material normalization warnings.
    out << "    <material name=\"Air\" state=\"gas\">\n"
        << "      <D type=\"density\" unit=\"g/cm3\" value=\"0.001205\"/>\n"
        << "      <fraction n=\"0.7552\" ref=\"N\"/>\n"
        << "      <fraction n=\"0.2320\" ref=\"O\"/>\n"
        << "      <fraction n=\"0.0128\" ref=\"Ar\"/>\n"
        << "    </material>\n";
}

} // anonymous namespace


// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

void write_gdml(const Geometry& geom, const std::string& filename,
                const SourceGeometry* source_geom,
                bool vacuum_world) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("write_gdml: cannot open file: " + filename);
    }

    // Collect all materials used in the geometry.
    std::vector<const Material*> all_mats;
    all_mats.push_back(geom.detector_material());
    for (const auto& att : geom.attenuators()) {
        all_mats.push_back(att.material);
    }

    // Source geometry materials.
    bool has_source_geom = source_geom && source_geom->has_source_effects()
                           && source_geom->is_configured();
    if (has_source_geom) {
        if (source_geom->source_material()) {
            all_mats.push_back(source_geom->source_material());
        }
        for (const auto& layer : source_geom->shields()) {
            all_mats.push_back(layer.material);
        }
    }

    // Collect element Z values.
    auto zs = collect_z_values(all_mats);

    // Detector dimensions.
    const double L  = geom.detector_length();
    const bool is_cyl = (geom.shape() == DetectorShape::Cylinder);
    const double R  = is_cyl ? geom.detector_radius() : 0.0;
    const double hx = is_cyl ? 0.0 : geom.detector_half_x();
    const double hy = is_cyl ? 0.0 : geom.detector_half_y();

    // World box half-dimensions (generous envelope around everything).
    const double det_outer_r = geom.outer_bounding_radius();  // detector only (for overlap avoidance)
    double outer_r = det_outer_r;
    auto [z_min_geom, z_max_geom] = geom.outer_z_extent();

    // Enlarge world to contain source geometry if present.
    if (has_source_geom) {
        double src_extent = source_geom->outermost_extent_radius();
        using SShape = SourceGeometry::Shape;
        Eigen::Vector3d src_center = Eigen::Vector3d::Zero();
        if (source_geom->shape() == SShape::Point)
            src_center = source_geom->point_position();
        else if (source_geom->shape() == SShape::Cylindrical)
            src_center = source_geom->cyl_center();
        else if (source_geom->shape() == SShape::Rectangular)
            src_center = source_geom->rect_center();
        // Marinelli: src_center stays at origin (axially symmetric around detector)

        if (source_geom->shape() == SShape::Marinelli) {
            // Marinelli: outer radius + shield thicknesses, z from z_bot to z_we
            if (src_extent > outer_r) outer_r = src_extent;
            if (source_geom->marinelli_z_bot() < z_min_geom)
                z_min_geom = source_geom->marinelli_z_bot();
            if (source_geom->marinelli_z_we() > z_max_geom)
                z_max_geom = source_geom->marinelli_z_we();
        } else {
            double src_r = std::sqrt(src_center.x()*src_center.x()
                                     + src_center.y()*src_center.y()) + src_extent;
            if (src_r > outer_r) outer_r = src_r;
            double src_z_min = src_center.z() - src_extent;
            double src_z_max = src_center.z() + src_extent;
            if (src_z_min < z_min_geom) z_min_geom = src_z_min;
            if (src_z_max > z_max_geom) z_max_geom = src_z_max;
        }
    }

    double world_hxy = outer_r * 3.0 + 5.0;
    double world_hz  = (std::max(std::abs(z_min_geom), std::abs(z_max_geom)) + 5.0) * 3.0;

    // Bore hole data.  The bore is part of the crystal polycone's profile, so
    // no separate bore solid or placement is needed.
    const bool has_bore = geom.has_bore_hole();
    double R_bore = 0.0, D_bore = 0.0;
    bool round_bore_tip = false;
    if (has_bore && is_cyl) {
        R_bore = geom.bore_hole()->radius;
        D_bore = geom.bore_hole()->depth;
        round_bore_tip = geom.bore_hole()->rounded_tip;
    }

    // Bulletized (rounded) front outer edge.  The fillet's ring circle has
    // radius `bullet_rho_c` and sits at local z = -L/2 + r_b.
    const double bullet_r   = is_cyl ? geom.bullet_radius() : 0.0;
    const bool has_bullet   = (bullet_r > 0.0);
    const double bullet_rho_c = R - bullet_r;
    // Material the fillet removes, by Pappus: 2 pi rho_c r_b^2 (1 - pi/4)
    // + pi r_b^3 / 3.  Reported in a comment for the validation logs.
    constexpr double kPi = 3.14159265358979323846;
    const double bullet_removed_vol =
        has_bullet ? (2.0 * kPi * bullet_rho_c * bullet_r * bullet_r
                          * (1.0 - kPi / 4.0)
                      + kPi * bullet_r * bullet_r * bullet_r / 3.0)
                   : 0.0;

    // Dead layer data.
    double dl_front = 0.0, dl_side = 0.0, dl_back = 0.0;
    if (geom.has_dead_layer()) {
        dl_front = geom.dead_layer()->front;
        dl_side  = geom.dead_layer()->side;
        dl_back  = geom.dead_layer()->back;
    }

    // Sanitized crystal material name.
    std::string crystal_mat_name = sanitize_name(geom.detector_material()->name());

    // Pre-compute GDML geometry parameters for each attenuator layer, mirroring the
    // ray-trace accumulation logic so GDML and MC ray-trace use identical geometry
    // (no volume overlaps, correct nesting of front endcaps).
    struct AttGDMLParams {
        double outer_r;       ///< Outer radius of cup
        double inner_r;       ///< Inner radius (crystal or previous layer outer)
        double outer_z;       ///< Total z-extent of outer solid
        double inner_z;       ///< z-extent of inner cavity solid
        double center_z;      ///< World z of solid center for physvol placement
        double inner_offset;  ///< z offset of inner solid center relative to outer solid center
        bool use_subtraction; ///< true when front endcap is present
    };

    std::vector<AttGDMLParams> att_gdml;
    {
        double acc_r     = is_cyl ? R : std::sqrt(hx*hx + hy*hy);
        if (geom.has_dead_layer()) acc_r += dl_side;
        double acc_z_min = 0.0;  // crystal front face
        double acc_z_max = L;    // crystal back face

        for (const auto& att : geom.attenuators()) {
            AttGDMLParams p;
            p.inner_r = acc_r;
            p.outer_r = acc_r + att.side_thickness;

            double front_z_min = acc_z_min - att.front_thickness;
            double shell_z_min = std::min(front_z_min, att.z_start);
            double shell_z_max = std::max(acc_z_max, att.z_end);

            p.outer_z   = shell_z_max - shell_z_min;
            p.inner_z   = acc_z_max  - acc_z_min;
            p.center_z  = (shell_z_min + shell_z_max) / 2.0;

            double inner_center_world = (acc_z_min + acc_z_max) / 2.0;
            p.inner_offset    = inner_center_world - p.center_z;
            p.use_subtraction = (att.front_thickness > 0.0 && p.inner_z > 0.0
                                 && !att.side_only);

            att_gdml.push_back(p);

            // Update accumulated inner geometry for the next layer
            acc_r     = p.outer_r;
            acc_z_min = shell_z_min;
            acc_z_max = shell_z_max;
        }
    }

    // -----------------------------------------------------------------------
    // XML header
    // -----------------------------------------------------------------------
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<gdml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
        << "      xsi:noNamespaceSchemaLocation=\""
        << "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd\">\n\n";

    // -----------------------------------------------------------------------
    // <define>: positions used in the structure section
    // -----------------------------------------------------------------------
    out << "  <define>\n"
        << "    <!-- Crystal center in world frame (front face at z=0, center at z=L/2) -->\n"
        << "    <position name=\"crystal_center\" x=\"0\" y=\"0\" z=\"" << fmt(L/2.0) << "\" unit=\"cm\"/>\n";

    // Attenuator inner-cavity offset positions (used in subtraction solids for cup-shaped endcaps).
    for (size_t i = 0; i < att_gdml.size(); ++i) {
        if (att_gdml[i].use_subtraction) {
            out << "    <position name=\"att" << i << "_inner_offset\""
                << " x=\"0\" y=\"0\" z=\"" << fmt(att_gdml[i].inner_offset) << "\""
                << " unit=\"cm\"/>\n";
        }
    }

    out << "    <rotation name=\"identity\" x=\"0\" y=\"0\" z=\"0\" unit=\"deg\"/>\n"
        << "  </define>\n\n";

    // -----------------------------------------------------------------------
    // <materials>
    // -----------------------------------------------------------------------
    out << "  <materials>\n";
    write_gdml_elements(out, zs);
    out << "\n";
    write_gdml_air(out);
    if (vacuum_world) {
        out << "    <material name=\"Vacuum\" state=\"gas\">\n"
            << "      <D type=\"density\" unit=\"g/cm3\" value=\"1.0e-25\"/>\n"
            << "      <fraction n=\"1.0\" ref=\"N\"/>\n"
            << "    </material>\n";
    }
    out << "\n";
    // Write all unique detector/attenuator materials.
    std::set<const Material*> written_mats;
    for (const Material* m : all_mats) {
        if (!m) continue;
        if (written_mats.insert(m).second) {
            write_gdml_material(out, m);
        }
    }
    out << "  </materials>\n\n";

    // -----------------------------------------------------------------------
    // <solids>
    // -----------------------------------------------------------------------
    out << "  <solids>\n";

    // World box.
    out << "    <!-- World envelope -->\n"
        << "    <box name=\"WorldBox\""
        << " x=\"" << fmt(2.0*world_hxy) << "\""
        << " y=\"" << fmt(2.0*world_hxy) << "\""
        << " z=\"" << fmt(2.0*world_hz) << "\""
        << " lunit=\"cm\"/>\n\n";

    if (is_cyl && (has_bore || has_bullet)) {
        // A coaxial and/or bulletized crystal is a surface of revolution, so it
        // exports as ONE native G4Polycone rather than a stack of boolean
        // solids.  That matters for more than tidiness: an earlier
        // subtraction-based version of this export produced ~0.3% stuck-track
        // warnings under cone-biased runs (G4 navigation struggles with tracks
        // skimming a subtracted surface), which both floods the log and risks
        // biasing the very comparison the export exists to support.  The
        // polycone has none.
        //
        // Profile in crystal-local z (front face at -L/2, back at +L/2):
        //   rmax  follows the front fillet from rho_c up to R over the first
        //         r_b of depth, then stays at R.
        //   rmin  is 0 until the bore starts, then the bore radius (via the
        //         tip hemisphere first, if the bore is round-tipped).
        // Circular arcs become chord segments; kArcSeg keeps the sagitta well
        // under a micron, i.e. ~1e-4 % of the crystal volume.
        constexpr int kArcSeg = 64;
        struct ZPlane { double z, rmin, rmax; };
        std::vector<ZPlane> profile;

        const double z_front = -L / 2.0;
        const double z_back  =  L / 2.0;

        // --- rmax: front fillet, then the straight side wall ---
        if (has_bullet) {
            // Arc centre at (rho_c, z_front + r_b); walk the quarter circle
            // from the front face outward to the side wall.
            for (int i = 0; i <= kArcSeg; ++i) {
                const double ang = (kPi / 2.0) * i / kArcSeg;   // 0 -> 90 deg
                const double z = z_front + bullet_r * (1.0 - std::cos(ang));
                const double r = bullet_rho_c + bullet_r * std::sin(ang);
                profile.push_back({z, 0.0, r});
            }
        } else {
            profile.push_back({z_front, 0.0, R});
        }
        profile.push_back({z_back, 0.0, R});

        // --- rmin: the bore, cut into the profile from the back ---
        if (has_bore) {
            const double z_bore_start = z_back - D_bore;   // apex of the bore
            std::vector<ZPlane> bore;
            if (round_bore_tip) {
                const double z_tip = z_bore_start + R_bore;  // hemisphere centre
                for (int i = 0; i <= kArcSeg; ++i) {
                    const double ang = (kPi / 2.0) * i / kArcSeg;
                    bore.push_back({z_tip - R_bore * std::cos(ang),
                                    R_bore * std::sin(ang), 0.0});
                }
            } else {
                bore.push_back({z_bore_start, 0.0, 0.0});      // flat bottom:
                bore.push_back({z_bore_start, R_bore, 0.0});   // rmin steps up
            }
            bore.push_back({z_back, R_bore, 0.0});

            // Merge onto the union of both profiles' z values, evaluating each
            // profile at the other's planes.  Evaluating rather than
            // truncating is what keeps a deep bore -- one whose closed end
            // reaches back into the front fillet -- from cutting the arc short
            // and interpolating rmax straight to the back face.
            auto rmax_at = [&](double z) {
                if (!has_bullet) return R;
                if (z >= z_front + bullet_r) return R;
                const double dz = (z_front + bullet_r) - z;
                return bullet_rho_c
                     + std::sqrt(std::max(0.0, bullet_r * bullet_r - dz * dz));
            };
            auto rmin_at = [&](double z) {
                if (z <= z_bore_start) return 0.0;
                const double z_tip = z_bore_start + R_bore;
                if (!round_bore_tip || z >= z_tip) return R_bore;
                const double dz = z_tip - z;
                return std::sqrt(std::max(0.0, R_bore * R_bore - dz * dz));
            };

            std::vector<ZPlane> merged;
            merged.reserve(profile.size() + bore.size());
            size_t oi = 0, bj = 0;
            while (oi < profile.size() || bj < bore.size()) {
                // On a tie take the outer plane first, so a flat bore's rmin
                // step stays the last thing to happen at that z.
                const bool take_outer = (bj >= bore.size())
                    || (oi < profile.size() && profile[oi].z <= bore[bj].z);
                if (take_outer) {
                    merged.push_back({profile[oi].z, rmin_at(profile[oi].z),
                                      profile[oi].rmax});
                    ++oi;
                } else {
                    merged.push_back({bore[bj].z, bore[bj].rmin,
                                      rmax_at(bore[bj].z)});
                    ++bj;
                }
            }
            // Both profiles end at the back face, so drop the duplicate.
            profile.clear();
            for (const auto& p : merged) {
                if (!profile.empty()) {
                    const ZPlane& q = profile.back();
                    if (q.z == p.z && q.rmin == p.rmin && q.rmax == p.rmax) continue;
                }
                profile.push_back(p);
            }
        }

        out << "    <!-- Detector crystal (polycone: "
            << (has_bullet ? "bulletized front edge" : "sharp front edge")
            << (has_bore ? (round_bore_tip ? ", round-tipped bore" : ", flat-bottomed bore")
                         : "")
            << ") -->\n";
        if (has_bullet) {
            out << "    <!-- Bulletizing radius " << fmt(bullet_r,4)
                << " cm about ring radius " << fmt(bullet_rho_c,4) << " cm;"
                << " removes " << fmt(bullet_removed_vol,4) << " cm^3 of crystal -->\n";
        }
        out << "    <polycone name=\"CrystalSolid\""
            << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\">\n";
        for (const auto& p : profile) {
            out << "      <zplane rmin=\"" << fmt(p.rmin) << "\""
                << " rmax=\"" << fmt(p.rmax) << "\""
                << " z=\"" << fmt(p.z) << "\"/>\n";
        }
        out << "    </polycone>\n";
    } else if (is_cyl) {
        // Plain cylinder: a single tube, exactly as before.
        out << "    <!-- Detector crystal outer solid -->\n"
            << "    <tube name=\"CrystalOuterTube\""
            << " rmin=\"0\" rmax=\"" << fmt(R) << "\""
            << " z=\"" << fmt(L) << "\""
            << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";
    } else {
        // Box detector.
        out << "    <!-- Detector crystal box solid -->\n"
            << "    <box name=\"CrystalSolid\""
            << " x=\"" << fmt(2.0*hx) << "\""
            << " y=\"" << fmt(2.0*hy) << "\""
            << " z=\"" << fmt(L) << "\""
            << " lunit=\"cm\"/>\n";
    }

    // Attenuator solids: cup shape with solid front endcap + hollow cylindrical sides.
    // Each layer uses the pre-computed att_gdml params (accumulated state matches ray-trace).
    for (size_t i = 0; i < att_gdml.size(); ++i) {
        const auto& att = geom.attenuators()[i];
        const auto& p   = att_gdml[i];

        out << "\n    <!-- Attenuator " << i << " ("
            << (att.material ? att.material->name() : "unknown")
            << ", front=" << fmt(att.front_thickness,3) << " cm"
            << ", side=" << fmt(att.side_thickness,3) << " cm) -->\n";

        if (p.use_subtraction) {
            // Cup: outer full cylinder minus shorter inner full cylinder.
            // Leaves a solid front disk whose thickness = (outer_z - inner_z).
            out << "    <tube name=\"AttOuterSolid" << i << "\""
                << " rmin=\"0\" rmax=\"" << fmt(p.outer_r) << "\""
                << " z=\"" << fmt(p.outer_z) << "\""
                << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n"
                << "    <tube name=\"AttInnerSolid" << i << "\""
                << " rmin=\"0\" rmax=\"" << fmt(p.inner_r) << "\""
                << " z=\"" << fmt(p.inner_z) << "\""
                << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n"
                << "    <subtraction name=\"AttSolid" << i << "\">\n"
                << "      <first ref=\"AttOuterSolid" << i << "\"/>\n"
                << "      <second ref=\"AttInnerSolid" << i << "\"/>\n"
                << "      <positionref ref=\"att" << i << "_inner_offset\"/>\n"
                << "    </subtraction>\n";
        } else {
            // No front endcap (front_thickness == 0): simple hollow annular tube.
            out << "    <tube name=\"AttSolid" << i << "\""
                << " rmin=\"" << fmt(p.inner_r) << "\""
                << " rmax=\"" << fmt(p.outer_r) << "\""
                << " z=\"" << fmt(p.outer_z) << "\""
                << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";
        }
    }

    // Source shielding solids.
    if (has_source_geom) {
        const auto& shields = source_geom->shields();
        using SShape = SourceGeometry::Shape;

        if (source_geom->shape() == SShape::Point) {
            // Spherical shells centered at the point source position.
            // Build nested spheres: innermost shield layer first.
            // Use a small rmin (1e-4 cm) for the innermost sphere to avoid
            // G4 navigation issues when GPS generates a particle at the exact
            // center of a solid sphere (rmin=0). The thin vacuum core has
            // negligible effect on physics.
            double r_inner = 0.0;
            for (size_t i = 0; i < shields.size(); ++i) {
                double r_outer = r_inner + shields[i].scalar_thickness();
                double gdml_rmin = (r_inner < 1e-6) ? 1e-4 : r_inner;
                out << "\n    <!-- Source shield " << i << " ("
                    << (shields[i].material ? shields[i].material->name() : "unknown")
                    << ", t=" << fmt(shields[i].scalar_thickness(),4) << " cm) -->\n"
                    << "    <sphere name=\"SrcShieldSolid" << i << "\""
                    << " rmin=\"" << fmt(gdml_rmin) << "\""
                    << " rmax=\"" << fmt(r_outer) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\""
                    << " starttheta=\"0\" deltatheta=\"3.1415927\""
                    << " aunit=\"rad\" lunit=\"cm\"/>\n";
                r_inner = r_outer;
            }

            // Source material: point sources have no volume, skip.
        } else if (source_geom->shape() == SShape::Cylindrical) {
            double inner_r = source_geom->cyl_radius();
            double inner_half_z = source_geom->cyl_half_length();

            // Source material fill (if present): tube matching source volume.
            // rmin = inner bore radius (0 = solid; a hollow/annular tube has a
            // non-attenuating void core).
            if (source_geom->source_material()) {
                out << "\n    <!-- Source material fill -->\n"
                    << "    <tube name=\"SrcMaterialSolid\""
                    << " rmin=\"" << fmt(source_geom->cyl_inner_radius()) << "\""
                    << " rmax=\"" << fmt(inner_r) << "\""
                    << " z=\"" << fmt(2.0 * inner_half_z) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";
            }

            // Cylindrical shielding shells.
            for (size_t i = 0; i < shields.size(); ++i) {
                double t_r = shields[i].tx;  // radial (tx == ty)
                double t_e = shields[i].tz;  // end caps
                double outer_r = inner_r + t_r;
                double outer_half_z = inner_half_z + t_e;

                // Inflate the subtracted solid by a tiny epsilon on both axes so
                // the shell's inner surface never coincides with the surface it
                // wraps (source material or previous shell). Coincident G4
                // boolean surfaces stall navigation and inflate the electron-
                // entry diagnostic; the 1e-4 cm gap thins the wall ~1 micron
                // (negligible) and, on a zero-thickness axis, lies outside the
                // outer solid so the open-face shape is unchanged.
                constexpr double kZeroDimEps = 1e-4;  // cm
                double sub_r = inner_r + kZeroDimEps;
                double sub_half_z = inner_half_z + kZeroDimEps;

                // Cup-shaped subtraction: outer full cylinder minus inner full cylinder.
                out << "\n    <!-- Source shield " << i << " ("
                    << (shields[i].material ? shields[i].material->name() : "unknown")
                    << ", t_radial=" << fmt(t_r,4) << " cm, t_end=" << fmt(t_e,4) << " cm) -->\n"
                    << "    <tube name=\"SrcShieldOuterSolid" << i << "\""
                    << " rmin=\"0\" rmax=\"" << fmt(outer_r) << "\""
                    << " z=\"" << fmt(2.0 * outer_half_z) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n"
                    << "    <tube name=\"SrcShieldInnerSolid" << i << "\""
                    << " rmin=\"0\" rmax=\"" << fmt(sub_r) << "\""
                    << " z=\"" << fmt(2.0 * sub_half_z) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n"
                    << "    <subtraction name=\"SrcShieldSolid" << i << "\">\n"
                    << "      <first ref=\"SrcShieldOuterSolid" << i << "\"/>\n"
                    << "      <second ref=\"SrcShieldInnerSolid" << i << "\"/>\n"
                    << "    </subtraction>\n";

                inner_r = outer_r;
                inner_half_z = outer_half_z;
            }
        } else if (source_geom->shape() == SShape::Rectangular) {
            const auto& hd = source_geom->rect_half_dims();
            double inner_hx = hd.x();
            double inner_hy = hd.y();
            double inner_hz = hd.z();

            // Source material fill (if present): box matching source volume.
            // A hollow box shell becomes a subtraction solid (outer box minus
            // inner void box); the void is strictly interior, so no
            // coincident-surface epsilon is needed. GPS `/gps/pos/confine
            // SrcMaterialPV` then rejects positions in the void.
            if (source_geom->source_material()) {
                const Eigen::Vector3d& ihd = source_geom->rect_inner_half_dims();
                if (ihd.minCoeff() > 1e-10) {
                    out << "\n    <!-- Source material fill (hollow box shell) -->\n"
                        << "    <box name=\"SrcMaterialOuterSolid\""
                        << " x=\"" << fmt(2.0 * inner_hx) << "\""
                        << " y=\"" << fmt(2.0 * inner_hy) << "\""
                        << " z=\"" << fmt(2.0 * inner_hz) << "\""
                        << " lunit=\"cm\"/>\n"
                        << "    <box name=\"SrcMaterialVoidSolid\""
                        << " x=\"" << fmt(2.0 * ihd.x()) << "\""
                        << " y=\"" << fmt(2.0 * ihd.y()) << "\""
                        << " z=\"" << fmt(2.0 * ihd.z()) << "\""
                        << " lunit=\"cm\"/>\n"
                        << "    <subtraction name=\"SrcMaterialSolid\">\n"
                        << "      <first ref=\"SrcMaterialOuterSolid\"/>\n"
                        << "      <second ref=\"SrcMaterialVoidSolid\"/>\n"
                        << "    </subtraction>\n";
                } else {
                    out << "\n    <!-- Source material fill -->\n"
                        << "    <box name=\"SrcMaterialSolid\""
                        << " x=\"" << fmt(2.0 * inner_hx) << "\""
                        << " y=\"" << fmt(2.0 * inner_hy) << "\""
                        << " z=\"" << fmt(2.0 * inner_hz) << "\""
                        << " lunit=\"cm\"/>\n";
                }
            }

            // Rectangular shielding shells: box subtraction solids.
            for (size_t i = 0; i < shields.size(); ++i) {
                double outer_hx = inner_hx + shields[i].tx;
                double outer_hy = inner_hy + shields[i].ty;
                double outer_hz = inner_hz + shields[i].tz;

                // Inflate the subtracted solid by a tiny epsilon on ALL axes so
                // the shell's inner face never coincides with the face it wraps
                // (the source-material box, or the previous shell). Coincident
                // G4 boolean surfaces stall navigation ("track stuck") and, for
                // extended sources with a shield, spuriously inflate the
                // crystal electron-entry diagnostic. The 1e-4 cm gap thins the
                // wall by ~1 micron (negligible) and, on zero-thickness axes,
                // lies outside the outer solid so the open-face shape is kept.
                constexpr double kZeroDimEps = 1e-4;  // cm
                double sub_hx = inner_hx + kZeroDimEps;
                double sub_hy = inner_hy + kZeroDimEps;
                double sub_hz = inner_hz + kZeroDimEps;

                out << "\n    <!-- Source shield " << i << " ("
                    << (shields[i].material ? shields[i].material->name() : "unknown")
                    << ", t_x=" << fmt(shields[i].tx,4) << " cm, t_y=" << fmt(shields[i].ty,4)
                    << " cm, t_z=" << fmt(shields[i].tz,4) << " cm) -->\n"
                    << "    <box name=\"SrcShieldOuterSolid" << i << "\""
                    << " x=\"" << fmt(2.0 * outer_hx) << "\""
                    << " y=\"" << fmt(2.0 * outer_hy) << "\""
                    << " z=\"" << fmt(2.0 * outer_hz) << "\""
                    << " lunit=\"cm\"/>\n"
                    << "    <box name=\"SrcShieldInnerSolid" << i << "\""
                    << " x=\"" << fmt(2.0 * sub_hx) << "\""
                    << " y=\"" << fmt(2.0 * sub_hy) << "\""
                    << " z=\"" << fmt(2.0 * sub_hz) << "\""
                    << " lunit=\"cm\"/>\n"
                    << "    <subtraction name=\"SrcShieldSolid" << i << "\">\n"
                    << "      <first ref=\"SrcShieldOuterSolid" << i << "\"/>\n"
                    << "      <second ref=\"SrcShieldInnerSolid" << i << "\"/>\n"
                    << "    </subtraction>\n";

                inner_hx = outer_hx;
                inner_hy = outer_hy;
                inner_hz = outer_hz;
            }
        } else if (source_geom->shape() == SShape::Sphere) {
            double inner_r = source_geom->sphere_radius();
            constexpr double kZeroDimEps = 1e-4;  // cm

            // Source material: ball or hollow shell. rmin = inner void radius
            // (clamped to >= 1e-4 cm to avoid G4 GPS issues for a particle at the
            // exact center of a solid sphere).
            if (source_geom->source_material()) {
                double rmin = std::max(source_geom->sphere_inner_radius(), 1e-4);
                out << "\n    <!-- Source material fill (sphere) -->\n"
                    << "    <sphere name=\"SrcMaterialSolid\""
                    << " rmin=\"" << fmt(rmin) << "\""
                    << " rmax=\"" << fmt(inner_r) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\""
                    << " starttheta=\"0\" deltatheta=\"3.1415927\""
                    << " aunit=\"rad\" lunit=\"cm\"/>\n";
            }

            // Spherical shielding shells. Each is a G4Sphere shell [rmin, rmax];
            // rmin is inflated by a tiny epsilon so the shell's inner face does
            // not coincide with the surface it wraps (source material or previous
            // shell) — coincident G4 boolean surfaces stall navigation. The
            // ~1 micron gap is negligible.
            for (size_t i = 0; i < shields.size(); ++i) {
                double outer_r = inner_r + shields[i].scalar_thickness();
                out << "\n    <!-- Source shield " << i << " ("
                    << (shields[i].material ? shields[i].material->name() : "unknown")
                    << ", t=" << fmt(shields[i].scalar_thickness(),4) << " cm) -->\n"
                    << "    <sphere name=\"SrcShieldSolid" << i << "\""
                    << " rmin=\"" << fmt(inner_r + kZeroDimEps) << "\""
                    << " rmax=\"" << fmt(outer_r) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\""
                    << " starttheta=\"0\" deltatheta=\"3.1415927\""
                    << " aunit=\"rad\" lunit=\"cm\"/>\n";
                inner_r = outer_r;
            }
        } else if (source_geom->shape() == SShape::Marinelli) {
            double well_r = source_geom->marinelli_well_inner_radius();
            double outer_r_src = source_geom->marinelli_outer_radius();
            double z_bk = source_geom->marinelli_z_bk();
            double z_we = source_geom->marinelli_z_we();
            double z_bot = source_geom->marinelli_z_bot();
            double fill_h = z_bk - z_bot;        // disk height
            double ring_d = z_we - z_bk;          // ring depth
            double total_h = z_we - z_bot;         // total z-extent

            // Source material: L-shaped solid = outer cylinder minus well cylinder.
            if (source_geom->source_material()) {
                out << "\n    <!-- Marinelli source material (L-shaped subtraction) -->\n"
                    << "    <tube name=\"SrcOuterCyl\""
                    << " rmin=\"0\" rmax=\"" << fmt(outer_r_src) << "\""
                    << " z=\"" << fmt(total_h) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n"
                    << "    <tube name=\"SrcWellCyl\""
                    << " rmin=\"0\" rmax=\"" << fmt(well_r) << "\""
                    << " z=\"" << fmt(ring_d) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n"
                    << "    <subtraction name=\"SrcMaterialSolid\">\n"
                    << "      <first ref=\"SrcOuterCyl\"/>\n"
                    << "      <second ref=\"SrcWellCyl\"/>\n"
                    << "      <position x=\"0\" y=\"0\" z=\"" << fmt(fill_h / 2.0) << "\" unit=\"cm\"/>\n"
                    << "    </subtraction>\n";
            }

            // Beaker wall: decompose into individual tube pieces.
            // Open top at the well opening (z_bk face — toward detector).
            for (size_t i = 0; i < shields.size(); ++i) {
                double t = shields[i].scalar_thickness();
                std::string mat_name = shields[i].material
                    ? sanitize_name(shields[i].material->name()) : "Air";

                out << "\n    <!-- Marinelli beaker wall " << i << " ("
                    << mat_name << ", t=" << fmt(t,4) << " cm) -->\n";

                // Outer cylindrical wall
                out << "    <tube name=\"SrcShieldOuter" << i << "\""
                    << " rmin=\"" << fmt(outer_r_src) << "\""
                    << " rmax=\"" << fmt(outer_r_src + t) << "\""
                    << " z=\"" << fmt(total_h) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";

                // Bottom disk
                out << "    <tube name=\"SrcShieldBottom" << i << "\""
                    << " rmin=\"0\" rmax=\"" << fmt(outer_r_src + t) << "\""
                    << " z=\"" << fmt(t) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";

                // Well inner wall
                out << "    <tube name=\"SrcShieldWellWall" << i << "\""
                    << " rmin=\"" << fmt(well_r - t) << "\""
                    << " rmax=\"" << fmt(well_r) << "\""
                    << " z=\"" << fmt(ring_d) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";

                // Well bottom annulus (separates well cavity from sample below).
                // rmin = detector outer bounding radius to avoid overlapping the crystal.
                double wb_rmin = det_outer_r + 0.001;  // tiny gap to avoid coincident surfaces
                if (wb_rmin > well_r - t) wb_rmin = well_r - t;  // clamp if detector nearly fills well
                out << "    <tube name=\"SrcShieldWellBottom" << i << "\""
                    << " rmin=\"" << fmt(wb_rmin) << "\" rmax=\"" << fmt(well_r - t) << "\""
                    << " z=\"" << fmt(t) << "\""
                    << " startphi=\"0\" deltaphi=\"6.2831853\" aunit=\"rad\" lunit=\"cm\"/>\n";
            }
        }
    }

    out << "  </solids>\n\n";

    // -----------------------------------------------------------------------
    // <structure>
    // -----------------------------------------------------------------------
    out << "  <structure>\n";

    // Dead-layer note (not explicitly modeled as separate GDML volume in this
    // first implementation; the crystal material provides attenuation).
    if (geom.has_dead_layer()) {
        out << "    <!-- NOTE: dead layer (front=" << fmt(dl_front,4) << " cm,"
            << " side=" << fmt(dl_side,4) << " cm,"
            << " back=" << fmt(dl_back,4) << " cm) is NOT explicitly modeled as a\n"
            << "         separate GDML volume.  The full crystal is scored.\n"
            << "         For accurate dead-layer comparison, score only the inner region. -->\n\n";
    }

    // Crystal logical volume (named "active_crystal" for SteppingAction scoring).
    const std::string crystal_solid_ref =
        (is_cyl && !has_bore && !has_bullet) ? "CrystalOuterTube" : "CrystalSolid";
    out << "    <!-- Active crystal: SteppingAction scores energy deposited here -->\n"
        << "    <volume name=\"active_crystal\">\n"
        << "      <materialref ref=\"" << crystal_mat_name << "\"/>\n"
        << "      <solidref ref=\"" << crystal_solid_ref << "\"/>\n";

    // The bore is part of the crystal polycone's profile, so the cavity belongs
    // to the world (vacuum for the validation exports -- which matches the MC,
    // where the bore region carries no material either).  It must NOT also be
    // filled with a daughter volume: that daughter would lie outside its mother
    // solid, and G4 navigation gets stuck on it.
    if (has_bore && is_cyl) {
        out << "\n      <!-- Bore cavity is part of the crystal polycone;"
            << " it is filled by the world material. -->\n";
    }

    out << "    </volume>\n\n";

    // Attenuator logical volumes.
    for (size_t i = 0; i < geom.attenuators().size(); ++i) {
        const auto& att = geom.attenuators()[i];
        std::string att_mat = att.material ? sanitize_name(att.material->name()) : "Air";
        out << "    <volume name=\"AttLV" << i << "\">\n"
            << "      <materialref ref=\"" << att_mat << "\"/>\n"
            << "      <solidref ref=\"AttSolid" << i << "\"/>\n"
            << "    </volume>\n\n";
    }

    // Source geometry logical volumes.
    if (has_source_geom) {
        if (source_geom->source_material()) {
            std::string src_mat = sanitize_name(source_geom->source_material()->name());
            out << "    <volume name=\"SrcMaterialLV\">\n"
                << "      <materialref ref=\"" << src_mat << "\"/>\n"
                << "      <solidref ref=\"SrcMaterialSolid\"/>\n"
                << "    </volume>\n\n";
        }

        const auto& shields = source_geom->shields();
        using SShape2 = SourceGeometry::Shape;
        if (source_geom->shape() == SShape2::Marinelli) {
            // Marinelli: separate logical volumes for each wall piece
            for (size_t i = 0; i < shields.size(); ++i) {
                std::string shield_mat = shields[i].material
                    ? sanitize_name(shields[i].material->name()) : "Air";
                out << "    <volume name=\"SrcShieldOuterLV" << i << "\">\n"
                    << "      <materialref ref=\"" << shield_mat << "\"/>\n"
                    << "      <solidref ref=\"SrcShieldOuter" << i << "\"/>\n"
                    << "    </volume>\n\n";
                out << "    <volume name=\"SrcShieldBottomLV" << i << "\">\n"
                    << "      <materialref ref=\"" << shield_mat << "\"/>\n"
                    << "      <solidref ref=\"SrcShieldBottom" << i << "\"/>\n"
                    << "    </volume>\n\n";
                out << "    <volume name=\"SrcShieldWellWallLV" << i << "\">\n"
                    << "      <materialref ref=\"" << shield_mat << "\"/>\n"
                    << "      <solidref ref=\"SrcShieldWellWall" << i << "\"/>\n"
                    << "    </volume>\n\n";
                out << "    <volume name=\"SrcShieldWellBottomLV" << i << "\">\n"
                    << "      <materialref ref=\"" << shield_mat << "\"/>\n"
                    << "      <solidref ref=\"SrcShieldWellBottom" << i << "\"/>\n"
                    << "    </volume>\n\n";
            }
        } else {
            for (size_t i = 0; i < shields.size(); ++i) {
                std::string shield_mat = shields[i].material
                    ? sanitize_name(shields[i].material->name()) : "Air";
                out << "    <volume name=\"SrcShieldLV" << i << "\">\n"
                    << "      <materialref ref=\"" << shield_mat << "\"/>\n"
                    << "      <solidref ref=\"SrcShieldSolid" << i << "\"/>\n"
                    << "    </volume>\n\n";
            }
        }
    }

    // World logical volume.
    const char* world_mat = vacuum_world ? "Vacuum" : "Air";
    out << "    <volume name=\"WorldLV\">\n"
        << "      <materialref ref=\"" << world_mat << "\"/>\n"
        << "      <solidref ref=\"WorldBox\"/>\n\n"
        << "      <!-- Crystal physical volume (front face at z=0, back at z=" << fmt(L,3) << " cm) -->\n"
        << "      <physvol name=\"active_crystal_PV\">\n"
        << "        <volumeref ref=\"active_crystal\"/>\n"
        << "        <positionref ref=\"crystal_center\"/>\n"
        << "      </physvol>\n";

    // Attenuator physical volumes.
    for (size_t i = 0; i < att_gdml.size(); ++i) {
        out << "\n      <physvol name=\"AttPV" << i << "\">\n"
            << "        <volumeref ref=\"AttLV" << i << "\"/>\n"
            << "        <position x=\"0\" y=\"0\" z=\"" << fmt(att_gdml[i].center_z) << "\" unit=\"cm\"/>\n"
            << "      </physvol>\n";
    }

    // Source geometry physical volumes.
    if (has_source_geom) {
        using SShape = SourceGeometry::Shape;

        if (source_geom->shape() == SShape::Marinelli) {
            // Marinelli: each piece placed at its own z-position (centered on z-axis).
            double z_bk = source_geom->marinelli_z_bk();
            double z_we = source_geom->marinelli_z_we();
            double z_bot = source_geom->marinelli_z_bot();
            double well_r = source_geom->marinelli_well_inner_radius();
            double total_h = z_we - z_bot;
            double fill_h = z_bk - z_bot;
            double ring_d = z_we - z_bk;

            // Source material (L-shaped solid): center of bounding outer cylinder.
            double src_z = (z_bot + z_we) / 2.0;
            if (source_geom->source_material()) {
                out << "\n      <physvol name=\"SrcMaterialPV\">\n"
                    << "        <volumeref ref=\"SrcMaterialLV\"/>\n"
                    << "        <position x=\"0\" y=\"0\" z=\"" << fmt(src_z) << "\" unit=\"cm\"/>\n"
                    << "      </physvol>\n";
            }

            // Beaker wall pieces (per shield layer).
            const auto& shields = source_geom->shields();
            for (size_t i = 0; i < shields.size(); ++i) {
                double t = shields[i].scalar_thickness();

                // Outer cylindrical wall: same center as source material
                out << "\n      <physvol name=\"SrcShieldOuterPV" << i << "\">\n"
                    << "        <volumeref ref=\"SrcShieldOuterLV" << i << "\"/>\n"
                    << "        <position x=\"0\" y=\"0\" z=\"" << fmt(src_z) << "\" unit=\"cm\"/>\n"
                    << "      </physvol>\n";

                // Bottom disk: at z_bot - t/2 (below the sample)
                double bottom_z = z_bot - t / 2.0;
                out << "      <physvol name=\"SrcShieldBottomPV" << i << "\">\n"
                    << "        <volumeref ref=\"SrcShieldBottomLV" << i << "\"/>\n"
                    << "        <position x=\"0\" y=\"0\" z=\"" << fmt(bottom_z) << "\" unit=\"cm\"/>\n"
                    << "      </physvol>\n";

                // Well inner wall: centered on the ring region
                double well_wall_z = z_bk + ring_d / 2.0;
                out << "      <physvol name=\"SrcShieldWellWallPV" << i << "\">\n"
                    << "        <volumeref ref=\"SrcShieldWellWallLV" << i << "\"/>\n"
                    << "        <position x=\"0\" y=\"0\" z=\"" << fmt(well_wall_z) << "\" unit=\"cm\"/>\n"
                    << "      </physvol>\n";

                // Well bottom annulus: at z_we + t/2 (caps the well from below)
                double well_bottom_z = z_we + t / 2.0;
                out << "      <physvol name=\"SrcShieldWellBottomPV" << i << "\">\n"
                    << "        <volumeref ref=\"SrcShieldWellBottomLV" << i << "\"/>\n"
                    << "        <position x=\"0\" y=\"0\" z=\"" << fmt(well_bottom_z) << "\" unit=\"cm\"/>\n"
                    << "      </physvol>\n";
            }
        } else {
            // Point / Cylindrical / Rectangular / Spherical sources.
            Eigen::Vector3d src_center;
            if (source_geom->shape() == SShape::Point) {
                src_center = source_geom->point_position();
            } else if (source_geom->shape() == SShape::Cylindrical) {
                src_center = source_geom->cyl_center();
            } else if (source_geom->shape() == SShape::Sphere) {
                src_center = source_geom->sphere_center();
            } else {
                src_center = source_geom->rect_center();
            }

            // Source material fill.
            if (source_geom->source_material() && source_geom->shape() != SShape::Point) {
                out << "\n      <physvol name=\"SrcMaterialPV\">\n"
                    << "        <volumeref ref=\"SrcMaterialLV\"/>\n"
                    << "        <position x=\"" << fmt(src_center.x()) << "\""
                    << " y=\"" << fmt(src_center.y()) << "\""
                    << " z=\"" << fmt(src_center.z()) << "\""
                    << " unit=\"cm\"/>\n"
                    << "      </physvol>\n";
            }

            // Shield layers.
            const auto& shields = source_geom->shields();
            for (size_t i = 0; i < shields.size(); ++i) {
                out << "\n      <physvol name=\"SrcShieldPV" << i << "\">\n"
                    << "        <volumeref ref=\"SrcShieldLV" << i << "\"/>\n"
                    << "        <position x=\"" << fmt(src_center.x()) << "\""
                    << " y=\"" << fmt(src_center.y()) << "\""
                    << " z=\"" << fmt(src_center.z()) << "\""
                    << " unit=\"cm\"/>\n"
                    << "      </physvol>\n";
            }
        }
    }

    out << "    </volume>\n";
    out << "  </structure>\n\n";

    // -----------------------------------------------------------------------
    // <setup>
    // -----------------------------------------------------------------------
    out << "  <setup name=\"Default\" version=\"1.0\">\n"
        << "    <world ref=\"WorldLV\"/>\n"
        << "  </setup>\n\n"
        << "</gdml>\n";

    out.flush();
}


void write_geant4_macro(const Eigen::Vector3d& source_pos,
                        double energy_keV,
                        uint64_t num_events,
                        const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("write_geant4_macro: cannot open file: " + filename);
    }

    // GEANT4 coordinate convention: our detector frame has z as detector axis,
    // with detector front face at z=0.  The source is typically at z < 0.
    // GEANT4 uses the same world frame we defined in the GDML (world origin =
    // our detector coordinate origin, which is the detector front face center).
    out << "# GEANT4 GPS macro generated by CeeLo\n"
        << "# Source at (" << fmt(source_pos.x(),4) << ", "
                           << fmt(source_pos.y(),4) << ", "
                           << fmt(source_pos.z(),4) << ") cm\n"
        << "# Energy: " << energy_keV << " keV\n"
        << "# Events: " << num_events << "\n\n"
        << "/gps/particle gamma\n"
        << "/gps/pos/type Point\n"
        << "/gps/pos/centre " << fmt(source_pos.x(),6) << " "
                              << fmt(source_pos.y(),6) << " "
                              << fmt(source_pos.z(),6) << " cm\n"
        << "/gps/ang/type iso\n"
        << "/gps/ene/type Mono\n"
        << "/gps/ene/mono " << fmt(energy_keV * 1e-3, 6) << " MeV\n\n"
        << "/run/beamOn " << num_events << "\n";

    out.flush();
}

void write_geant4_macro_marinelli(const SourceGeometry& source_geom,
                                   double energy_keV,
                                   uint64_t num_events,
                                   const std::string& filename)
{
    if (source_geom.shape() != SourceGeometry::Shape::Marinelli) {
        throw std::runtime_error(
            "write_geant4_macro_marinelli: source geometry is not Marinelli");
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error(
            "write_geant4_macro_marinelli: cannot open file: " + filename);
    }

    double z_bk = source_geom.marinelli_z_bk();
    double z_we = source_geom.marinelli_z_we();
    double z_bot = source_geom.marinelli_z_bot();
    double outer_r = source_geom.marinelli_outer_radius();
    double total_h = z_we - z_bot;
    double half_z = total_h / 2.0;
    double center_z = (z_bot + z_we) / 2.0;

    out << "# GEANT4 GPS macro generated by CeeLo (Marinelli volumetric source)\n"
        << "# Bounding cylinder: R=" << fmt(outer_r, 4) << " cm, halfZ="
        << fmt(half_z, 4) << " cm, center z=" << fmt(center_z, 4) << " cm\n"
        << "# Confined to SrcMaterialPV (L-shaped sample volume)\n"
        << "# Energy: " << energy_keV << " keV\n"
        << "# Events: " << num_events << "\n\n"
        << "/gps/particle gamma\n"
        << "/gps/pos/type Volume\n"
        << "/gps/pos/shape Cylinder\n"
        << "/gps/pos/radius " << fmt(outer_r, 6) << " cm\n"
        << "/gps/pos/halfz " << fmt(half_z, 6) << " cm\n"
        << "/gps/pos/centre 0 0 " << fmt(center_z, 6) << " cm\n"
        << "/gps/pos/confine SrcMaterialPV\n"
        << "/gps/ang/type iso\n"
        << "/gps/ene/type Mono\n"
        << "/gps/ene/mono " << fmt(energy_keV * 1e-3, 6) << " MeV\n\n"
        << "/run/beamOn " << num_events << "\n";

    out.flush();
}

void write_geant4_macro_cylindrical(const Eigen::Vector3d& center,
                                     double radius, double half_length,
                                     double energy_keV,
                                     uint64_t num_events,
                                     const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error(
            "write_geant4_macro_cylindrical: cannot open file: " + filename);
    }

    out << "# GEANT4 GPS macro generated by CeeLo (cylindrical volume source)\n"
        << "# Cylinder: R=" << fmt(radius, 4) << " cm, halfZ="
        << fmt(half_length, 4) << " cm, center z=" << fmt(center.z(), 4) << " cm\n"
        << "# Confined to SrcMaterialPV\n"
        << "# Energy: " << energy_keV << " keV\n"
        << "# Events: " << num_events << "\n\n"
        << "/gps/particle gamma\n"
        << "/gps/pos/type Volume\n"
        << "/gps/pos/shape Cylinder\n"
        << "/gps/pos/radius " << fmt(radius, 6) << " cm\n"
        << "/gps/pos/halfz " << fmt(half_length, 6) << " cm\n"
        << "/gps/pos/centre " << fmt(center.x(), 6) << " "
                              << fmt(center.y(), 6) << " "
                              << fmt(center.z(), 6) << " cm\n"
        << "/gps/pos/confine SrcMaterialPV\n"
        << "/gps/ang/type iso\n"
        << "/gps/ene/type Mono\n"
        << "/gps/ene/mono " << fmt(energy_keV * 1e-3, 6) << " MeV\n\n"
        << "/run/beamOn " << num_events << "\n";

    out.flush();
}

void write_geant4_macro_rectangular(const Eigen::Vector3d& center,
                                     const Eigen::Vector3d& half_dims,
                                     double energy_keV,
                                     uint64_t num_events,
                                     const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error(
            "write_geant4_macro_rectangular: cannot open file: " + filename);
    }

    out << "# GEANT4 GPS macro generated by CeeLo (rectangular volume source)\n"
        << "# Box half-dims: " << fmt(half_dims.x(), 4) << " x "
        << fmt(half_dims.y(), 4) << " x " << fmt(half_dims.z(), 4) << " cm\n"
        << "# Center: (" << fmt(center.x(), 4) << ", "
        << fmt(center.y(), 4) << ", " << fmt(center.z(), 4) << ") cm\n"
        << "# Confined to SrcMaterialPV\n"
        << "# Energy: " << energy_keV << " keV\n"
        << "# Events: " << num_events << "\n\n"
        << "/gps/particle gamma\n"
        << "/gps/pos/type Volume\n"
        << "/gps/pos/shape Para\n"
        << "/gps/pos/halfx " << fmt(half_dims.x(), 6) << " cm\n"
        << "/gps/pos/halfy " << fmt(half_dims.y(), 6) << " cm\n"
        << "/gps/pos/halfz " << fmt(half_dims.z(), 6) << " cm\n"
        << "/gps/pos/centre " << fmt(center.x(), 6) << " "
                              << fmt(center.y(), 6) << " "
                              << fmt(center.z(), 6) << " cm\n"
        << "/gps/pos/confine SrcMaterialPV\n"
        << "/gps/ang/type iso\n"
        << "/gps/ene/type Mono\n"
        << "/gps/ene/mono " << fmt(energy_keV * 1e-3, 6) << " MeV\n\n"
        << "/run/beamOn " << num_events << "\n";

    out.flush();
}

void write_geant4_macro_spherical(const Eigen::Vector3d& center,
                                   double radius,
                                   double energy_keV,
                                   uint64_t num_events,
                                   const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error(
            "write_geant4_macro_spherical: cannot open file: " + filename);
    }

    out << "# GEANT4 GPS macro generated by CeeLo (spherical volume source)\n"
        << "# Sphere: R=" << fmt(radius, 4) << " cm, center=("
        << fmt(center.x(), 4) << ", " << fmt(center.y(), 4) << ", "
        << fmt(center.z(), 4) << ") cm\n"
        << "# Confined to SrcMaterialPV (a hollow shell rejects the void core)\n"
        << "# Energy: " << energy_keV << " keV\n"
        << "# Events: " << num_events << "\n\n"
        << "/gps/particle gamma\n"
        << "/gps/pos/type Volume\n"
        << "/gps/pos/shape Sphere\n"
        << "/gps/pos/radius " << fmt(radius, 6) << " cm\n"
        << "/gps/pos/centre " << fmt(center.x(), 6) << " "
                              << fmt(center.y(), 6) << " "
                              << fmt(center.z(), 6) << " cm\n"
        << "/gps/pos/confine SrcMaterialPV\n"
        << "/gps/ang/type iso\n"
        << "/gps/ene/type Mono\n"
        << "/gps/ene/mono " << fmt(energy_keV * 1e-3, 6) << " MeV\n\n"
        << "/run/beamOn " << num_events << "\n";

    out.flush();
}

} // namespace ceelo
