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

/// @file main.cc
/// @brief GEANT4 validation harness for CeeLo.
///
/// Usage:
///   ceelo_g4val <detector.gdml> <run.mac> <output.csv> [options]
///
/// Options:
///   --histogram      Write _histogram.csv and _summary.csv diagnostic files
///   --sphere-flux R  Write _sphereflux.csv: photon TRACK-LENGTH flux (phi=L/V)
///                    inside the scoring sphere (the volume named
///                    active_crystal), binned in (E, mu), all generations,
///                    split direct vs scattered. This is the observable CeeLo's
///                    NextEventEstimator produces -- use it for environment
///                    (air/ground/room) scatter validation.
///   --fluence        Write _fluence.csv: photon fluence incident on the detector
///                    (uncollided peak + air-scattered continuum). Air-scatter study.
///   --cone-bias      Enable cone importance sampling (matches our MC's approach)
///                    Parses source position and energy from the macro file,
///                    detector geometry from the GDML.
///   --cone-scale <f> Widen the auto cone half-angle by factor f (air-scatter
///                    study: capture off-axis air scatter-in). Requires --cone-bias.
///                    Omit --cone-bias entirely for fully isotropic emission.
///   --em-standard    Use G4EmStandardPhysics (default EM, no Doppler broadening)
///   --em-option4     Use option4 as a systematic cross-check instead of the
///                    FullSpecApprox campaign default (Livermore EM)
///   --component-labels  Add zero/FEP/escape/Compton/pair weighted histograms
///   --seed N         Deterministic Geant4 random seed
///   --diagnostics    Write per-interaction diagnostic histograms (11 CSVs)
///   --entry-diag     Write crystal entry diagnostic CSVs (energy, z-pos, summary)
///   --fep-window <keV>  Full-energy half-window. Default is CeeLo's
///                    kDefaultFepWindowKeV; pass 1.5 to match the references
///                    tests/data/geant4_reference/README.md lists at 1.5 keV.
///   --lowcut         Set production cut to 0.001 mm (default 0.7 mm) for thin
///                    detector studies where default cuts may suppress physics
///   --pixe           Particle-induced x-ray emission: electron-impact
///                    ionization x-rays at every energy (without it option4
///                    makes them only below 100 keV, through PenIoni)
///   --fep-precision <target>  Run until FEP relative precision <= target
///                    (e.g., 0.01 for 1%). Safety cap: 100M events.
///   --vis            Open interactive OpenGL visualization (no simulation).
///                    The macro file should contain /vis/ commands (e.g., vis.mac).
///   --vis-batch      Same, but headless: runs the /vis/ macro and exits without
///                    opening a session, so it works over ssh and in CI. Use a
///                    file-writing driver that needs no tracking -- VRML2FILE
///                    dumps G4's own polyhedron for every solid, which is the
///                    honest check that G4 read a solid the way the exporter
///                    intended:
///                      /vis/open VRML2FILE
///                      /vis/drawVolume
///                      /vis/viewer/flush      -> writes g4_NN.wrl
///                    See macros/vis_crystal_vrml.mac, and render_vrml.py to
///                    turn the .wrl into a PNG.  RayTracer does NOT work from
///                    this executable -- see the run manager below for why.
///
/// Physics list: G4EmLivermorePhysics by default for formal environment anchors.

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "io/EnvFitArchive.h"

#include "FTFP_BERT.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <exception>

/// Parse source position (cm) and energy (MeV) from GPS macro.
static bool parse_macro(const std::string& macro_file,
                        G4ThreeVector& source_pos_mm,
                        double& energy_MeV,
                        int& num_events)
{
    std::ifstream mac(macro_file);
    if (!mac) return false;

    bool got_pos = false, got_energy = false, got_events = false;
    std::string line;
    while (std::getline(mac, line)) {
        if (line.find("/gps/pos/centre") != std::string::npos) {
            std::istringstream iss(line);
            std::string cmd;
            double x, y, z;
            std::string unit;
            if (iss >> cmd >> x >> y >> z >> unit) {
                double scale = (unit == "mm") ? 1.0
                             : (unit == "m") ? 1000.0 : 10.0;
                source_pos_mm = G4ThreeVector(x * scale, y * scale, z * scale);
                got_pos = true;
            }
        }
        if (line.find("/gps/ene/mono") != std::string::npos ||
            line.find("/gps/energy") != std::string::npos) {
            std::istringstream iss(line);
            std::string cmd;
            double val;
            std::string unit;
            if (iss >> cmd >> val >> unit) {
                energy_MeV = (unit == "keV") ? val * 0.001
                           : (unit == "eV") ? val * 1.0e-6 : val;
                got_energy = true;
            }
        }
        if (line.find("/run/beamOn") != std::string::npos) {
            std::istringstream iss(line);
            std::string cmd;
            int n;
            if (iss >> cmd >> n) {
                num_events = n;
                got_events = true;
            }
        }
    }
    return got_pos && got_energy && got_events;
}

struct CeeFitGeometryMetadata {
    std::string geometry_tier;
    std::string material;
    double energy_keV = 0.0;
    double distance_cm = 0.0;
    double source_height_cm = 0.0;
    double detector_height_cm = 0.0;
    double scoring_radius_cm = 0.0;
    double world_half_extent_cm = 0.0;
};

static bool parse_ceefit_geometry_metadata(
    const std::string& gdml_file, CeeFitGeometryMetadata& result)
{
    std::ifstream source(gdml_file);
    if (!source) return false;
    const std::string content((std::istreambuf_iterator<char>(source)),
                              std::istreambuf_iterator<char>());
    constexpr const char* marker = "CEEFIT-GEOMETRY-v1";
    const std::size_t begin = content.find(marker);
    if (begin == std::string::npos) return false;
    const std::size_t end = content.find("-->", begin);
    if (end == std::string::npos) return false;
    std::istringstream fields(
        content.substr(begin + std::char_traits<char>::length(marker),
                       end - begin));
    bool tier = false, material = false, energy = false, distance = false;
    bool source_height = false, detector_height = false, radius = false;
    bool world = false;
    std::string token;
    try {
        while (fields >> token) {
            const std::size_t equals = token.find('=');
            if (equals == std::string::npos) continue;
            const std::string key = token.substr(0, equals);
            const std::string value = token.substr(equals + 1);
            if (key == "geometry_tier") {
                result.geometry_tier = value; tier = true;
            } else if (key == "material") {
                result.material = value; material = true;
            } else if (key == "energy_keV") {
                result.energy_keV = std::stod(value); energy = true;
            } else if (key == "distance_cm") {
                result.distance_cm = std::stod(value); distance = true;
            } else if (key == "source_height_cm") {
                result.source_height_cm = std::stod(value); source_height = true;
            } else if (key == "detector_height_cm") {
                result.detector_height_cm = std::stod(value); detector_height = true;
            } else if (key == "scoring_radius_cm") {
                result.scoring_radius_cm = std::stod(value); radius = true;
            } else if (key == "world_half_extent_cm") {
                result.world_half_extent_cm = std::stod(value); world = true;
            }
        }
    } catch (const std::exception&) {
        return false;
    }
    return tier && material && energy && distance && source_height
        && detector_height && radius && world;
}

static std::string file_sha256(const std::string& path)
{
    std::ifstream source(path, std::ios::binary);
    if (!source) return "";
    const std::string bytes((std::istreambuf_iterator<char>(source)),
                            std::istreambuf_iterator<char>());
    return ceelo::env_fit_sha256_hex(bytes.data(), bytes.size());
}

static bool close_metadata(double left, double right)
{
    return std::isfinite(left) && std::isfinite(right)
        && std::abs(left - right)
               <= 1.0e-9 * std::max({std::abs(left), std::abs(right), 1.0});
}

static std::string xml_attribute(const std::string& tag,
                                 const std::string& attribute)
{
    const std::string marker = attribute + "=\"";
    std::size_t begin = tag.find(marker);
    if (begin == std::string::npos) return "";
    begin += marker.size();
    const std::size_t end = tag.find('"', begin);
    return end == std::string::npos ? "" : tag.substr(begin, end - begin);
}

static double length_to_m(double value, const std::string& unit)
{
    if (unit == "m") return value;
    if (unit == "cm" || unit.empty()) return value * 0.01;
    if (unit == "mm") return value * 0.001;
    return 0.0;
}

/// Read the actual world solid from the source GDML.  Box dimensions are full
/// lengths, so the smallest half-length is the containment limit.
static bool get_world_half_extent_m_from_gdml(const std::string& gdml_file,
                                              double& half_extent_m)
{
    std::ifstream source(gdml_file);
    if (!source) return false;
    const std::string content((std::istreambuf_iterator<char>(source)),
                              std::istreambuf_iterator<char>());
    std::size_t position = 0;
    while ((position = content.find("<box", position)) != std::string::npos) {
        const std::size_t end = content.find("/>", position);
        if (end == std::string::npos) return false;
        const std::string tag = content.substr(position, end - position);
        const std::string name = xml_attribute(tag, "name");
        if (name.find("world") != std::string::npos ||
            name.find("World") != std::string::npos) {
            try {
                const double x = std::stod(xml_attribute(tag, "x"));
                const double y = std::stod(xml_attribute(tag, "y"));
                const double z = std::stod(xml_attribute(tag, "z"));
                const double scale = length_to_m(1.0, xml_attribute(tag, "lunit"));
                half_extent_m = 0.5 * std::min({x, y, z}) * scale;
                return std::isfinite(half_extent_m) && half_extent_m > 0.0;
            } catch (const std::exception&) {
                return false;
            }
        }
        position = end + 2;
    }
    return false;
}

/// Confirm that the volume used by the sphere scorer has the radius supplied
/// for track-length/volume normalization.
static bool get_scoring_radius_cm_from_gdml(const std::string& gdml_file,
                                            double& radius_cm)
{
    std::ifstream source(gdml_file);
    if (!source) return false;
    const std::string content((std::istreambuf_iterator<char>(source)),
                              std::istreambuf_iterator<char>());
    std::size_t position = 0;
    while ((position = content.find("<sphere", position)) != std::string::npos) {
        const std::size_t end = content.find("/>", position);
        if (end == std::string::npos) return false;
        const std::string tag = content.substr(position, end - position);
        const std::string name = xml_attribute(tag, "name");
        if (name.find("score") != std::string::npos ||
            name.find("Score") != std::string::npos) {
            try {
                const double rmax = std::stod(xml_attribute(tag, "rmax"));
                const double metres = length_to_m(rmax, xml_attribute(tag, "lunit"));
                radius_cm = metres * 100.0;
                return std::isfinite(radius_cm) && radius_cm > 0.0;
            } catch (const std::exception&) {
                return false;
            }
        }
        position = end + 2;
    }
    return false;
}

static bool macro_declares_isotropic_gps(const std::string& macro_file)
{
    std::ifstream source(macro_file);
    if (!source) return false;
    bool found = false;
    std::string line;
    while (std::getline(source, line)) {
        std::istringstream fields(line);
        std::string command;
        std::string value;
        if (!(fields >> command)) continue;
        if (command.rfind("/random/", 0) == 0) return false;
        if (!(fields >> value) || command != "/gps/ang/type") continue;
        if (value != "iso" || found) return false;
        found = true;
    }
    return found;
}

/// Get bounding radius and z-extent of the outermost solid in the GDML.
/// Looks for the "world" or largest volume.  For our GDML, the active_crystal
/// or outermost attenuator defines the detector extent.
static bool get_detector_bounds_from_gdml(const std::string& gdml_file,
                                          double& radius_mm,
                                          double& z_min_mm,
                                          double& z_max_mm)
{
    // Parse GDML to find outermost detector solid dimensions.
    // Our GDML convention: detector front face at z=0, extends to z=L.
    // We search G4LogicalVolumeStore for volumes (after Construct()) to find
    // the bounding radius.  But at this point Construct() hasn't been called.
    //
    // Alternative: re-parse the GDML ourselves looking for solid dimensions.
    // Our GDML always has an active_crystal solid.  For cylinder detectors,
    // it's a G4Tubs.  For attenuators, they're subtraction solids.
    //
    // Simplest approach: scan the GDML file for the tube/box parameters
    // of the outermost solid.
    std::ifstream f(gdml_file);
    if (!f) return false;

    double max_r = 0.0;
    double max_z = 0.0;
    double min_z = 0.0;

    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    // Exports since Oct 2026 carry the geometry's outer bounding cylinder
    //  (the one CeeLo's own cone is sized from).  Prefer it: the solid scan
    //  below takes z_min from "att_" positions the exporter no longer writes,
    //  so an endcap or clutter standing in front of the crystal face was left
    //  outside the cone.
    {
        const size_t b = content.find("ceelo_cone_bounds_cm");
        if (b != std::string::npos) {
            const std::string tag = content.substr(b, content.find("-->", b) - b);
            auto attr = [&tag](const std::string& key) {
                const size_t p = tag.find(key + "=\"");
                return std::stod(tag.substr(p + key.size() + 2));
            };
            radius_mm = attr("r") * 10.0;
            z_min_mm  = attr("zmin") * 10.0;
            z_max_mm  = attr("zmax") * 10.0;
            return radius_mm > 0.0 && z_max_mm > z_min_mm;
        }
    }

    // Helper to parse a named attribute from an XML tag string
    auto parse_attr = [](const std::string& tag, const std::string& attr) -> double {
        size_t p = tag.find(attr + "=\"");
        if (p == std::string::npos) return 0.0;
        p += attr.size() + 2;
        return std::stod(tag.substr(p));
    };

    // Helper to parse the "name" attribute from an XML tag string
    auto parse_name = [](const std::string& tag) -> std::string {
        size_t p = tag.find("name=\"");
        if (p == std::string::npos) return "";
        p += 6;
        size_t q = tag.find('"', p);
        if (q == std::string::npos) return "";
        return tag.substr(p, q - p);
    };

    // Find all <tube> elements and get the largest rmax.
    // Skip the WorldBox and any world-related solids.
    // Format: <tube name="..." rmin="..." rmax="3.810000" z="7.620000" .../>
    // Units are cm in our GDML
    size_t pos = 0;
    while ((pos = content.find("<tube", pos)) != std::string::npos) {
        size_t end = content.find("/>", pos);
        if (end == std::string::npos) break;
        std::string tag = content.substr(pos, end - pos);

        std::string name = parse_name(tag);
        // Skip world-related solids
        if (name.find("World") != std::string::npos ||
            name.find("world") != std::string::npos) {
            pos = end + 2;
            continue;
        }

        double rmax = parse_attr(tag, "rmax");
        double z = parse_attr(tag, "z");

        if (rmax > max_r) max_r = rmax;
        if (z > max_z) max_z = z;

        pos = end + 2;
    }

    // Also check <polycone> elements.  Coaxial and/or bulletized crystals are
    // exported as a single polycone (a surface of revolution) rather than a
    // tube, so the bounds come from the widest rmax and the z-plane span.
    // Unlike <tube>, whose z attribute is a full height about its own centre,
    // polycone z-planes are already in the solid's local frame.
    pos = 0;
    while ((pos = content.find("<polycone", pos)) != std::string::npos) {
        size_t end = content.find("</polycone>", pos);
        if (end == std::string::npos) break;
        std::string block = content.substr(pos, end - pos);

        std::string name = parse_name(block.substr(0, block.find('>')));
        if (name.find("World") != std::string::npos ||
            name.find("world") != std::string::npos) {
            pos = end + 1;
            continue;
        }

        double z_lo = 0.0, z_hi = 0.0;
        bool first = true;
        size_t zp = 0;
        while ((zp = block.find("<zplane", zp)) != std::string::npos) {
            size_t zend = block.find("/>", zp);
            if (zend == std::string::npos) break;
            std::string tag = block.substr(zp, zend - zp);

            double rmax = parse_attr(tag, "rmax");
            double z    = parse_attr(tag, "z");
            if (rmax > max_r) max_r = rmax;
            if (first) { z_lo = z_hi = z; first = false; }
            else { z_lo = std::min(z_lo, z); z_hi = std::max(z_hi, z); }

            zp = zend + 2;
        }
        // Express as a full height, matching how <tube> z is treated above.
        if (!first && (z_hi - z_lo) > max_z) max_z = z_hi - z_lo;

        pos = end + 1;
    }

    // Also check <box> elements for box detectors.
    // Skip the WorldBox solid.
    pos = 0;
    while ((pos = content.find("<box", pos)) != std::string::npos) {
        size_t end = content.find("/>", pos);
        if (end == std::string::npos) break;
        std::string tag = content.substr(pos, end - pos);

        std::string name = parse_name(tag);
        // Skip world-related solids
        if (name.find("World") != std::string::npos ||
            name.find("world") != std::string::npos) {
            pos = end + 2;
            continue;
        }

        double x = parse_attr(tag, "x");
        double y = parse_attr(tag, "y");
        double z = parse_attr(tag, "z");

        // Box half-diagonal as bounding radius
        double diag = std::sqrt(x * x + y * y) / 2.0;
        if (diag > max_r) max_r = diag;
        if (z > max_z) max_z = z;

        pos = end + 2;
    }

    if (max_r <= 0.0 || max_z <= 0.0) return false;

    // Our GDML convention: detector front face at z=0, extends to z=L.
    // Attenuators may extend slightly in front (z < 0).
    // Parse crystal_center position to find z offset.
    // For simplicity, use z_min=0 for bare detectors.

    // Check for attenuator z positions
    // Look for position elements
    pos = 0;
    double att_z_min = 0.0;
    while ((pos = content.find("<position", pos)) != std::string::npos) {
        size_t end = content.find("/>", pos);
        if (end == std::string::npos) break;
        std::string tag = content.substr(pos, end - pos);

        auto parse_attr = [&](const std::string& attr) -> double {
            size_t p = tag.find(attr + "=\"");
            if (p == std::string::npos) return 0.0;
            p += attr.size() + 2;
            return std::stod(tag.substr(p));
        };

        double z = parse_attr("z");
        // Attenuator positions can be negative
        if (z < att_z_min && tag.find("att_") != std::string::npos) {
            att_z_min = z;
        }

        pos = end + 2;
    }

    // Convert cm to mm (GDML uses cm, G4 internal uses mm)
    radius_mm = max_r * 10.0;
    z_min_mm = att_z_min * 10.0;  // could be negative for attenuators
    z_max_mm = max_z * 10.0;

    return true;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <detector.gdml> <run.mac> <output.csv>"
                  << " [--histogram] [--cone-bias] [--bin-width <keV>]\n";
        return 1;
    }

    std::string gdml_file   = argv[1];
    std::string macro_file  = argv[2];
    std::string output_csv  = argv[3];

    bool histogram_mode = false;
    bool cone_bias = false;
    bool em_standard = false;
    bool em_option4 = false;
    bool component_labels = false;
    bool diagnostics_mode = false;
    bool entry_diag_mode = false;
    bool low_cut = false;
    bool pixe = false;
    bool vis_mode = false;
    bool fluence_mode = false; // --fluence: photon-fluence-at-detector scorer
    // --sphere-flux <R_cm>: environment-scatter track-length flux in a scoring
    // sphere (the volume named active_crystal), binned in (E, mu).
    double sphere_flux_radius_cm = 0.0;
    double cone_scale = 1.0;   // --cone-scale: widen the auto cone (air-scatter)
    bool vis_batch = false;
    double bin_width_keV = 1.0;
    double fep_precision_target = 0.0; // 0 = disabled
    double fep_window_keV = 0.0;       // 0 = use ceelo::kDefaultFepWindowKeV
    bool correlated_gamma = false; // --correlated-gamma: enable G4 gamma-gamma angular correlation
    long random_seed = 246813579;
    for (int i = 4; i < argc; ++i) {
        if (std::string(argv[i]) == "--histogram") histogram_mode = true;
        if (std::string(argv[i]) == "--cone-bias") cone_bias = true;
        if (std::string(argv[i]) == "--fluence") fluence_mode = true;
        else if (std::string(argv[i]) == "--sphere-flux" && i + 1 < argc)
            sphere_flux_radius_cm = std::atof(argv[++i]);
        if (std::string(argv[i]) == "--cone-scale" && i + 1 < argc) {
            cone_scale = std::atof(argv[++i]);
        }
        if (std::string(argv[i]) == "--em-standard") em_standard = true;
        if (std::string(argv[i]) == "--em-option4") em_option4 = true;
        if (std::string(argv[i]) == "--component-labels") component_labels = true;
        if (std::string(argv[i]) == "--diagnostics") diagnostics_mode = true;
        if (std::string(argv[i]) == "--entry-diag") entry_diag_mode = true;
        if (std::string(argv[i]) == "--lowcut") low_cut = true;
        if (std::string(argv[i]) == "--pixe") pixe = true;
        if (std::string(argv[i]) == "--vis") vis_mode = true;
        if (std::string(argv[i]) == "--vis-batch") vis_batch = true;
        if (std::string(argv[i]) == "--correlated-gamma") correlated_gamma = true;
        if (std::string(argv[i]) == "--bin-width" && i + 1 < argc) {
            bin_width_keV = std::atof(argv[++i]);
        }
        if (std::string(argv[i]) == "--fep-precision" && i + 1 < argc) {
            fep_precision_target = std::atof(argv[++i]);
        }
        if (std::string(argv[i]) == "--fep-window" && i + 1 < argc) {
            fep_window_keV = std::atof(argv[++i]);
        }
        if (std::string(argv[i]) == "--seed" && i + 1 < argc) {
            random_seed = std::stol(argv[++i]);
        }
    }

    // Build run manager (Serial mode).
    //
    // NOTE for visualization: do NOT use the RayTracer driver from this
    // executable.  In an MT-enabled GEANT4, G4RayTracerViewer unconditionally
    // instantiates G4TheMTRayTracer, whose StoreUserActions() dereferences
    // G4RunManagerFactory::GetMTMasterRunManager() -- null here -- and
    // segfaults; switching this manager to MT only trades that for RayTracer's
    // worker threads having no primary generator, since this harness sets user
    // actions directly rather than through a G4VUserActionInitialization.  Use
    // the file-writing drivers (VRML2FILE, HepRepFile) instead: they traverse
    // the geometry and write G4's own polyhedron for each solid without
    // tracking anything.  See macros/vis_crystal_vrml.mac.
    auto* run_manager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Serial);
    G4Random::setTheSeed(random_seed);
    std::cout << "Geant4 random seed: " << random_seed << "\n";

    // Detector geometry from exported GDML.
    auto* det = new DetectorConstruction(gdml_file);
    run_manager->SetUserInitialization(det);

    // Physics list selection.
    auto* physics = new FTFP_BERT;
    std::string physics_list_name = "G4EmLivermorePhysics";
    if (em_standard) {
        physics->ReplacePhysics(new G4EmStandardPhysics());
        physics_list_name = "G4EmStandardPhysics";
        std::cout << "Using G4EmStandardPhysics (no Doppler broadening)\n";
    } else if (em_option4) {
        physics->ReplacePhysics(new G4EmStandardPhysics_option4());
        physics_list_name = "G4EmStandardPhysics_option4";
        std::cout << "Using G4EmStandardPhysics_option4 systematic cross-check\n";
    } else {
        physics->ReplacePhysics(new G4EmLivermorePhysics());
        std::cout << "Using G4EmLivermorePhysics (campaign default)\n";
    }
    // Radioactive decay: enables `/gps/particle ion` full-decay cascade sources
    // for true-coincidence-summing validation. Harmless for photon primaries
    // (they do not radioactively decay), so always registered.
    physics->RegisterPhysics(new G4RadioactiveDecayPhysics());
    run_manager->SetUserInitialization(physics);

    // Configure atomic deexcitation: fluorescence + Auger, ignore production cuts.
    // Without IgnoreCut, G4 suppresses fluorescence photons whose range in the
    // creation material is below the production threshold (default 0.7 mm).
    // Pb Kα (75 keV) has ~0.15 mm range in Pb, so fluorescence from Pb
    // attenuators would be incorrectly suppressed.
    auto* em_params = G4EmParameters::Instance();
    em_params->SetFluo(true);
    em_params->SetAuger(true);
    em_params->SetDeexcitationIgnoreCut(true);
    if (pixe) {
        em_params->SetPixe(true);
        std::cout << "PIXE on (electron-impact ionization x-rays at every energy)\n";
    }

    // Low production cut for thin detector studies (e.g., CZT).
    // Default 0.7 mm may suppress secondaries relevant to thin detectors.
    if (low_cut) {
        physics->SetDefaultCutValue(0.001 * mm);
        std::cout << "Low production cut: 0.001 mm (default 0.7 mm)\n";
    }

    // User actions.
    auto* run_action  = new RunAction(output_csv, gdml_file, macro_file);
    auto* primary_gen = new PrimaryGeneratorAction();
    auto* event_action = new EventAction(run_action, primary_gen);
    auto* step_action  = new SteppingAction(event_action, det->GetActiveCrystalLV());

    // Full-energy window.  Defaults to ceelo::kDefaultFepWindowKeV so the harness
    //  and CeeLo agree; pass --fep-window to regenerate a reference set at
    //  another window (most committed references are scored at 1.5 keV -- see
    //  tests/data/geant4_reference/README.md).
    if (fep_window_keV > 0.0) {
        event_action->SetFepWindowKeV(fep_window_keV);
    }
    std::cout << "FEP window: " << event_action->GetFepWindowKeV()
              << " keV half-width\n";

    // Enable histogram mode if requested.
    if (histogram_mode) {
        run_action->EnableHistogram(true);
        G4ThreeVector parsed_position;
        double parsed_energy_MeV = 12.0;
        int parsed_events = 0;
        (void)parse_macro(
            macro_file, parsed_position, parsed_energy_MeV, parsed_events);
        run_action->SetMaxEnergy(1.05 * parsed_energy_MeV);
        run_action->SetBinWidth(bin_width_keV);
        event_action->EnableHistogram(true);
        std::cout << "Histogram mode enabled (bin width: " << bin_width_keV << " keV)\n";
    }

    // Enable fluence-at-detector scorer if requested.
    if (fluence_mode) {
        run_action->EnableFluence(true);
        run_action->SetBinWidth(bin_width_keV);
        event_action->EnableFluence(true);
        step_action->EnableFluence(true);
        std::cout << "Fluence-at-detector scorer enabled (bin width: "
                  << bin_width_keV << " keV)\n";
    }

    if (component_labels) {
        if (!histogram_mode) {
            std::cerr << "ERROR: --component-labels requires --histogram\n";
            return 1;
        }
        run_action->EnableComponentLabels(true);
        event_action->EnableComponentLabels(true);
        step_action->EnableDiagnostics(true);
        std::cout << "Weighted component labels enabled\n";
    }

    // Enable per-interaction diagnostics if requested.
    if (diagnostics_mode) {
        run_action->EnableDiagnostics(true);
        event_action->EnableDiagnostics(true);
        step_action->EnableDiagnostics(true);
        std::cout << "Per-interaction diagnostics enabled\n";
    }

    // Enable crystal entry diagnostics if requested.
    if (entry_diag_mode) {
        run_action->EnableEntryDiag(true);
        event_action->EnableEntryDiag(true);
        step_action->EnableEntryDiag(true);
        // Crystal dimensions for surface classification (3"x3" NaI default)
        // TODO: parse from GDML if needed for other detector sizes
        step_action->SetCrystalDimensions(7.62, 3.81);
        std::cout << "Crystal entry diagnostics enabled\n";
    }

    // Enable cone bias if requested.
    if (cone_bias) {
        G4ThreeVector source_pos_mm;
        double energy_MeV = 0.0;
        int num_events = 0;

        if (!parse_macro(macro_file, source_pos_mm, energy_MeV, num_events)) {
            std::cerr << "ERROR: --cone-bias requires macro with "
                      << "/gps/pos/centre, /gps/ene/mono, /run/beamOn\n";
            return 1;
        }

        double det_r_mm = 0.0, det_z_min_mm = 0.0, det_z_max_mm = 0.0;
        if (!get_detector_bounds_from_gdml(gdml_file, det_r_mm,
                                            det_z_min_mm, det_z_max_mm)) {
            std::cerr << "ERROR: --cone-bias could not parse detector bounds "
                      << "from GDML\n";
            return 1;
        }

        primary_gen->SetConeScale(cone_scale);
        primary_gen->SetDirectionSeed(static_cast<std::uint64_t>(random_seed));
        primary_gen->EnableConeBias(source_pos_mm, energy_MeV,
                                    det_r_mm, det_z_min_mm, det_z_max_mm);

        run_action->SetConeSolidAngleFraction(
            primary_gen->GetConeSolidAngleFraction());

        std::cout << "Cone bias enabled:\n"
                  << "  Source: (" << source_pos_mm.x() << ", "
                  << source_pos_mm.y() << ", "
                  << source_pos_mm.z() << ") mm\n"
                  << "  Energy: " << energy_MeV << " MeV\n"
                  << "  Det R:  " << det_r_mm << " mm\n"
                  << "  Det Z:  [" << det_z_min_mm << ", "
                  << det_z_max_mm << "] mm\n";
    }

    // Environment-scatter sphere flux scorer.
    if (sphere_flux_radius_cm > 0.0) {
        if (physics_list_name != "G4EmLivermorePhysics") {
            std::cerr << "ERROR: formal --sphere-flux anchors require "
                         "G4EmLivermorePhysics; use another scorer for EM "
                         "systematic studies\n";
            delete run_manager;
            return 1;
        }
        if (random_seed <= 0) {
            std::cerr << "ERROR: formal --sphere-flux anchors require a positive "
                         "explicitly recorded RNG seed\n";
            delete run_manager;
            return 1;
        }
        const double production_cut_mm =
            physics->GetDefaultCutValue() / CLHEP::mm;
        if (!std::isfinite(production_cut_mm) || production_cut_mm <= 0.0
            || production_cut_mm > 0.7) {
            std::cerr << "ERROR: formal --sphere-flux production cut must be in "
                         "(0, 0.7] mm\n";
            delete run_manager;
            return 1;
        }
        if (cone_bias) {
            std::cerr << "ERROR: --sphere-flux requires unbiased isotropic GPS; "
                         "the sphere scorer does not apply cone event weights\n";
            delete run_manager;
            return 1;
        }
        if (fep_precision_target > 0.0 || vis_mode) {
            std::cerr << "ERROR: --sphere-flux formal anchors require one fixed "
                         "non-interactive /run/beamOn arm\n";
            delete run_manager;
            return 1;
        }
        if (!macro_declares_isotropic_gps(macro_file)) {
            std::cerr << "ERROR: --sphere-flux requires exactly one "
                         "/gps/ang/type iso declaration\n";
            delete run_manager;
            return 1;
        }
        double world_half_extent_m = 0.0;
        if (!get_world_half_extent_m_from_gdml(gdml_file, world_half_extent_m)) {
            std::cerr << "ERROR: --sphere-flux could not determine the world "
                         "half extent from GDML\n";
            delete run_manager;
            return 1;
        }
        double gdml_scoring_radius_cm = 0.0;
        if (!get_scoring_radius_cm_from_gdml(gdml_file, gdml_scoring_radius_cm)
            || !std::isfinite(sphere_flux_radius_cm)
            || std::abs(gdml_scoring_radius_cm - sphere_flux_radius_cm)
                   > 1.0e-10 * std::max(gdml_scoring_radius_cm,
                                        sphere_flux_radius_cm)) {
            std::cerr << "ERROR: --sphere-flux radius does not match the GDML "
                         "scoring sphere\n";
            delete run_manager;
            return 1;
        }
        CeeFitGeometryMetadata fit_geometry;
        if (!parse_ceefit_geometry_metadata(gdml_file, fit_geometry)) {
            std::cerr << "ERROR: --sphere-flux requires CEEFIT-GEOMETRY-v1 "
                         "metadata from env_scatter_export\n";
            delete run_manager;
            return 1;
        }
        G4ThreeVector source_position_mm;
        double source_energy_MeV = 0.0;
        int source_events = 0;
        if (!parse_macro(macro_file, source_position_mm, source_energy_MeV,
                         source_events)
            || source_events <= 0
            || !close_metadata(source_energy_MeV * 1000.0,
                               fit_geometry.energy_keV)
            || !close_metadata(source_position_mm.mag() / CLHEP::cm,
                               fit_geometry.distance_cm)
            || !close_metadata(source_position_mm.y() / CLHEP::cm,
                               fit_geometry.source_height_cm
                                   - fit_geometry.detector_height_cm)
            || !close_metadata(fit_geometry.scoring_radius_cm,
                               gdml_scoring_radius_cm)
            || !close_metadata(fit_geometry.world_half_extent_cm * 0.01,
                               world_half_extent_m)) {
            std::cerr << "ERROR: CEEFIT geometry metadata does not close to "
                         "the exact GDML/macro inputs\n";
            delete run_manager;
            return 1;
        }
        if (fit_geometry.geometry_tier != "air"
            && fit_geometry.geometry_tier != "field"
            && fit_geometry.geometry_tier != "room") {
            std::cerr << "ERROR: invalid CEEFIT geometry tier\n";
            delete run_manager;
            return 1;
        }
        if ((fit_geometry.geometry_tier == "air"
             && fit_geometry.material != "none")
            || (fit_geometry.geometry_tier != "air"
                && fit_geometry.material != "soil"
                && fit_geometry.material != "concrete")) {
            std::cerr << "ERROR: invalid CEEFIT geometry material\n";
            delete run_manager;
            return 1;
        }
        const std::string gdml_sha256 = file_sha256(gdml_file);
        const std::string macro_sha256 = file_sha256(macro_file);
        if (gdml_sha256.size() != 64 || macro_sha256.size() != 64) {
            std::cerr << "ERROR: could not hash the exact sphere-flux GDML "
                         "and macro inputs\n";
            delete run_manager;
            return 1;
        }
        run_action->EnableSphereFlux(true, sphere_flux_radius_cm);
        run_action->SetSphereRunMetadata(
            physics_list_name,
            world_half_extent_m,
            "gps-isotropic-4pi",
            G4Random::getTheEngine()->name(),
            random_seed,
            production_cut_mm,
            gdml_sha256,
            macro_sha256,
            fit_geometry.geometry_tier,
            fit_geometry.material,
            fit_geometry.distance_cm,
            fit_geometry.source_height_cm,
            fit_geometry.detector_height_cm,
            gdml_scoring_radius_cm);
        step_action->SetRunAction(run_action);
        event_action->SetSphereActions(step_action, run_action);
    }

    run_manager->SetUserAction(primary_gen);
    run_manager->SetUserAction(run_action);
    run_manager->SetUserAction(event_action);
    run_manager->SetUserAction(step_action);

    // Enable gamma-gamma angular correlation (G4PolarizationTransition) while
    // still in PreInit -- the deex parameter is locked once Initialize() runs.
    if (correlated_gamma) {
        G4NuclearLevelData::GetInstance()->GetParameters()->SetCorrelatedGamma(true);
        std::cout << "Correlated gamma emission ENABLED (W(theta) angular correlation)\n";
    }

    // Initialize.
    run_manager->Initialize();

    // Execute GPS macro (or in cone-bias mode, the macro's /run/beamOn triggers events).
    G4UImanager* ui = G4UImanager::GetUIpointer();

    if (vis_batch) {
        // Headless visualization: run the /vis/ macro and exit.  No
        // G4UIExecutive, so there is no window and no session to close --
        // usable over ssh and from scripts.  Use a file-writing driver
        // (VRML2FILE); it writes G4's own polyhedron for each solid, which is
        // what makes it a real check of the exported geometry rather than of
        // the GDML parser alone.  RayTracer is not usable here -- see the run
        // manager comment above.
        auto* vis_manager = new G4VisExecutive;
        vis_manager->Initialize();
        ui->ApplyCommand("/control/execute " + macro_file);
        delete vis_manager;
    } else if (vis_mode) {
        // Interactive visualization mode.
        // G4UIExecutive MUST be created first so Qt application context exists
        // before G4VisExecutive tries to register Qt-based OpenGL drivers.
        auto* ui_exec = new G4UIExecutive(argc, argv);

        auto* vis_manager = new G4VisExecutive;
        vis_manager->Initialize();

        // Execute the macro (which should contain /vis/open and drawing commands)
        ui->ApplyCommand("/control/execute " + macro_file);

        // Start interactive session (user can rotate, zoom, run commands)
        ui_exec->SessionStart();
        delete vis_manager;
        delete ui_exec;
    } else if (fep_precision_target > 0.0) {
        // Precision-driven mode: run the macro first to configure GPS,
        // then iterate with smaller beamOn calls until precision is met.
        // The macro's /run/beamOn sets up the initial batch.
        run_action->PreserveAcrossRuns();
        ui->ApplyCommand("/control/execute " + macro_file);

        constexpr uint64_t kSafetyMaxEvents = 100'000'000;
        constexpr int kBatchSize = 10000;

        while (!run_action->CheckFepPrecision(fep_precision_target) &&
               run_action->GetTotalEvents() < kSafetyMaxEvents) {
            ui->ApplyCommand("/run/beamOn " + std::to_string(kBatchSize));
        }

        std::cout << "FEP precision mode: stopped at "
                  << run_action->GetTotalEvents() << " events\n";
    } else {
        ui->ApplyCommand("/control/execute " + macro_file);
    }

    delete run_manager;
    return 0;
}
