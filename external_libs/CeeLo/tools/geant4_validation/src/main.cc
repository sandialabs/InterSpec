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
///   --cone-bias      Enable cone importance sampling (matches our MC's approach)
///                    Parses source position and energy from the macro file,
///                    detector geometry from the GDML.
///   --em-standard    Use G4EmStandardPhysics (default EM, no Doppler broadening)
///                    instead of option4
///   --diagnostics    Write per-interaction diagnostic histograms (11 CSVs)
///   --entry-diag     Write crystal entry diagnostic CSVs (energy, z-pos, summary)
///   --lowcut         Set production cut to 0.001 mm (default 0.7 mm) for thin
///                    detector studies where default cuts may suppress physics
///   --fep-precision <target>  Run until FEP relative precision <= target
///                    (e.g., 0.01 for 1%). Safety cap: 100M events.
///   --vis            Open interactive OpenGL visualization (no simulation).
///                    The macro file should contain /vis/ commands (e.g., vis.mac).
///
/// Physics list: G4EmStandardPhysics_option4 (most accurate EM) by default.

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include "FTFP_BERT.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

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
                double scale = (unit == "mm") ? 1.0 : 10.0; // cm -> mm
                source_pos_mm = G4ThreeVector(x * scale, y * scale, z * scale);
                got_pos = true;
            }
        }
        if (line.find("/gps/ene/mono") != std::string::npos) {
            std::istringstream iss(line);
            std::string cmd;
            double val;
            std::string unit;
            if (iss >> cmd >> val >> unit) {
                energy_MeV = (unit == "keV") ? val * 0.001 : val;
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
    bool diagnostics_mode = false;
    bool entry_diag_mode = false;
    bool low_cut = false;
    bool vis_mode = false;
    double bin_width_keV = 1.0;
    double fep_precision_target = 0.0; // 0 = disabled
    bool correlated_gamma = false; // --correlated-gamma: enable G4 gamma-gamma angular correlation
    for (int i = 4; i < argc; ++i) {
        if (std::string(argv[i]) == "--histogram") histogram_mode = true;
        if (std::string(argv[i]) == "--cone-bias") cone_bias = true;
        if (std::string(argv[i]) == "--em-standard") em_standard = true;
        if (std::string(argv[i]) == "--diagnostics") diagnostics_mode = true;
        if (std::string(argv[i]) == "--entry-diag") entry_diag_mode = true;
        if (std::string(argv[i]) == "--lowcut") low_cut = true;
        if (std::string(argv[i]) == "--vis") vis_mode = true;
        if (std::string(argv[i]) == "--correlated-gamma") correlated_gamma = true;
        if (std::string(argv[i]) == "--bin-width" && i + 1 < argc) {
            bin_width_keV = std::atof(argv[++i]);
        }
        if (std::string(argv[i]) == "--fep-precision" && i + 1 < argc) {
            fep_precision_target = std::atof(argv[++i]);
        }
    }

    // Build run manager (Serial mode).
    auto* run_manager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Serial);

    // Detector geometry from exported GDML.
    auto* det = new DetectorConstruction(gdml_file);
    run_manager->SetUserInitialization(det);

    // Physics list selection.
    auto* physics = new FTFP_BERT;
    if (em_standard) {
        physics->ReplacePhysics(new G4EmStandardPhysics());
        std::cout << "Using G4EmStandardPhysics (no Doppler broadening)\n";
    } else {
        physics->ReplacePhysics(new G4EmStandardPhysics_option4());
        std::cout << "Using G4EmStandardPhysics_option4 (full Doppler broadening)\n";
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

    // Enable histogram mode if requested.
    if (histogram_mode) {
        run_action->EnableHistogram(true);
        run_action->SetBinWidth(bin_width_keV);
        event_action->EnableHistogram(true);
        std::cout << "Histogram mode enabled (bin width: " << bin_width_keV << " keV)\n";
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

    if (vis_mode) {
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
