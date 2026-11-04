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

// Geant4 anchor driver for a centered point source in a homogeneous spherical
// elemental shield. Output matches the CeeLo source campaign CSV schema.

#include "cross_sections/CrossSectionData.h"
#include "materials/Material.h"

#include "FTFP_BERT.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessTable.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VProcess.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "Randomize.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

enum Component : std::size_t {
    Direct,
    Elastic,
    SingleCompton,
    DoubleCompton,
    TripleCompton,
    QuadrupleCompton,
    FivePlusCompton,
    Pair511,
    PairAnnihilationContinuum,
    PairAnnihilationInFlight,
    PhotoelectronBrems,
    ComptonElectronBrems,
    PairElectronBrems,
    PhotonuclearSecondary,
    FluorescenceKAlpha,
    FluorescenceKBeta,
    FluorescenceL,
    FluorescenceContinuum,
    OtherSecondaryPhoton,
    Absorbed,
    NumComponents
};

constexpr const char* kComponentNames[] = {
    "uncollided", "elastic_line", "single_compton", "double_compton",
    "triple_compton", "quadruple_compton", "five_plus_compton",
    "pair_511", "pair_annihilation_continuum",
    "pair_annihilation_in_flight",
    "photoelectron_brems", "compton_electron_brems",
    "pair_electron_brems", "photonuclear_secondary",
    "fluorescence_k_alpha", "fluorescence_k_beta", "fluorescence_l",
    "fluorescence_continuum",
    "other_secondary_photon", "absorbed"};

enum class TrackOrigin : std::uint8_t {
    Primary,
    Photoelectron,
    ComptonElectron,
    PairElectron,
    PairAnnihilation,
    PairAnnihilationInFlight,
    Photonuclear,
    Fluorescence,
    Other
};

struct Config {
    double energy_keV = 661.657;
    int Z = 82;
    double tau = 1.0;
    double areal_density_g_cm2 = 0.0;
    bool areal_density_was_set = false;
    // Shield hydrogen mass fraction; the remainder is element Z. Matches the
    // CeeLo campaign's --hydrogen-mass-fraction two-component convention.
    double hydrogen_mass_fraction = 0.0;
    std::uint64_t events = 100000;
    long seed = 246813579;
    int bins = 256;
    bool option4 = false;
    bool source_rayleigh = true;
    // Atomic deexcitation (K and L x-rays, Auger) and, with pixe, electron-
    // impact ionization x-rays at every energy -- the CeeLo campaign's
    // source fluorescence and impact ionization.
    bool fluorescence = false;
    bool pixe = false;
    std::string output = "g4_source_escape.csv";
};

struct Tally {
    double sum_w = 0.0;
    double sum_w2 = 0.0;
    void add(double weight = 1.0) {
        sum_w += weight;
        sum_w2 += weight * weight;
    }
};

class Detector final : public G4VUserDetectorConstruction {
public:
    explicit Detector(const Config& config) : config_(config) {}

    G4VPhysicalVolume* Construct() override {
        auto* nist = G4NistManager::Instance();
        auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
        auto* element = nist->FindOrBuildElement(config_.Z);
        const double f_H = config_.hydrogen_mass_fraction;
        auto* shield = new G4Material(
            "EquivalentElement", 1.0 * g / cm3, f_H > 0.0 ? 2 : 1);
        if (f_H > 0.0)
            shield->AddElement(nist->FindOrBuildElement(1), f_H);
        shield->AddElement(element, 1.0 - f_H);

        const double radius =
            std::max(1e-7, config_.areal_density_g_cm2) * cm;
        const double world_half = radius + 20.0 * cm;
        auto* world_solid = new G4Box(
            "WorldSolid", world_half, world_half, world_half);
        auto* world_logic = new G4LogicalVolume(
            world_solid, vacuum, "World");
        auto* world = new G4PVPlacement(
            nullptr, {}, world_logic, "World", nullptr, false, 0);
        auto* shield_solid = new G4Orb("ShieldSolid", radius);
        auto* shield_logic = new G4LogicalVolume(
            shield_solid, shield, "Shield");
        new G4PVPlacement(
            nullptr, {}, shield_logic, "Shield", world_logic, false, 0);
        return world;
    }

private:
    Config config_;
};

class Primary final : public G4VUserPrimaryGeneratorAction {
public:
    explicit Primary(double energy_keV)
        : gun_(1), energy_(energy_keV * keV) {
        gun_.SetParticleDefinition(G4Gamma::GammaDefinition());
        gun_.SetParticlePosition({});
        gun_.SetParticleEnergy(energy_);
    }

    void GeneratePrimaries(G4Event* event) override {
        const double cosine = 2.0 * G4UniformRand() - 1.0;
        const double phi = twopi * G4UniformRand();
        const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
        gun_.SetParticleMomentumDirection(
            {sine * std::cos(phi), sine * std::sin(phi), cosine});
        gun_.GeneratePrimaryVertex(event);
    }

private:
    G4ParticleGun gun_;
    double energy_;
};

struct EventState {
    int compton = 0;
    int rayleigh = 0;
    int pair = 0;
    bool primary_escaped = false;
};

class Run;

class Event final : public G4UserEventAction {
public:
    explicit Event(Run* run) : run_(run) {}
    void BeginOfEventAction(const G4Event*) override {
        state_ = {};
        origins_.clear();
    }
    void EndOfEventAction(const G4Event*) override;
    EventState& state() { return state_; }
    void register_track(const G4Track& track);
    void escape(const G4Track& track);

private:
    Run* run_;
    EventState state_;
    std::unordered_map<int, TrackOrigin> origins_;
};

class Run final : public G4UserRunAction {
public:
    explicit Run(const Config& config) : config_(config) {
        for (auto& values : spectra_)
            values.resize(static_cast<std::size_t>(config.bins));
        // The shield element's lines (CeeLo's tables) and their families.
        const auto& xs = ceelo::CrossSectionData::instance();
        if (const ceelo::FluorescenceData* k = xs.fluorescence(config.Z))
            for (int i = 0; i < k->num_lines; ++i)
                lines_.emplace_back(
                    k->line_energy_keV[i],
                    ceelo::k_line_family(k->line_filling_subshell[i])
                            == ceelo::KLineFamily::Alpha
                        ? FluorescenceKAlpha : FluorescenceKBeta);
        if (const ceelo::LFluorescenceData* l = xs.l_fluorescence(config.Z))
            for (const auto& sub : l->sub)
                for (int i = 0; i < sub.num_lines; ++i)
                    lines_.emplace_back(sub.line_energy_keV[i], FluorescenceL);
    }

    /// An escaping x-ray within 0.5 keV of one of the shield's lines (GEANT4's
    /// line energies differ from CeeLo's by a few eV) is that line's family's
    /// atom; anything else (a Compton-scattered x-ray) is continuum.
    Component fluorescence_family(double energy_keV) const {
        double best = 0.5;
        Component family = FluorescenceContinuum;
        for (const auto& [line, f] : lines_)
            if (std::abs(line - energy_keV) < best) {
                best = std::abs(line - energy_keV);
                family = f;
            }
        return family;
    }

    void tally(Component component, double energy_keV, bool atom) {
        areas_[component].add();
        if (!atom && energy_keV >= 0.0 && energy_keV < config_.energy_keV) {
            const int bin = std::min(
                config_.bins - 1,
                static_cast<int>(
                    energy_keV / config_.energy_keV * config_.bins));
            spectra_[component][static_cast<std::size_t>(bin)].add();
        }
    }

    void EndOfRunAction(const G4Run*) override {
        std::ofstream stream(config_.output);
        if (!stream)
            throw std::runtime_error("cannot open " + config_.output);
        stream << std::setprecision(12)
            << "cell_kind,energy_keV,Z,tau,areal_density_g_cm2,cos_theta,"
               "detector,component,bin_lo_keV,bin_hi_keV,sum_w,sum_w2,"
               "events,seed,analytic_direct,rayleigh_enabled,"
               "compton_doppler_enabled,fluorescence_enabled,physics_list,"
               "hydrogen_mass_fraction\n";
        for (std::size_t c = 0; c < NumComponents; ++c) {
            const auto component = static_cast<Component>(c);
            const bool fluorescence_atom = component == FluorescenceKAlpha
                || component == FluorescenceKBeta || component == FluorescenceL;
            const bool line = component == Direct || component == Elastic
                || component == Pair511 || component == Absorbed
                || fluorescence_atom;
            if (line) {
                const double energy = component == Pair511 ? 510.99895
                    : (component == Absorbed ? -1.0
                       : (fluorescence_atom ? 0.0 : config_.energy_keV));
                stream << "source," << config_.energy_keV << ','
                       << config_.Z << ',' << config_.tau << ','
                       << config_.areal_density_g_cm2 << ",,,"
                       << kComponentNames[c] << ',' << energy << ',' << energy
                       << ',' << areas_[c].sum_w << ',' << areas_[c].sum_w2
                       << ',' << config_.events << ',' << config_.seed << ','
                       << std::exp(-config_.tau) << ','
                       << config_.source_rayleigh << ",1," << config_.fluorescence << ','
                       << (config_.option4 ? "option4" : "livermore")
                       << ',' << config_.hydrogen_mass_fraction
                       << '\n';
            } else {
                for (int bin = 0; bin < config_.bins; ++bin) {
                    const double lo =
                        config_.energy_keV * bin / config_.bins;
                    const double hi =
                        config_.energy_keV * (bin + 1) / config_.bins;
                    const Tally& tally =
                        spectra_[c][static_cast<std::size_t>(bin)];
                    stream << "source," << config_.energy_keV << ','
                           << config_.Z << ',' << config_.tau << ','
                           << config_.areal_density_g_cm2 << ",,,"
                           << kComponentNames[c] << ',' << lo << ',' << hi
                           << ',' << tally.sum_w << ',' << tally.sum_w2
                           << ',' << config_.events << ',' << config_.seed
                           << ',' << std::exp(-config_.tau) << ','
                           << config_.source_rayleigh << ",1," << config_.fluorescence << ','
                           << (config_.option4 ? "option4" : "livermore")
                           << ',' << config_.hydrogen_mass_fraction
                           << '\n';
                }
            }
        }
    }

private:
    Config config_;
    std::vector<std::pair<double, Component>> lines_;
    std::array<Tally, NumComponents> areas_;
    std::array<std::vector<Tally>, NumComponents> spectra_;
};

void Event::EndOfEventAction(const G4Event*) {
    if (!state_.primary_escaped)
        run_->tally(Absorbed, -1.0, true);
}

void Event::register_track(const G4Track& track) {
    const int track_id = track.GetTrackID();
    if (track_id == 1) {
        origins_[track_id] = TrackOrigin::Primary;
        return;
    }
    const auto parent = origins_.find(track.GetParentID());
    TrackOrigin origin =
        parent == origins_.end() ? TrackOrigin::Other : parent->second;
    const G4VProcess* creator = track.GetCreatorProcess();
    const std::string name = creator ? creator->GetProcessName() : "";
    if (name.find("photonNuclear") != std::string::npos) {
        origin = TrackOrigin::Photonuclear;
    } else if (track.GetDefinition() == G4Gamma::GammaDefinition()) {
        if (name.find("annihil") != std::string::npos) {
            origin = std::abs(
                track.GetKineticEnergy() / keV - 510.99895) < 0.05
                ? TrackOrigin::PairAnnihilation
                : TrackOrigin::PairAnnihilationInFlight;
        } else if (name == "phot" || name == "compt" || name == "eIoni"
                   || name == "conv" || name == "hIoni") {
            // Atomic deexcitation after a vacancy (photoabsorption, Compton,
            // electron-impact ionization).
            origin = TrackOrigin::Fluorescence;
        }
        // Brems photons inherit the physical origin assigned to the parent
        // electron or positron.
    } else {
        if (name == "phot")
            origin = TrackOrigin::Photoelectron;
        else if (name == "compt")
            origin = TrackOrigin::ComptonElectron;
        else if (name == "conv")
            origin = TrackOrigin::PairElectron;
    }
    origins_[track_id] = origin;
}

void Event::escape(const G4Track& track) {
    const double energy = track.GetKineticEnergy() / keV;
    if (track.GetTrackID() == 1) {
        state_.primary_escaped = true;
        Component component = Direct;
        if (state_.compton == 1)
            component = SingleCompton;
        else if (state_.compton == 2)
            component = DoubleCompton;
        else if (state_.compton == 3)
            component = TripleCompton;
        else if (state_.compton == 4)
            component = QuadrupleCompton;
        else if (state_.compton >= 5)
            component = FivePlusCompton;
        else if (state_.rayleigh > 0)
            component = Elastic;
        run_->tally(
            component, energy, component == Direct || component == Elastic);
        return;
    }
    const auto found = origins_.find(track.GetTrackID());
    const TrackOrigin origin =
        found == origins_.end() ? TrackOrigin::Other : found->second;
    Component component = OtherSecondaryPhoton;
    bool atom = false;
    switch (origin) {
    case TrackOrigin::PairAnnihilation:
        atom = std::abs(energy - 510.99895) < 0.05;
        component = atom ? Pair511 : PairAnnihilationContinuum;
        break;
    case TrackOrigin::PairAnnihilationInFlight:
        component = PairAnnihilationInFlight;
        break;
    case TrackOrigin::Photoelectron:
        component = PhotoelectronBrems;
        break;
    case TrackOrigin::ComptonElectron:
        component = ComptonElectronBrems;
        break;
    case TrackOrigin::PairElectron:
        component = PairElectronBrems;
        break;
    case TrackOrigin::Photonuclear:
        component = PhotonuclearSecondary;
        break;
    case TrackOrigin::Fluorescence:
        // Lines below 10 keV are outside the CeeLo model and campaign.
        if (energy < 10.0)
            return;
        component = run_->fluorescence_family(energy);
        atom = component != FluorescenceContinuum;
        break;
    default:
        break;
    }
    run_->tally(component, energy, atom);
}

class Tracking final : public G4UserTrackingAction {
public:
    explicit Tracking(Event* event) : event_(event) {}

    void PreUserTrackingAction(const G4Track* track) override {
        event_->register_track(*track);
    }

private:
    Event* event_;
};

class Stepping final : public G4UserSteppingAction {
public:
    explicit Stepping(Event* event) : event_(event) {}

    void UserSteppingAction(const G4Step* step) override {
        const G4Track* track = step->GetTrack();
        const auto* pre = step->GetPreStepPoint()->GetPhysicalVolume();
        const auto* post = step->GetPostStepPoint()->GetPhysicalVolume();
        if (!pre || pre->GetName() != "Shield")
            return;
        if (track->GetTrackID() == 1) {
            const auto* process =
                step->GetPostStepPoint()->GetProcessDefinedStep();
            const std::string name = process ? process->GetProcessName() : "";
            if (name == "compt")
                ++event_->state().compton;
            else if (name == "Rayl")
                ++event_->state().rayleigh;
            else if (name == "conv")
                ++event_->state().pair;
        }
        if (post && post->GetName() == "World"
            && track->GetDefinition() == G4Gamma::GammaDefinition()) {
            event_->escape(*track);
            const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
        }
    }

private:
    Event* event_;
};

Config parse(int argc, char** argv) {
    Config config;
    auto value = [&](int& index) -> std::string {
        if (++index >= argc)
            throw std::invalid_argument("missing option value");
        return argv[index];
    };
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--energy-kev")
            config.energy_keV = std::stod(value(i));
        else if (arg == "--z")
            config.Z = std::stoi(value(i));
        else if (arg == "--tau")
            config.tau = std::stod(value(i));
        else if (arg == "--areal-density") {
            config.areal_density_g_cm2 = std::stod(value(i));
            config.areal_density_was_set = true;
        }
        else if (arg == "--hydrogen-mass-fraction")
            config.hydrogen_mass_fraction = std::stod(value(i));
        else if (arg == "--events")
            config.events = std::stoull(value(i));
        else if (arg == "--seed")
            config.seed = std::stol(value(i));
        else if (arg == "--bins")
            config.bins = std::stoi(value(i));
        else if (arg == "--output")
            config.output = value(i);
        else if (arg == "--em-option4")
            config.option4 = true;
        else if (arg == "--source-rayleigh")
            config.source_rayleigh = true;
        else if (arg == "--no-source-rayleigh")
            config.source_rayleigh = false;
        else if (arg == "--fluorescence")
            config.fluorescence = true;
        else if (arg == "--pixe")
            config.pixe = true;
        else
            throw std::invalid_argument("unknown option: " + arg);
    }
    if (config.energy_keV < 10.0 || config.energy_keV > 12000.0
        || config.Z < 1 || config.Z > 92 || config.tau < 0.0)
        throw std::invalid_argument("invalid energy, Z, or tau");
    if (config.hydrogen_mass_fraction < 0.0
        || config.hydrogen_mass_fraction >= 1.0)
        throw std::invalid_argument(
            "hydrogen mass fraction must be in [0,1)");
    std::vector<ceelo::MaterialComponent> composition;
    if (config.hydrogen_mass_fraction > 0.0)
        composition.push_back(
            {std::uint8_t{1}, config.hydrogen_mass_fraction});
    composition.push_back(
        {static_cast<std::uint8_t>(config.Z),
         1.0 - config.hydrogen_mass_fraction});
    ceelo::Material material("equivalent", 1.0, composition);
    const auto xs =
        material.macroscopic_xs(config.energy_keV * 1e-3);
    const double mu_mass = (
        xs.mu_pe + xs.mu_cs + xs.mu_pp
        + (config.source_rayleigh ? xs.mu_rs : 0.0))
        / material.density();
    if (config.areal_density_was_set) {
        if (config.areal_density_g_cm2 < 0.0)
            throw std::invalid_argument(
                "areal density must be nonnegative");
        config.tau = config.areal_density_g_cm2 * mu_mass;
    } else {
        config.areal_density_g_cm2 =
            config.tau == 0.0 ? 0.0 : config.tau / mu_mass;
    }
    return config;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse(argc, argv);
        G4Random::setTheSeed(config.seed);
        auto* manager = G4RunManagerFactory::CreateRunManager(
            G4RunManagerType::Serial);
        manager->SetUserInitialization(new Detector(config));
        auto* physics = new FTFP_BERT(0);
        physics->SetVerboseLevel(0);
        if (config.option4)
            physics->ReplacePhysics(new G4EmStandardPhysics_option4());
        else
            physics->ReplacePhysics(new G4EmLivermorePhysics());
        auto* em_parameters = G4EmParameters::Instance();
        em_parameters->SetVerbose(0);
        em_parameters->SetFluo(config.fluorescence);
        em_parameters->SetAuger(config.fluorescence);
        em_parameters->SetDeexcitationIgnoreCut(config.fluorescence);
        em_parameters->SetPixe(config.pixe);
        // Keep gamma processes separate so the coherent process can be
        // deactivated explicitly for the far-field source campaign.
        em_parameters->SetGeneralProcessActive(false);
        physics->SetDefaultCutValue(0.001 * mm);
        manager->SetUserInitialization(physics);
        auto* run = new Run(config);
        auto* event = new Event(run);
        manager->SetUserAction(new Primary(config.energy_keV));
        manager->SetUserAction(run);
        manager->SetUserAction(event);
        manager->SetUserAction(new Tracking(event));
        manager->SetUserAction(new Stepping(event));
        manager->Initialize();
        if (!config.source_rayleigh) {
            G4ProcessTable::GetProcessTable()->SetProcessActivation(
                "Rayl", G4Gamma::GammaDefinition(), false);
        }
        manager->BeamOn(static_cast<G4int>(config.events));
        delete manager;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "g4_source_escape_campaign: " << error.what() << '\n';
        return 2;
    }
}
