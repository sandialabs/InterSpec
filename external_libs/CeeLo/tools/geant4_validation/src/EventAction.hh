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

#include "G4UserEventAction.hh"
#include "SteppingAction.hh"  // for G4InteractionRecord

#include <cstdint>
#include <vector>

class RunAction;
class PrimaryGeneratorAction;

/// Accumulates energy deposited in the active crystal per event.
/// At end-of-event, classifies the event as FEP or any-detection and
/// notifies RunAction.
///
/// When histogram mode is enabled, also tracks per-particle-type deposits
/// and fills a fixed-bin energy histogram.
///
/// When diagnostics mode is enabled, collects per-interaction records
/// and forwards them to RunAction for histogram filling.
class EventAction : public G4UserEventAction {
public:
    EventAction(RunAction* run_action,
                PrimaryGeneratorAction* primary_gen = nullptr);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    /// Called by SteppingAction for each step inside the active crystal.
    void AddEnergyDeposit(double edep_MeV) { energy_deposited_MeV_ += edep_MeV; }

    /// Called by SteppingAction: energy deposited by e-/e+.
    void AddElectronDeposit(double edep_MeV) { edep_electron_MeV_ += edep_MeV; }

    /// Called by SteppingAction: energy deposited by gammas.
    void AddGammaDeposit(double edep_MeV) { edep_gamma_MeV_ += edep_MeV; }

    /// Called by SteppingAction: add interaction record for primary gamma.
    /// Assigns compton_index by counting prior Compton records.
    void AddInteraction(G4InteractionRecord rec);

    /// Called by SteppingAction: record a secondary gamma escaping the crystal.
    void AddFluorEscape(double escape_energy_keV) {
        n_fluor_escaped_++;
        sum_fluor_escape_keV_ += escape_energy_keV;
        fluor_escape_energies_.push_back(escape_energy_keV);
    }

    /// Called by SteppingAction: record a photon entering the crystal.
    void AddCrystalEntry(const G4CrystalEntryRecord& rec) {
        crystal_entries_.push_back(rec);
    }

    /// Called by SteppingAction: record an e-/e+ entering crystal from outside.
    void SetExternalElectronEntry() { has_external_electron_ = true; }

    /// Called by SteppingAction: the primary gamma had a physics interaction
    /// outside the active crystal before first reaching it (u/s source-class
    /// decomposition; mirrors CeeLo's SourceEffectDiag).
    void SetPrimaryInteractedOutsideCrystal() {
        primary_interacted_outside_ = true;
    }

    /// Called by SteppingAction when the primary gamma steps inside the
    /// crystal; freezes the u/s classification at its first-pass value.
    void SetPrimaryReachedCrystal() { primary_reached_crystal_ = true; }
    bool PrimaryReachedCrystal() const { return primary_reached_crystal_; }

    /// Primary photon energy set at the start of each event (MeV).
    double primary_energy_MeV = 0.0;

    /// Enable diagnostic histogram mode.
    void EnableHistogram(bool on) { histogram_enabled_ = on; }

    /// Enable per-interaction diagnostics.
    void EnableDiagnostics(bool on) { diagnostics_enabled_ = on; }

    /// Enable crystal entry diagnostics.
    void EnableEntryDiag(bool on) { entry_diag_enabled_ = on; }

private:
    RunAction* run_action_;
    PrimaryGeneratorAction* primary_gen_;
    double     energy_deposited_MeV_ = 0.0;
    double     edep_electron_MeV_ = 0.0;
    double     edep_gamma_MeV_ = 0.0;
    bool       histogram_enabled_ = false;
    bool       diagnostics_enabled_ = false;
    bool       entry_diag_enabled_ = false;

    // Per-event interaction records (cleared each event)
    std::vector<G4InteractionRecord> interactions_;
    int compton_count_ = 0;

    // Per-event fluorescence escape tracking
    int n_fluor_escaped_ = 0;
    double sum_fluor_escape_keV_ = 0.0;
    std::vector<double> fluor_escape_energies_;

    // Per-event crystal entry records
    std::vector<G4CrystalEntryRecord> crystal_entries_;

    // Per-event external electron entry flag
    bool has_external_electron_ = false;

    // Per-event u/s source-class flags
    bool primary_interacted_outside_ = false;
    bool primary_reached_crystal_ = false;
};
