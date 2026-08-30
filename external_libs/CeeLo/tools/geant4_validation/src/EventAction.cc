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

#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"

EventAction::EventAction(RunAction* run_action,
                         PrimaryGeneratorAction* primary_gen)
    : run_action_(run_action)
    , primary_gen_(primary_gen)
{}

void EventAction::BeginOfEventAction(const G4Event* event) {
    energy_deposited_MeV_ = 0.0;
    edep_electron_MeV_ = 0.0;
    edep_gamma_MeV_ = 0.0;
    n_fluor_escaped_ = 0;
    sum_fluor_escape_keV_ = 0.0;
    fluor_escape_energies_.clear();

    if (diagnostics_enabled_) {
        interactions_.clear();
        compton_count_ = 0;
    }

    if (entry_diag_enabled_) {
        crystal_entries_.clear();
    }

    has_external_electron_ = false;
    primary_interacted_outside_ = false;
    primary_reached_crystal_ = false;

    // Record the primary photon energy for FEP classification.
    if (event->GetNumberOfPrimaryVertex() > 0) {
        const G4PrimaryVertex* vtx = event->GetPrimaryVertex(0);
        if (vtx->GetNumberOfParticle() > 0) {
            // GetKineticEnergy() returns MeV (Geant4 native units); no conversion needed.
            primary_energy_MeV = vtx->GetPrimary(0)->GetKineticEnergy();
        }
    }
}

void EventAction::AddInteraction(G4InteractionRecord rec) {
    if (rec.type == G4InteractionRecord::kCompton) {
        rec.compton_index = compton_count_++;
    } else if (rec.type == G4InteractionRecord::kPhotoelectric) {
        rec.compton_index = compton_count_;  // # Compton scatters before this PE
    }
    interactions_.push_back(rec);
}

void EventAction::EndOfEventAction(const G4Event*) {
    // Half-width of the full-energy-peak window, in MeV.  MUST match the window
    //  CeeLo scores with (ceelo::kDefaultFepWindowKeV, physics/FepWindow.h) or
    //  the comparison is apples-to-oranges: a reference generated at a wider
    //  window counts events CeeLo does not.  Settable so a reference set can be
    //  regenerated at whatever window CeeLo currently uses -- and so the value
    //  can be recorded alongside the data it produced.
    const double kFepTol = fep_window_keV_ * 1.0e-3;

    // Get importance weight from primary generator (1.0 if no cone bias)
    double weight = 1.0;
    if (primary_gen_) {
        weight = primary_gen_->GetEventWeight();
    }

    const bool any_deposit = energy_deposited_MeV_ > 0.0;
    const bool is_fep = any_deposit &&
        std::abs(energy_deposited_MeV_ - primary_energy_MeV) < kFepTol;

    if (any_deposit) {
        run_action_->RecordAny(energy_deposited_MeV_, weight);
        if (is_fep) {
            run_action_->RecordFEP(energy_deposited_MeV_, weight);
        }
    }

    // u/s source-class decomposition (mirrors CeeLo SourceEffectDiag).
    run_action_->RecordSrcClass(primary_interacted_outside_, any_deposit,
                                is_fep, weight);

    // Always add to weight sum (for proper normalization)
    run_action_->RecordWeight(weight);

    // Histogram mode: record per-event energy deposit + per-particle breakdown.
    if (histogram_enabled_) {
        run_action_->RecordHistogram(energy_deposited_MeV_,
                                     edep_electron_MeV_,
                                     edep_gamma_MeV_,
                                     weight);
    }

    // Diagnostics mode: forward interaction records to RunAction.
    if (diagnostics_enabled_) {
        run_action_->RecordDiagnostics(interactions_,
                                        energy_deposited_MeV_,
                                        primary_energy_MeV);

        // Forward fluorescence escape data
        for (double e_keV : fluor_escape_energies_) {
            run_action_->RecordFluorEscape(e_keV);
        }
    }

    // Crystal entry diagnostic mode.
    if (entry_diag_enabled_ && !crystal_entries_.empty()) {
        run_action_->RecordCrystalEntries(crystal_entries_);
    }

    // External electron entry tracking: report to RunAction if any e-/e+
    // crossed into the crystal from outside this event.
    if (has_external_electron_) {
        bool gamma_entered = !crystal_entries_.empty();
        run_action_->RecordExternalElectron(energy_deposited_MeV_ > 0.0,
                                             gamma_entered);
    }
}
