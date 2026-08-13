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

#include "SteppingAction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4Gamma.hh"

SteppingAction::SteppingAction(EventAction* event_action,
                               G4LogicalVolume* active_crystal_lv)
    : event_action_(event_action)
    , active_crystal_lv_(active_crystal_lv)
{}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    // Lazy initialization: Construct() is called during Initialize() *after*
    // this SteppingAction is constructed, so the pointer passed at construction
    // time is always null.  Look it up from the store on the first step.
    if (!active_crystal_lv_) {
        active_crystal_lv_ =
            G4LogicalVolumeStore::GetInstance()->GetVolume("active_crystal");
        if (!active_crystal_lv_) return;
    }

    // Only score energy deposited in the active crystal logical volume.
    const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
    G4LogicalVolume* lv = touchable->GetVolume()->GetLogicalVolume();

    // --- Crystal entry diagnostic: detect gammas entering crystal from outside ---
    if (entry_diag_enabled_ && lv == active_crystal_lv_) {
        const G4Track* trk = step->GetTrack();
        G4StepPoint* pre_pt = step->GetPreStepPoint();
        if (pre_pt->GetStepStatus() == fGeomBoundary) {
            if (trk->GetDefinition() == G4Gamma::Definition()) {
                G4ThreeVector pos = pre_pt->GetPosition() / CLHEP::cm;
                double energy_keV = pre_pt->GetKineticEnergy() / CLHEP::keV;

                G4CrystalEntryRecord rec;
                rec.position = pos;
                rec.energy_keV = energy_keV;

                // Classify surface: face (z~0), side (r~R), back (z~L)
                double z = pos.z();
                if (z < 0.05)
                    rec.surface = G4CrystalEntryRecord::kFace;
                else if (z > crystal_length_ - 0.05)
                    rec.surface = G4CrystalEntryRecord::kBack;
                else
                    rec.surface = G4CrystalEntryRecord::kSide;

                event_action_->AddCrystalEntry(rec);
            } else {
                // Detect e-/e+ entering crystal from outside (boundary crossing)
                const G4String& pname = trk->GetDefinition()->GetParticleName();
                if (pname == "e-" || pname == "e+") {
                    event_action_->SetExternalElectronEntry();
                }
            }
        }
    }

    // --- Source-class (u/s) decomposition: tag the primary photon's first
    // physics interaction OUTSIDE the active crystal. For bare-detector
    // source-effect configs (11/12) every non-crystal volume is source
    // material/shielding, so this reproduces CeeLo's
    // SourceEffectDiag u/s classification. (For configs with detector
    // attenuators the tag also covers attenuator interactions.)
    if (lv != active_crystal_lv_) {
        const G4Track* trk = step->GetTrack();
        if (trk->GetTrackID() == 1 &&
            trk->GetDefinition() == G4Gamma::Definition() &&
            !event_action_->PrimaryReachedCrystal()) {
            const G4VProcess* proc =
                step->GetPostStepPoint()->GetProcessDefinedStep();
            if (proc) {
                const G4String& pn = proc->GetProcessName();
                if (pn == "compt" || pn == "phot" || pn == "Rayl" ||
                    pn == "conv") {
                    event_action_->SetPrimaryInteractedOutsideCrystal();
                }
            }
        }
        return;
    }

    // Primary gamma stepping inside the crystal: later interactions outside
    // (e.g. crystal-escape then backscatter off the source box) no longer
    // count toward the u/s source classification — CeeLo classifies on
    // the first source-geometry pass only.
    {
        const G4Track* trk = step->GetTrack();
        if (trk->GetTrackID() == 1 &&
            trk->GetDefinition() == G4Gamma::Definition())
            event_action_->SetPrimaryReachedCrystal();
    }

    // Geant4 native energy unit is MeV; no conversion needed.
    double edep = step->GetTotalEnergyDeposit(); // MeV
    if (edep > 0.0) {
        event_action_->AddEnergyDeposit(edep);

        // Tag deposit by particle type for diagnostic breakdown.
        const G4String& pname =
            step->GetTrack()->GetParticleDefinition()->GetParticleName();
        if (pname == "e-" || pname == "e+") {
            event_action_->AddElectronDeposit(edep);
        } else if (pname == "gamma") {
            event_action_->AddGammaDeposit(edep);
        }
        // Other particles (protons from nuclear interactions, etc.) are rare
        // and lumped into the total but not tagged separately.
    }

    // --- Per-interaction diagnostics ---
    if (!diagnostics_enabled_) return;

    const G4Track* track = step->GetTrack();

    // Track secondary gammas (fluorescence X-rays, etc.) escaping the crystal.
    if (track->GetDefinition() == G4Gamma::Definition() && track->GetTrackID() > 1) {
        G4StepPoint* pre_pt = step->GetPreStepPoint();
        G4StepPoint* post_pt = step->GetPostStepPoint();
        G4LogicalVolume* pre_lv_sec = pre_pt->GetTouchable()->GetVolume()
                                          ->GetLogicalVolume();
        G4LogicalVolume* post_lv_sec = nullptr;
        if (post_pt->GetTouchable()->GetVolume())
            post_lv_sec = post_pt->GetTouchable()->GetVolume()->GetLogicalVolume();

        if (pre_lv_sec == active_crystal_lv_ && post_lv_sec != active_crystal_lv_) {
            double escape_keV = post_pt->GetKineticEnergy() * 1000.0; // MeV -> keV
            event_action_->AddFluorEscape(escape_keV);
        }
    }

    // Only record interaction details for gamma tracks (not electron/positron steps)
    if (track->GetDefinition() != G4Gamma::Definition()) return;

    // Only primary gamma (trackID == 1)
    if (track->GetTrackID() != 1) return;

    // Get the process that defined this step
    const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (!proc) return;

    const G4String& procName = proc->GetProcessName();

    G4InteractionRecord rec;
    rec.position = step->GetPostStepPoint()->GetPosition();
    rec.energy_before_MeV = step->GetPreStepPoint()->GetKineticEnergy();
    rec.energy_after_MeV = step->GetPostStepPoint()->GetKineticEnergy();

    // Scattering angle from momentum directions
    G4ThreeVector d_pre = step->GetPreStepPoint()->GetMomentumDirection();
    G4ThreeVector d_post = step->GetPostStepPoint()->GetMomentumDirection();
    rec.cos_theta = d_pre.dot(d_post);

    if (procName == "compt") {
        rec.type = G4InteractionRecord::kCompton;
    } else if (procName == "phot") {
        rec.type = G4InteractionRecord::kPhotoelectric;
        rec.energy_after_MeV = 0.0;
        rec.cos_theta = 1.0;
    } else if (procName == "Rayl") {
        rec.type = G4InteractionRecord::kRayleigh;
    } else if (procName == "conv") {
        rec.type = G4InteractionRecord::kPairProduction;
        rec.energy_after_MeV = 0.0;
        rec.cos_theta = 1.0;
    } else {
        return; // Transportation, etc. — not a physics interaction
    }

    // compton_index will be assigned by EventAction::AddInteraction()
    rec.compton_index = 0;
    event_action_->AddInteraction(rec);
}
