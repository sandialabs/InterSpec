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

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"

#include <cstdint>

class EventAction;
class G4LogicalVolume;

/// Interaction record for per-interaction diagnostic histograms.
struct G4InteractionRecord {
    enum Type : uint8_t { kCompton, kPhotoelectric, kRayleigh, kPairProduction };
    Type type;
    G4ThreeVector position;
    double energy_before_MeV;  ///< Photon KE before this interaction
    double energy_after_MeV;   ///< Photon KE after (0 for PE)
    double cos_theta;          ///< Cosine of scattering angle
    int compton_index;         ///< 0-based Compton scatter index (or # prior Comptons for PE)
};

/// Record of a photon entering the crystal from outside.
struct G4CrystalEntryRecord {
    G4ThreeVector position;      ///< Entry position (cm)
    double energy_keV;           ///< Photon kinetic energy at entry
    enum Surface : uint8_t { kFace, kSide, kBack };
    Surface surface;
};

/// Accumulates energy deposits in the "active_crystal" volume only.
/// Steps in air (bore hole) or attenuator volumes are ignored.
///
/// When diagnostics are enabled, also records per-interaction data
/// for the primary gamma in the active crystal.
class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction(EventAction* event_action, G4LogicalVolume* active_crystal_lv);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

    void EnableDiagnostics(bool on) { diagnostics_enabled_ = on; }
    void EnableEntryDiag(bool on) { entry_diag_enabled_ = on; }

    /// Set crystal dimensions for entry surface classification.
    void SetCrystalDimensions(double length_cm, double radius_cm) {
        crystal_length_ = length_cm;
        crystal_radius_ = radius_cm;
    }

private:
    EventAction*     event_action_;
    G4LogicalVolume* active_crystal_lv_;
    bool diagnostics_enabled_ = false;
    bool entry_diag_enabled_ = false;
    double crystal_length_ = 7.62;  // cm
    double crystal_radius_ = 3.81;  // cm
};
