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

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Loads a GDML file exported by CeeLo and provides a pointer to the
/// "active_crystal" logical volume for SteppingAction scoring.
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit DetectorConstruction(const std::string& gdml_filename);
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    /// Returns the logical volume named "active_crystal" (set after Construct()).
    G4LogicalVolume* GetActiveCrystalLV() const { return active_crystal_lv_; }

private:
    std::string      gdml_filename_;
    G4LogicalVolume* active_crystal_lv_ = nullptr;
};
