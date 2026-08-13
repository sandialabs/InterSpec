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

#include "DetectorConstruction.hh"

#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>

DetectorConstruction::DetectorConstruction(const std::string& gdml_filename)
    : gdml_filename_(gdml_filename)
{}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    G4GDMLParser parser;
    parser.Read(gdml_filename_, /* validate = */ false);

    G4VPhysicalVolume* world = parser.GetWorldVolume();

    // Find the "active_crystal" logical volume for scoring (names are stripped
    // by G4GDMLParser::Read before returning).
    active_crystal_lv_ = parser.GetVolume("active_crystal");
    if (!active_crystal_lv_) {
        // Fallback: scan the store.
        G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
        for (auto* lv : *store) {
            if (lv->GetName() == "active_crystal") {
                active_crystal_lv_ = lv;
                break;
            }
        }
    }

    if (!active_crystal_lv_) {
        G4Exception("DetectorConstruction::Construct", "GDML001",
                    FatalException,
                    "Could not find logical volume named 'active_crystal' in GDML file.");
    }

    return world;
}
