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

/// @file Geant4Export.h
/// @brief GDML geometry export and GEANT4 GPS macro generation.
///
/// Provides two free functions:
///   - write_gdml(): exports the detector geometry as a GDML file that
///     GEANT4 can directly import via G4GDMLParser.
///   - write_geant4_macro(): writes a GPS macro file that sets up a point
///     source for the GEANT4 validation harness.
///
/// The exported GDML names the active detector volume "active_crystal" so
/// that the GEANT4 SteppingAction can identify it for energy tallying.
///
/// Dead layers are not explicitly modeled in the GDML (they are noted in a
/// comment); the full crystal volume is exported as "active_crystal".
/// Attenuator shells are exported as separate physical volumes.

#include <string>
#include <cstdint>

#include <Eigen/Core>

namespace ceelo {

class Geometry;
class SourceGeometry;

/// Write a GDML geometry file describing the detector geometry.
///
/// The file can be loaded into GEANT4 with:
///   G4GDMLParser parser;
///   parser.Read(filename);
///   G4VPhysicalVolume* world = parser.GetWorldVolume();
///
/// @param geom      The detector geometry (crystal, bore, dead layers, attenuators).
/// @param filename  Output file path (e.g., "detector.gdml").
/// @param source_geom  Optional source geometry (material + shielding shells).
///                      If non-null, source material and shielding volumes are
///                      included in the GDML using nested placement.
/// @param vacuum_world  If true, the world volume uses G4_Galactic (vacuum)
///                      instead of Air. Use this when MC does not include air
///                      attenuation, so G4 and MC are consistent.
void write_gdml(const Geometry& geom, const std::string& filename,
                const SourceGeometry* source_geom = nullptr,
                bool vacuum_world = false);

/// Write a GEANT4 GPS macro file for a point-source monoenergetic simulation.
///
/// The macro sets up:
///   - Particle: gamma
///   - Position: the given source position (detector frame, cm)
///   - Energy: monoenergetic at energy_keV
///   - Run: num_events events
///
/// @param source_pos   Source position in the detector frame (cm).
/// @param energy_keV   Photon energy in keV.
/// @param num_events   Number of events to simulate.
/// @param filename     Output file path (e.g., "run.mac").
void write_geant4_macro(const Eigen::Vector3d& source_pos,
                        double energy_keV,
                        uint64_t num_events,
                        const std::string& filename);

/// Write a GEANT4 GPS macro for a volumetric Marinelli beaker source.
///
/// Uses /gps/pos/type Volume with /gps/pos/confine to restrict source
/// points to the SrcMaterialPV physical volume (the L-shaped sample).
/// GEANT4 samples uniformly in the bounding cylinder and rejects points
/// outside the logical volume.
///
/// @param source_geom  The Marinelli source geometry (must be Shape::Marinelli).
/// @param energy_keV   Photon energy in keV.
/// @param num_events   Number of events to simulate.
/// @param filename     Output file path (e.g., "run_marinelli.mac").
void write_geant4_macro_marinelli(const SourceGeometry& source_geom,
                                   double energy_keV,
                                   uint64_t num_events,
                                   const std::string& filename);

/// Write a GEANT4 GPS macro for a cylindrical volume source.
/// Uses GPS /gps/pos/type Volume with /gps/pos/shape Cylinder.
void write_geant4_macro_cylindrical(const Eigen::Vector3d& center,
                                     double radius, double half_length,
                                     double energy_keV,
                                     uint64_t num_events,
                                     const std::string& filename);

/// Write a GEANT4 GPS macro for a rectangular (box) volume source.
/// Uses GPS /gps/pos/type Volume with /gps/pos/shape Para.
void write_geant4_macro_rectangular(const Eigen::Vector3d& center,
                                     const Eigen::Vector3d& half_dims,
                                     double energy_keV,
                                     uint64_t num_events,
                                     const std::string& filename);

/// Write a GEANT4 GPS macro for a spherical volume source. Uses GPS
/// /gps/pos/type Volume with /gps/pos/shape Sphere, generating in the outer
/// ball and confining to SrcMaterialPV — for a hollow shell (rmin > 0) the
/// void points are rejected by the confinement, so `radius` is the OUTER radius.
void write_geant4_macro_spherical(const Eigen::Vector3d& center,
                                   double radius,
                                   double energy_keV,
                                   uint64_t num_events,
                                   const std::string& filename);

} // namespace ceelo
