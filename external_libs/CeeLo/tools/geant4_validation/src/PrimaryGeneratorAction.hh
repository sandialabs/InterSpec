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

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"

#include <cstdint>
#include <random>

class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;

/// Primary generator with optional cone biasing.
///
/// In default (GPS) mode, uses GEANT4's General Particle Source driven by
/// macro commands.  In cone-bias mode, samples photon directions uniformly
/// within a bounding cone toward the detector, giving ~1/(cone solid angle
/// fraction) more detector hits per event.
///
/// The importance weight (= cone solid angle fraction = (1-cos_theta)/2)
/// is stored per event and must be applied to all tallies.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

    /// Enable cone-bias mode.  Call before Initialize().
    /// @param source_pos  Source position (mm, G4 internal units)
    /// @param energy_MeV  Photon energy in MeV
    /// @param det_radius_mm  Bounding radius of outermost geometry (mm)
    /// @param det_z_min_mm   Front face z in mm (typically 0 for our GDML)
    /// @param det_z_max_mm   Back face z in mm
    void EnableConeBias(const G4ThreeVector& source_pos,
                        double energy_MeV,
                        double det_radius_mm,
                        double det_z_min_mm,
                        double det_z_max_mm);

    /// Importance weight for the current event (1.0 in GPS mode).
    double GetEventWeight() const { return event_weight_; }

    /// Is cone bias enabled?
    bool IsConeBiasEnabled() const { return cone_bias_enabled_; }

    /// Cone solid angle as fraction of 4π (0 if not cone-biased).
    double GetConeSolidAngleFraction() const { return cone_weight_; }

    /// Multiply the auto-computed cone half-angle by this factor (default 1.0).
    /// Used by the air-scatter study to widen the cone so off-axis air that
    /// scatters photons into the detector is sampled. The importance weight
    /// (= cone solid-angle fraction) scales accordingly, so tallies stay
    /// unbiased. Call BEFORE EnableConeBias().
    void SetConeScale(double s) { cone_scale_ = s; }

    /// Seed the cone-direction generator (separate from G4Random).  main()
    /// passes the run's --seed, so seeded processes sample independent
    /// directions; a fixed seed made every parallel cone run repeat one
    /// direction sequence.
    void SetDirectionSeed(std::uint64_t seed) { rng_.seed(seed); }

private:
    // GPS mode (default)
    G4GeneralParticleSource* gps_ = nullptr;

    // Cone-bias mode
    bool cone_bias_enabled_ = false;
    G4ParticleGun* gun_ = nullptr;
    G4ThreeVector source_pos_;
    G4ThreeVector cone_axis_;  // Unit vector from source toward detector center
    double energy_MeV_ = 0.0;
    double cos_theta_max_ = 0.0;  // cos(cone half-angle)
    double cone_weight_ = 0.0;    // (1 - cos_theta_max) / 2
    double cone_scale_ = 1.0;     // multiplier on the auto half-angle
    double event_weight_ = 1.0;

    std::mt19937_64 rng_{42};
};
