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

#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : gps_(new G4GeneralParticleSource())
{}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete gps_;
    delete gun_;
}

void PrimaryGeneratorAction::EnableConeBias(
    const G4ThreeVector& source_pos,
    double energy_MeV,
    double det_radius_mm,
    double det_z_min_mm,
    double det_z_max_mm)
{
    cone_bias_enabled_ = true;
    source_pos_ = source_pos;
    energy_MeV_ = energy_MeV;

    // Cone axis: from source toward detector center
    double det_center_z = 0.5 * (det_z_min_mm + det_z_max_mm);
    cone_axis_ = G4ThreeVector(0.0, 0.0, det_center_z) - source_pos;
    cone_axis_ = cone_axis_.unit();

    // Cone half-angle: must subtend the full detector bounding cylinder.
    // Use the apparent radius from the source position.
    double dx = source_pos.x();
    double dy = source_pos.y();
    double radial_offset = std::sqrt(dx * dx + dy * dy);
    double max_apparent_radius = det_radius_mm + radial_offset;

    // Distance along cone axis from source to nearest face
    double dz_front = std::abs(source_pos.z() - det_z_min_mm);
    double dz_back  = std::abs(source_pos.z() - det_z_max_mm);
    double dz = std::min(dz_front, dz_back);

    // Half-angle with 5% margin
    double half_angle = std::atan2(max_apparent_radius, dz) * 1.05;
    if (half_angle > M_PI) half_angle = M_PI;

    cos_theta_max_ = std::cos(half_angle);
    cone_weight_ = (1.0 - cos_theta_max_) / 2.0;

    // Set up particle gun
    gun_ = new G4ParticleGun(1);
    gun_->SetParticleDefinition(G4Gamma::Gamma());
    gun_->SetParticleEnergy(energy_MeV * MeV);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    if (!cone_bias_enabled_) {
        // Default GPS mode
        event_weight_ = 1.0;
        gps_->GeneratePrimaryVertex(event);
        return;
    }

    // Cone-bias mode: sample direction uniformly within cone
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double cos_theta = 1.0 - uniform(rng_) * (1.0 - cos_theta_max_);
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = 2.0 * M_PI * uniform(rng_);

    // Direction in cone-axis frame (cone_axis_ = local z)
    G4ThreeVector local_dir(sin_theta * std::cos(phi),
                            sin_theta * std::sin(phi),
                            cos_theta);

    // Rotate from local frame (z = cone_axis_) to world frame
    // Build orthonormal basis: cone_axis_ = w, u, v
    G4ThreeVector w = cone_axis_;
    G4ThreeVector u, v;
    if (std::abs(w.x()) < 0.9) {
        u = G4ThreeVector(1, 0, 0).cross(w).unit();
    } else {
        u = G4ThreeVector(0, 1, 0).cross(w).unit();
    }
    v = w.cross(u);

    G4ThreeVector world_dir = local_dir.x() * u
                            + local_dir.y() * v
                            + local_dir.z() * w;

    gun_->SetParticlePosition(source_pos_);
    gun_->SetParticleMomentumDirection(world_dir);
    gun_->GeneratePrimaryVertex(event);

    event_weight_ = cone_weight_;
}
