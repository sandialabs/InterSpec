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

/// @file Material.h
/// @brief Material definitions for photon transport.
///
/// A Material represents a substance with known elemental composition and density.
/// It computes macroscopic cross-sections (Σ = ρ * N_A * Σ_i(w_i/A_i * σ_i))
/// for use in photon transport.
///
/// Material objects are intended to be constructed once and then used as const
/// throughout the simulation (thread-safe by design).

#include <string>
#include <vector>
#include <cstdint>

namespace ceelo {

/// A single element in a material composition.
struct MaterialComponent {
    uint8_t Z;              ///< Atomic number
    double mass_fraction;   ///< Mass fraction (should sum to 1 across all components)
};

/// Macroscopic cross-sections for a material at a given energy.
/// Units: 1/cm (inverse centimeters).
struct MacroscopicXS {
    double mu_pe;       ///< Photoelectric
    double mu_cs;       ///< Compton scattering
    double mu_rs;       ///< Rayleigh scattering
    double mu_pp;       ///< Pair production

    double mu_total() const { return mu_pe + mu_cs + mu_rs + mu_pp; }
};

/// Represents a material with known composition and density.
class Material {
public:
    /// Construct a material from its name, density, and elemental composition.
    /// @param name  Human-readable name (e.g., "NaI(Tl)")
    /// @param density_g_per_cm3  Bulk density in g/cm³
    /// @param composition  Elements and their mass fractions (must sum to ~1.0)
    Material(std::string name, double density_g_per_cm3,
             std::vector<MaterialComponent> composition);

    /// Compute macroscopic cross-sections at a given energy.
    /// The formula is: Σ_type = ρ * N_A * Σ_i (w_i / A_i) * σ_type_i(E)
    /// where w_i is mass fraction, A_i is atomic weight, σ_type_i is microscopic XS.
    /// @param energy_MeV  Photon energy in MeV
    /// @return Macroscopic cross-sections in 1/cm
    MacroscopicXS macroscopic_xs(double energy_MeV) const;

    /// Total linear attenuation coefficient mu_total in 1/cm.
    double mu_total(double energy_MeV) const;

    /// Convenience: mass attenuation coefficient mu/rho in cm²/g.
    double mass_attenuation(double energy_MeV) const;

    /// Select which element in the material was hit, given that an interaction occurred.
    /// Returns the Z of the interacting element, sampled proportionally to each
    /// element's contribution to the specified cross-section type.
    /// @param energy_MeV  Photon energy
    /// @param interaction_type  0=PE, 1=CS, 2=RS, 3=PP
    /// @param xi  Random number in [0, 1)
    /// @return Atomic number Z of the selected element
    int select_element(double energy_MeV, int interaction_type, double xi) const;

    // Accessors
    const std::string& name() const { return name_; }
    double density() const { return density_g_per_cm3_; }
    const std::vector<MaterialComponent>& composition() const { return composition_; }

    /// Number density for element with index i (atoms/cm³).
    /// n_i = ρ * N_A * w_i / A_i
    double number_density(size_t component_index) const;

private:
    std::string name_;
    double density_g_per_cm3_;
    std::vector<MaterialComponent> composition_;

    // Pre-computed: ρ * N_A * w_i / A_i for each component (atoms/cm³ / barn-to-cm²)
    // Actually: ρ * N_A * w_i / A_i  in atoms/(cm³), then σ (barn) * 1e-24 → Σ (1/cm)
    std::vector<double> number_densities_;  // atoms/cm³ for each component
};

// Avogadro's number
constexpr double kAvogadro = 6.02214076e23;  // atoms/mol

// Barn to cm² conversion
constexpr double kBarnToCm2 = 1.0e-24;

// --- Built-in material factory functions ---
// Defined in BuiltinMaterials.cpp

/// NaI (Tl-doped sodium iodide) — most common scintillator
Material make_NaI();

/// LaBr3 (Ce-doped lanthanum bromide) — high-resolution scintillator
Material make_LaBr3();

/// CeBr3 (cerium bromide) — similar to LaBr3
Material make_CeBr3();

/// HPGe (high-purity germanium)
Material make_HPGe();

/// CZT (cadmium zinc telluride, Cd0.9Zn0.1Te)
Material make_CZT();

/// CLYC (Cs2LiYCl6, Ce-doped)
Material make_CLYC();

// --- Common shielding materials ---

Material make_Lead();
Material make_Copper();
Material make_Iron();
Material make_Aluminum();
Material make_Tin();
Material make_Tungsten();
Material make_StainlessSteel304();

// --- Environmental / source materials ---

/// Dry soil (average continental crust composition), density 1.5 g/cm³.
/// Elemental composition from IAEA-TECDOC-1011 (1998) / CRC Handbook:
///   O(46%), Si(28%), Al(8%), Fe(5%), Ca(4%), Mg(2%), Na(2%), K(2%), Ti(0.5%), H(0.5%).
/// Use for extended environmental source geometries (cylinder, box of soil).
Material make_Soil();

/// Water (H₂O), density 1.0 g/cm³.
Material make_Water();

/// Polyethylene (C₂H₄)ₙ, density 0.94 g/cm³ (HDPE).
/// Common Marinelli beaker container material.
Material make_Polyethylene();

/// Pyrex (Corning 7740 borosilicate glass), density 2.23 g/cm³.
/// Composition: SiO₂(80.6%), B₂O₃(13%), Na₂O(4%), Al₂O₃(2.3%), K₂O(0.1%).
Material make_Pyrex();

/// Cellulose (C₆H₁₀O₅)ₙ, density 0.5 g/cm³ (loose powder/insulation).
Material make_Cellulose();

/// Dry air at sea level, density 0.001205 g/cm³ (20°C, 1 atm).
/// Composition: N(75.53%), O(23.18%), Ar(1.28%), C(0.01%) by mass.
Material make_Air();

} // namespace ceelo
