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

#include "materials/Material.h"

namespace ceelo {

// Helper: compute mass fraction from stoichiometry
// Given pairs of (Z, count), compute mass fractions.
// E.g., NaI: {(11, 1), (53, 1)} → w_Na = 22.99/(22.99+126.90), w_I = 126.90/(22.99+126.90)
namespace {
struct Stoich { uint8_t Z; int count; double atomic_weight; };

std::vector<MaterialComponent> from_stoichiometry(std::initializer_list<Stoich> elements) {
    double total_mass = 0.0;
    for (const auto& e : elements) {
        total_mass += e.count * e.atomic_weight;
    }

    std::vector<MaterialComponent> result;
    result.reserve(elements.size());
    for (const auto& e : elements) {
        double w = (e.count * e.atomic_weight) / total_mass;
        result.push_back({e.Z, w});
    }
    return result;
}
} // anonymous namespace

// ---- Detector Materials ----

Material make_NaI() {
    // NaI (sodium iodide), density 3.67 g/cm³
    // Na: Z=11, A=22.990
    // I:  Z=53, A=126.904
    return Material("NaI", 3.67,
        from_stoichiometry({{11, 1, 22.990}, {53, 1, 126.904}}));
}

Material make_LaBr3() {
    // LaBr3 (lanthanum bromide), density 5.06 g/cm³
    // La: Z=57, A=138.905
    // Br: Z=35, A=79.904
    return Material("LaBr3", 5.06,
        from_stoichiometry({{57, 1, 138.905}, {35, 3, 79.904}}));
}

Material make_CeBr3() {
    // CeBr3 (cerium bromide), density 5.1 g/cm³
    // Ce: Z=58, A=140.116
    // Br: Z=35, A=79.904
    return Material("CeBr3", 5.1,
        from_stoichiometry({{58, 1, 140.116}, {35, 3, 79.904}}));
}

Material make_HPGe() {
    // High-purity germanium, density 5.323 g/cm³
    // Ge: Z=32, A=72.630
    return Material("HPGe", 5.323,
        {{32, 1.0}});
}

Material make_CZT() {
    // Cadmium zinc telluride: Cd₀.₉Zn₀.₁Te, density 5.78 g/cm³
    // Cd: Z=48, A=112.414
    // Zn: Z=30, A=65.380
    // Te: Z=52, A=127.600
    // Formula: Cd₀.₉Zn₀.₁Te → mass = 0.9*112.414 + 0.1*65.380 + 127.600
    double m_Cd = 0.9 * 112.414;
    double m_Zn = 0.1 * 65.380;
    double m_Te = 127.600;
    double total = m_Cd + m_Zn + m_Te;
    return Material("CZT", 5.78,
        {{48, m_Cd / total},
         {30, m_Zn / total},
         {52, m_Te / total}});
}

Material make_CLYC() {
    // Cs₂LiYCl₆ (CLYC), density 3.31 g/cm³
    // Cs: Z=55, A=132.905
    // Li: Z=3,  A=6.941
    // Y:  Z=39, A=88.906
    // Cl: Z=17, A=35.453
    return Material("CLYC", 3.31,
        from_stoichiometry({
            {55, 2, 132.905},  // Cs₂
            {3,  1, 6.941},    // Li
            {39, 1, 88.906},   // Y
            {17, 6, 35.453}    // Cl₆
        }));
}

// ---- Shielding Materials ----

Material make_Lead() {
    return Material("Pb", 11.35, {{82, 1.0}});
}

Material make_Copper() {
    return Material("Cu", 8.96, {{29, 1.0}});
}

Material make_Iron() {
    return Material("Fe", 7.874, {{26, 1.0}});
}

Material make_Aluminum() {
    return Material("Al", 2.70, {{13, 1.0}});
}

Material make_Tin() {
    return Material("Sn", 7.265, {{50, 1.0}});
}

Material make_Tungsten() {
    return Material("W", 19.3, {{74, 1.0}});
}

Material make_StainlessSteel304() {
    // SS304: ~70% Fe, ~19% Cr, ~10% Ni, ~1% Mn (approximate)
    // Fe: Z=26, Cr: Z=24, Ni: Z=28, Mn: Z=25
    return Material("SS304", 8.00,
        {{26, 0.70},    // Fe
         {24, 0.19},    // Cr
         {28, 0.10},    // Ni
         {25, 0.01}});  // Mn
}

// ---- Environmental / source materials ----

Material make_Soil() {
    // Dry soil — average continental crust composition.
    // Source: IAEA-TECDOC-1011 (1998) Table 1 / CRC Handbook (rounded).
    // Density: 1.5 g/cm³ (dry, loose — typical for garden/field soil).
    //
    // Elements and mass fractions:
    //   O  (Z= 8): 0.460   Si (Z=14): 0.276   Al (Z=13): 0.083
    //   Fe (Z=26): 0.050   Ca (Z=20): 0.037   Mg (Z=12): 0.021
    //   Na (Z=11): 0.021   K  (Z=19): 0.026   Ti (Z=22): 0.006
    //   H  (Z= 1): 0.020   (sum = 1.000)
    return Material("Soil", 1.5,
        {{ 8, 0.460},   // O
         {14, 0.276},   // Si
         {13, 0.083},   // Al
         {26, 0.050},   // Fe
         {20, 0.037},   // Ca
         {12, 0.021},   // Mg
         {11, 0.021},   // Na
         {19, 0.026},   // K
         {22, 0.006},   // Ti
         { 1, 0.020}}); // H
}

Material make_Water() {
    // Water (H₂O), density 1.0 g/cm³
    // H: Z=1, A=1.008  O: Z=8, A=15.999
    return Material("Water", 1.0,
        from_stoichiometry({{1, 2, 1.008}, {8, 1, 15.999}}));
}

Material make_Polyethylene() {
    // Polyethylene (C₂H₄)ₙ, density 0.94 g/cm³ (high-density HDPE)
    // C: Z=6, A=12.011  H: Z=1, A=1.008
    return Material("Polyethylene", 0.94,
        from_stoichiometry({{6, 2, 12.011}, {1, 4, 1.008}}));
}

Material make_Pyrex() {
    // Pyrex (Corning 7740 borosilicate glass), density 2.23 g/cm³.
    // Oxide composition: SiO₂ 80.6%, B₂O₃ 13.0%, Na₂O 4.0%, Al₂O₃ 2.3%, K₂O 0.1%.
    // Elemental mass fractions (computed from oxide stoichiometry):
    //   O  (Z= 8): 0.5396   Si (Z=14): 0.3768   B (Z= 5): 0.0404
    //   Na (Z=11): 0.0297   Al (Z=13): 0.0122   K (Z=19): 0.0013
    return Material("Pyrex", 2.23,
        {{ 8, 0.5396},   // O
         {14, 0.3768},   // Si
         { 5, 0.0404},   // B
         {11, 0.0297},   // Na
         {13, 0.0122},   // Al
         {19, 0.0013}}); // K
}

Material make_Cellulose() {
    // Cellulose (C₆H₁₀O₅)ₙ, density ~0.5 g/cm³ (loose cellulose powder/insulation).
    // C: Z=6, A=12.011  H: Z=1, A=1.008  O: Z=8, A=15.999
    return Material("Cellulose", 0.5,
        from_stoichiometry({{6, 6, 12.011}, {1, 10, 1.008}, {8, 5, 15.999}}));
}

Material make_Air() {
    // Dry air at sea level — ICRU Report 37 / NIST standard composition.
    // Density: 0.001205 g/cm³ (20°C, 1 atm).
    // Components by mass (matches G4_AIR):
    //   N  (Z= 7): 0.7553   O (Z= 8): 0.2318   Ar (Z=18): 0.0128   C (Z= 6): 0.0001
    return Material("Air", 0.001205,
        {{ 7, 0.7553},   // N
         { 8, 0.2318},   // O
         {18, 0.0128},   // Ar
         { 6, 0.0001}}); // C (trace)
}

} // namespace ceelo
