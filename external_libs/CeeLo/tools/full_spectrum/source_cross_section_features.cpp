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

#include "full_spectrum/FullSpectrumModel.h"

#include "cross_sections/CrossSectionData.h"
#include "physics/ElectronCsda.h"

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

int main(int argc, char** argv) {
    // --with-hydrogen-column: each input row is "E Z T f_H" and every output
    // quantity is the two-component mixture (f_H hydrogen + remainder Z):
    // photon processes mix by mass fraction, stopping power adds, ranges and
    // radiation length combine reciprocally (Bragg additivity). Without the
    // flag each input row is "E Z T" (a single element).
    // --emit-process-features: instead of the raw cross-section columns, emit
    // the FINISHED 7-D "process" feature vector for each secondary component,
    // straight out of detail::source_process_features -- the same function the
    // C++ evaluator calls. The Python fitter reads these rather than
    // recomputing the arithmetic, which makes fitter/evaluator lockstep
    // structural instead of a convention that can silently drift (the basis
    // normalization drifting that way once cost an 83% shape-TV error).
    // Implies --with-hydrogen-column input (E Z T f_H).
    // --smooth-z-blend: report mu_pe/mu_cs/mu_rs/mu_pp from the same Gaussian
    // effective-Z cardinal blend the evaluator uses (SmoothElementCrossSection,
    // width 0.42) rather than the exact element. At an integer Z the two differ
    // by the ~0.7% neighbour contamination the blend carries, which is enough
    // to matter near an absorption edge. The Python fitter's mechanistic bases
    // must be built from the SAME numbers the evaluator will divide back out,
    // or the fitted conditional yield absorbs the difference at fit time and
    // re-introduces it at evaluation time.
    // --emit-component-bases: emit the FINISHED mechanistic baseline of each
    // secondary component, straight out of detail::source_component_bases --
    // the same function the evaluator multiplies back in. The Python fitter
    // divides these OUT at fit time, so a second copy of the arithmetic there
    // would be a purely multiplicative bias that no fit metric can see.
    // Implies --with-hydrogen-column input (E Z T f_H).
    // --fluorescence-base-terms N: which vacancy sources the fluorescence
    // bases include, a bit set (1 direct photoabsorption, 2 once-Compton then
    // absorbed, 4 absorbed annihilation photon; default 7, what the evaluator
    // uses) -- for kernel-span studies only.
    //
    // The emitted bases assume the campaign's physics: Rayleigh ON (use
    // --no-rayleigh to gate it, matching SourceModelPayload::enable_rayleigh)
    // and photonuclear OFF, which is the payload default. A payload with
    // enable_photonuclear = true adds mu_pn to mu0_total and would need that
    // term here too.
    bool with_hydrogen = false;
    bool emit_process_features = false;
    bool emit_component_bases = false;
    bool smooth_z_blend = false;
    bool enable_rayleigh = true;
    int fluorescence_terms = ceelo::full_spectrum::detail::kFluorescenceAllTerms;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--fluorescence-base-terms" && i + 1 < argc)
            fluorescence_terms = std::stoi(argv[++i]);
        else if (option == "--with-hydrogen-column")
            with_hydrogen = true;
        else if (option == "--smooth-z-blend")
            smooth_z_blend = true;
        else if (option == "--no-rayleigh")
            enable_rayleigh = false;
        else if (option == "--emit-process-features") {
            emit_process_features = true;
            with_hydrogen = true;
        } else if (option == "--emit-component-bases") {
            emit_component_bases = true;
            with_hydrogen = true;
        } else {
            std::cerr << "unknown option: " << argv[i] << '\n';
            return 2;
        }
    }
    if (emit_process_features && emit_component_bases) {
        std::cerr << "--emit-process-features and --emit-component-bases are "
                     "separate modes\n";
        return 2;
    }
    constexpr double electron_mass_keV = 510.99895;

    if (emit_component_bases) {
        using ceelo::full_spectrum::detail::kSourceComponentCount;
        // Must match SourceModelPayload::z_blend_width; the evaluator builds
        // its SmoothElementCrossSection from that same value.
        constexpr double kZBlendWidth = 0.42;
        const ceelo::full_spectrum::SmoothElementCrossSection
            cross_sections(kZBlendWidth);
        ceelo::full_spectrum::SourceSecondarySurrogatePayload options;
        static const std::array<const char*, kSourceComponentCount> names{{
            "pair_511_atom",
            "pair_annihilation_continuum",
            "pair_annihilation_in_flight",
            "photoelectron_brems",
            "compton_electron_brems",
            "pair_electron_brems",
            "fluorescence_k_alpha_atom",
            "fluorescence_k_beta_atom",
            "fluorescence_l_atom",
        }};
        std::cout << "energy_keV,Z,areal_density_g_cm2,"
                     "hydrogen_mass_fraction,component,base\n"
                  << std::setprecision(17);
        double energy = 0.0, thickness = 0.0, fraction = 0.0;
        int z = 0;
        while (std::cin >> energy >> z >> thickness >> fraction) {
            if (!(energy > 0.0) || z < 1 || z > 92 || thickness < 0.0
                || fraction < 0.0 || fraction >= 1.0) {
                std::cerr << "invalid E, Z, areal-density, or f_H input\n";
                return 2;
            }
            const double effective_Z = static_cast<double>(z);
            const auto mu0 =
                cross_sections.evaluate(energy, effective_Z, fraction);
            const double mu0_total =
                mu0[0] + mu0[1] + mu0[3] + (enable_rayleigh ? mu0[2] : 0.0);
            const auto bases =
                ceelo::full_spectrum::detail::source_component_bases(
                    energy, thickness, effective_Z, fraction, mu0, mu0_total,
                    cross_sections, options, fluorescence_terms);
            for (std::size_t index = 0; index < kSourceComponentCount;
                 ++index) {
                std::cout << energy << ',' << z << ',' << thickness << ','
                          << fraction << ',' << names[index] << ','
                          << bases.base[index] << '\n';
            }
        }
        return 0;
    }

    if (emit_process_features) {
        using ceelo::full_spectrum::kSourceRbfMaxFeatures;
        using Process = ceelo::full_spectrum::detail::SourceSecondaryProcess;
        // Must match SourceModelPayload::z_blend_width; the evaluator builds
        // its SmoothElementCrossSection from that same value.
        constexpr double kZBlendWidth = 0.42;
        const ceelo::full_spectrum::SmoothElementCrossSection
            cross_sections(kZBlendWidth);
        const std::array<std::pair<const char*, Process>, 9> components{{
            {"pair_511_atom", Process::Pair511Atom},
            {"pair_annihilation_continuum",
             Process::PairAnnihilationContinuum},
            {"pair_annihilation_in_flight",
             Process::PairAnnihilationInFlight},
            {"photoelectron_brems", Process::PhotoelectronBrems},
            {"compton_electron_brems", Process::ComptonElectronBrems},
            {"pair_electron_brems", Process::PairElectronBrems},
            {"fluorescence_k_alpha_atom", Process::FluorescenceKAlpha},
            {"fluorescence_k_beta_atom", Process::FluorescenceKBeta},
            {"fluorescence_l_atom", Process::FluorescenceL},
        }};
        std::cout << "energy_keV,Z,areal_density_g_cm2,"
                     "hydrogen_mass_fraction,component";
        for (std::size_t index = 0; index < kSourceRbfMaxFeatures; ++index)
            std::cout << ",f" << index;
        std::cout << '\n' << std::setprecision(17);
        double energy = 0.0, thickness = 0.0, fraction = 0.0;
        int z = 0;
        while (std::cin >> energy >> z >> thickness >> fraction) {
            if (!(energy > 0.0) || z < 1 || z > 92 || thickness < 0.0
                || fraction < 0.0 || fraction >= 1.0) {
                std::cerr << "invalid E, Z, areal-density, or f_H input\n";
                return 2;
            }
            const double effective_Z = static_cast<double>(z);
            const auto mu0 =
                cross_sections.evaluate(energy, effective_Z, fraction);
            const double tau_total =
                thickness * (mu0[0] + mu0[1] + mu0[2] + mu0[3]);
            for (const auto& entry : components) {
                std::array<double, kSourceRbfMaxFeatures> feature{};
                ceelo::full_spectrum::detail::source_process_features(
                    energy, thickness, effective_Z, fraction, tau_total,
                    mu0, cross_sections, kZBlendWidth, entry.second,
                    feature);
                std::cout << energy << ',' << z << ',' << thickness << ','
                          << fraction << ',' << entry.first;
                for (double value : feature)
                    std::cout << ',' << value;
                std::cout << '\n';
            }
        }
        return 0;
    }
    std::cout
        << "energy_keV,Z,areal_density_g_cm2,"
           "mu_pe,mu_cs,mu_rs,mu_pp,"
           "electron_collision_stopping_MeV_cm2_g,"
           "electron_brems_rate_above_10keV_per_gcm2,"
           "moller_delta_rate_above_50keV_per_gcm2,"
           "radiation_length_g_cm2,electron_range_g_cm2,"
           "positron_range_g_cm2,"
           "reference3_keV,reference8_keV,"
           "tau_abs3,tau_scatter3,tau_abs8,tau_scatter8\n";
    std::cout << std::setprecision(17);

    double energy_keV = 0.0;
    int z_value = 0;
    double areal_density_g_cm2 = 0.0;
    double f_H = 0.0;
    const auto mix_xs = [](
        const ceelo::full_spectrum::PartialMassAttenuation& blend,
        const ceelo::full_spectrum::PartialMassAttenuation& hydrogen,
        double fraction) {
        ceelo::full_spectrum::PartialMassAttenuation out = blend;
        out.pe += fraction * (hydrogen.pe - blend.pe);
        out.cs += fraction * (hydrogen.cs - blend.cs);
        out.rs += fraction * (hydrogen.rs - blend.rs);
        out.pp += fraction * (hydrogen.pp - blend.pp);
        return out;
    };
    const auto reciprocal_mix = [](
        double blend, double hydrogen, double fraction) {
        return blend / (1.0 + fraction * (blend / hydrogen - 1.0));
    };
    while (true) {
        if (!(std::cin >> energy_keV >> z_value >> areal_density_g_cm2))
            break;
        if (with_hydrogen && !(std::cin >> f_H))
            break;
        if (!(energy_keV > 0.0) || z_value < 1 || z_value > 92
            || areal_density_g_cm2 < 0.0 || f_H < 0.0 || f_H >= 1.0) {
            std::cerr << "invalid E, Z, areal-density, or f_H input\n";
            return 2;
        }
        const ceelo::full_spectrum::SmoothElementCrossSection blended(0.42);
        const auto smooth_at = [&](double energy) {
            const auto mixed = blended.evaluate(
                energy, static_cast<double>(z_value), f_H);
            ceelo::full_spectrum::PartialMassAttenuation out;
            out.pe = mixed[0];
            out.cs = mixed[1];
            out.rs = mixed[2];
            out.pp = mixed[3];
            return out;
        };
        const auto element_xs =
            ceelo::full_spectrum::element_mass_cross_sections(
                z_value, energy_keV);
        const auto incident =
            smooth_z_blend
                ? smooth_at(energy_keV)
                : mix_xs(
                    element_xs,
                    ceelo::full_spectrum::element_mass_cross_sections(
                        1, energy_keV),
                    f_H);
        const double element_radiation_length =
            ceelo::ElectronCsda::radiation_length_gcm2_element(
                z_value,
                ceelo::CrossSectionData::instance().atomic_weight(z_value));
        const double radiation_length = reciprocal_mix(
            element_radiation_length,
            ceelo::ElectronCsda::radiation_length_gcm2_element(
                1, ceelo::CrossSectionData::instance().atomic_weight(1)),
            f_H);
        const double element_stopping =
            ceelo::ElectronCsda::stopping_power_MeV_cm2_g(
                z_value,
                ceelo::CrossSectionData::instance().atomic_weight(z_value),
                energy_keV, false);
        const double electron_stopping = element_stopping
            + f_H * (ceelo::ElectronCsda::stopping_power_MeV_cm2_g(
                         1,
                         ceelo::CrossSectionData::instance().atomic_weight(1),
                         energy_keV, false)
                     - element_stopping);
        std::vector<ceelo::MaterialComponent> composition;
        if (f_H > 0.0)
            composition.push_back({std::uint8_t{1}, f_H});
        composition.push_back(
            {static_cast<std::uint8_t>(z_value), 1.0 - f_H});
        ceelo::Material element("equivalent-element", 1.0, composition);
        const double brems_rate =
            ceelo::ElectronCsda::instance()
                .brems_emission_rate_per_gcm2(
                    element, energy_keV, false);
        const double delta_rate =
            ceelo::ElectronCsda::instance()
                .delta_electron_rate_per_gcm2(
                    element, energy_keV, 50.0, false);
        const double electron_range = reciprocal_mix(
            ceelo::ElectronCsda::instance().range_gcm2(
                z_value, energy_keV, false),
            ceelo::ElectronCsda::instance().range_gcm2(
                1, energy_keV, false),
            f_H);
        const double positron_range = reciprocal_mix(
            ceelo::ElectronCsda::instance().range_gcm2(
                z_value, energy_keV, true),
            ceelo::ElectronCsda::instance().range_gcm2(
                1, energy_keV, true),
            f_H);
        const double reference3 = energy_keV
            / (1.0 + 3.0 * energy_keV / electron_mass_keV);
        const double reference8 = energy_keV
            / (1.0 + 8.0 * energy_keV / electron_mass_keV);
        const auto xs3 =
            smooth_z_blend
                ? smooth_at(reference3)
                : mix_xs(
                    ceelo::full_spectrum::element_mass_cross_sections(
                        z_value, reference3),
                    ceelo::full_spectrum::element_mass_cross_sections(
                        1, reference3),
                    f_H);
        const auto xs8 =
            smooth_z_blend
                ? smooth_at(reference8)
                : mix_xs(
                    ceelo::full_spectrum::element_mass_cross_sections(
                        z_value, reference8),
                    ceelo::full_spectrum::element_mass_cross_sections(
                        1, reference8),
                    f_H);
        std::cout
            << energy_keV << ',' << z_value << ','
            << areal_density_g_cm2 << ','
            << incident.pe << ',' << incident.cs << ','
            << incident.rs << ',' << incident.pp << ','
            << electron_stopping << ','
            << brems_rate << ',' << delta_rate << ','
            << radiation_length << ',' << electron_range << ','
            << positron_range << ','
            << reference3 << ',' << reference8 << ','
            << areal_density_g_cm2 * (xs3.pe + xs3.pp) << ','
            << areal_density_g_cm2 * (xs3.cs + xs3.rs) << ','
            << areal_density_g_cm2 * (xs8.pe + xs8.pp) << ','
            << areal_density_g_cm2 * (xs8.cs + xs8.rs) << '\n';
    }
    return 0;
}
