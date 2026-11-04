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

// Build Ba-133, Cs-137, Co-60 and Eu-152 spectra from SandiaDecay emissions
// and the differentiable FullSpecApprox source and detector kernels.

#include "full_spectrum/FullSpectrumModel.h"
#include "cascade/BetaSpectrum.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = ceelo::full_spectrum;

namespace {

std::vector<fs::EmissionLine> lines_for(
    const SandiaDecay::SandiaDecayDataBase& database,
    const std::string& label,
    double age_seconds)
{
    const SandiaDecay::Nuclide* nuclide = database.nuclide(label);
    if (!nuclide)
        throw std::runtime_error("Unknown SandiaDecay nuclide: " + label);
    SandiaDecay::NuclideMixture mixture;
    if (age_seconds < 0.0)
        mixture.addNuclideInPromptEquilibrium(nuclide, 1.0);
    else
        mixture.addAgedNuclideByActivity(nuclide, 1.0, age_seconds);
    const auto gammas = mixture.gammas(
        0.0, SandiaDecay::NuclideMixture::OrderByEnergy,
        /*includeAnnihilation=*/false);
    std::vector<fs::EmissionLine> lines;
    for (const auto& gamma : gammas) {
        if (gamma.energy >= 10.0 && gamma.energy <= 12000.0
            && gamma.numPerSecond > 1e-8) {
            // The mixture was normalized to one Bq of parent activity, so
            // photons/second is photons/parent-decay for this fixed-age list.
            lines.push_back({gamma.energy, gamma.numPerSecond});
        }
    }
    return lines;
}

std::string read_text(const std::string& path)
{
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("cannot open payload: " + path);
    std::ostringstream text;
    text << stream.rdbuf();
    return text.str();
}

std::vector<fs::ChargedParticleLine> charged_lines_for(
    const SandiaDecay::SandiaDecayDataBase& database,
    const std::string& label,
    int quadrature_bins)
{
    const auto emissions =
        ceelo::cascade_adapter::build_radioactive_emissions(database, label);
    std::vector<fs::ChargedParticleLine> lines;
    for (const auto& branch : emissions.beta_branches) {
        if (!(branch.endpoint_keV > 0.0)
            || !(branch.yield_per_decay > 0.0)) {
            continue;
        }
        const ceelo::BetaSpectrumSampler spectrum(branch);
        const double width = branch.endpoint_keV / quadrature_bins;
        std::vector<double> weights(
            static_cast<std::size_t>(quadrature_bins));
        double normalization = 0.0;
        for (int index = 0; index < quadrature_bins; ++index) {
            const double energy = (index + 0.5) * width;
            const double weight =
                spectrum.normalized_density_per_keV(energy) * width;
            weights[static_cast<std::size_t>(index)] = weight;
            normalization += weight;
        }
        if (!(normalization > 0.0))
            continue;
        for (int index = 0; index < quadrature_bins; ++index) {
            lines.push_back({
                (index + 0.5) * width,
                branch.yield_per_decay
                    * weights[static_cast<std::size_t>(index)]
                    / normalization,
                branch.is_positron});
        }
    }
    for (const auto& line : emissions.conversion_electrons) {
        if (line.energy_keV > 0.0 && line.yield_per_decay > 0.0) {
            lines.push_back({
                line.energy_keV, line.yield_per_decay, false});
        }
    }
    return lines;
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::string output = "full_spectrum_nuclides.csv";
        double age_seconds = -1.0; // prompt equilibrium
        double areal_density = 0.0;
        double effective_Z = 26.0;
        double decays = 1.0e6;
        int beta_quadrature_bins = 64;
        int source_quadrature_bins = 8;
        std::string detector_payload_path;
        std::string charged_payload_path;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto value = [&]() -> const char* {
                if (++i >= argc)
                    throw std::invalid_argument("missing option value");
                return argv[i];
            };
            if (arg == "--output")
                output = value();
            else if (arg == "--age-seconds")
                age_seconds = std::stod(value());
            else if (arg == "--areal-density")
                areal_density = std::stod(value());
            else if (arg == "--effective-Z")
                effective_Z = std::stod(value());
            else if (arg == "--decays")
                decays = std::stod(value());
            else if (arg == "--beta-quadrature-bins")
                beta_quadrature_bins = std::stoi(value());
            else if (arg == "--source-quadrature-bins")
                source_quadrature_bins = std::stoi(value());
            else if (arg == "--detector-payload")
                detector_payload_path = value();
            else if (arg == "--charged-particle-payload")
                charged_payload_path = value();
            else
                throw std::invalid_argument("unknown option: " + arg);
        }

        if (beta_quadrature_bins < 8)
            throw std::invalid_argument(
                "beta quadrature needs at least eight bins");
        const SandiaDecay::SandiaDecayDataBase database(
            SANDIA_DECAY_XML_PATH);
        fs::FullSpectrumPayload payload;
        // This research setting is the measured quality/runtime compromise
        // used by the current ~seconds-scale 100-line benchmark.
        payload.source_quadrature_bins = source_quadrature_bins;
        payload.distance_cm = 100.0;
        payload.incidence_cos_theta = 1.0;
        payload.detector = detector_payload_path.empty()
            ? fs::coaxial_hpge_payload()
            : fs::DetectorPulseResponseModel::deserialize(
                  read_text(detector_payload_path)).payload();
        const std::optional<fs::ChargedParticleSourceFluxModel>
            charged_source = charged_payload_path.empty()
                ? std::nullopt
                : std::optional<fs::ChargedParticleSourceFluxModel>(
                    fs::ChargedParticleSourceFluxModel::deserialize(
                        read_text(charged_payload_path)));
        const fs::FastSpectrumModel model = charged_source
            ? fs::FastSpectrumModel(payload, *charged_source)
            : fs::FastSpectrumModel(payload);
        std::vector<double> channels(4097);
        for (std::size_t i = 0; i < channels.size(); ++i)
            channels[i] = 12000.0 * static_cast<double>(i) / 4096.0;
        const fs::ShieldParameters<double> shield{
            areal_density, effective_Z};
        // Representative coaxial-HPGe resolution. Calibration/tails remain
        // outside v1 as stated in the campaign plan.
        const fs::ResolutionParameters<double> resolution{
            0.21, 0.0023, 0.0};

        std::ofstream stream(output);
        if (!stream)
            throw std::runtime_error("cannot open output: " + output);
        stream << "# FullSpecApprox SandiaDecay spectra\n"
               << "# geometry: 100 cm, front-on, vacuum\n"
               << "# age_seconds: " << age_seconds << "\n"
               << "# detector_payload: "
               << (detector_payload_path.empty()
                       ? "mechanistic-pilot" : detector_payload_path)
               << "\n"
               << "# charged_particle_payload: "
               << (charged_payload_path.empty()
                       ? "disabled" : charged_payload_path)
               << "\n"
               << "nuclide,channel,lower_keV,upper_keV,expected_counts\n";
        for (const std::string label :
             {"Ba133", "Cs137", "Co60", "Eu152"}) {
            const auto lines = lines_for(database, label, age_seconds);
            const auto charged_lines = charged_source
                ? charged_lines_for(
                    database, label, beta_quadrature_bins)
                : std::vector<fs::ChargedParticleLine>{};
            const auto spectrum = model.predict(
                lines, charged_lines, decays, shield, resolution, channels);
            for (std::size_t channel = 0;
                 channel < spectrum.size(); ++channel) {
                stream << label << ',' << channel << ','
                       << channels[channel] << ','
                       << channels[channel + 1] << ','
                       << spectrum[channel] << '\n';
            }
            std::cout << label << ": " << lines.size()
                      << " gamma lines, " << charged_lines.size()
                      << " charged quadrature lines\n";
        }
        std::cout << "wrote " << output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "full_spectrum_nuclides: " << error.what() << '\n';
        return 2;
    }
}
