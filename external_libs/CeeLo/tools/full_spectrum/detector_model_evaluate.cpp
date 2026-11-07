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

#include "full_spectrum/ArrivalResponse.h"
#include "full_spectrum/FullSpectrumModel.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ceelo::full_spectrum::SpectralComponent;

const char* component_name(SpectralComponent component)
{
    switch (component) {
    case SpectralComponent::DetectorZeroDeposit:
        return "zero_deposit";
    case SpectralComponent::DetectorFullEnergy:
        return "full_energy";
    case SpectralComponent::DetectorPairFullEnergy:
        return "pair_full_energy";
    case SpectralComponent::DetectorSingleEscape:
        return "single_escape";
    case SpectralComponent::DetectorDoubleEscape:
        return "double_escape";
    case SpectralComponent::DetectorAnnihilation511:
        return "annihilation_511";
    case SpectralComponent::DetectorKAlphaEscape:
        return "k_alpha_escape";
    case SpectralComponent::DetectorKBetaEscape:
        return "k_beta_escape";
    case SpectralComponent::DetectorKAlphaAbsorbed:
        return "k_alpha_absorbed";
    case SpectralComponent::DetectorKBetaAbsorbed:
        return "k_beta_absorbed";
    case SpectralComponent::DetectorPhotoelectricEscape:
        return "photoelectric_escape";
    case SpectralComponent::SingleCompton:
        return "single_compton";
    case SpectralComponent::MultipleCompton:
        return "multiple_compton";
    case SpectralComponent::DetectorPairNoCompton:
        return "pair_no_compton";
    case SpectralComponent::DetectorPairSingleCompton:
        return "pair_single_compton";
    case SpectralComponent::DetectorPairMultipleCompton:
        return "pair_multiple_compton";
    default:
        return "unsupported";
    }
}

struct Config {
    std::string payload;
    std::string output = "detector_model.csv";
    double energy_keV = 1332.492;
    double cos_theta = 1.0;
    // An arrival class (iso, front, back, side, mu-band:LO,HI) instead of
    // cos_theta: the payload averaged over the band, as the detector
    // campaign's --incidence samples it.
    std::string incidence;
    std::size_t nodes_per_hemisphere = 16;
    int bins = 256;
};

Config parse(int argc, char** argv)
{
    Config config;
    auto value = [&](int& index) {
        if (++index >= argc)
            throw std::invalid_argument("missing option value");
        return std::string(argv[index]);
    };
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--payload")
            config.payload = value(index);
        else if (option == "--output")
            config.output = value(index);
        else if (option == "--energy-kev")
            config.energy_keV = std::stod(value(index));
        else if (option == "--cos-theta")
            config.cos_theta = std::stod(value(index));
        else if (option == "--incidence")
            config.incidence = value(index);
        else if (option == "--nodes-per-hemisphere")
            config.nodes_per_hemisphere = std::stoul(value(index));
        else if (option == "--bins")
            config.bins = std::stoi(value(index));
        else
            throw std::invalid_argument("unknown option: " + option);
    }
    if (config.payload.empty() || !(config.energy_keV > 0.0)
        || config.cos_theta < -1.0 || config.cos_theta > 1.0
        || config.bins < 1 || config.nodes_per_hemisphere < 1) {
        throw std::invalid_argument(
            "invalid detector-model evaluation configuration");
    }
    return config;
}

std::string read_text(const std::string& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open " + path);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Config config = parse(argc, argv);
        const auto model =
            ceelo::full_spectrum::DetectorPulseResponseModel::deserialize(
                read_text(config.payload));
        std::vector<double> edges(
            static_cast<std::size_t>(config.bins + 1));
        for (int index = 0; index <= config.bins; ++index) {
            edges[static_cast<std::size_t>(index)] =
                config.energy_keV * index / config.bins;
        }
        const bool arrival = !config.incidence.empty();
        ceelo::full_spectrum::ArrivalBand band;
        double mean_area_cm2 = 0.0;
        ceelo::full_spectrum::SpectralMeasure<double> measure;
        if (arrival) {
            band = ceelo::full_spectrum::arrival_band(config.incidence);
            const auto& anchor =
                model.payload().detected_rate_surface.anchor;
            mean_area_cm2 = ceelo::full_spectrum::arrival_quadrature(
                band, anchor.envelope_radius_cm, anchor.envelope_length_cm,
                config.nodes_per_hemisphere).mean_projected_area_cm2;
            measure = ceelo::full_spectrum::evaluate_arrival(
                model, config.energy_keV, band, edges,
                config.nodes_per_hemisphere);
        } else {
            measure = model.evaluate(
                config.energy_keV, config.cos_theta, edges);
        }

        // The campaign's detector column: the payload's own label (the
        // hand-set pilots carry none).
        const std::string detector = model.payload().detector_name.empty()
            ? "pilot" : model.payload().detector_name;
        std::ofstream output(config.output);
        if (!output)
            throw std::runtime_error("cannot open " + config.output);
        // The campaign's layout, arrival columns included (the envelope's
        // exact mean projected area, sigma 0).
        output << std::setprecision(17)
            << "cell_kind,energy_keV,Z,tau,areal_density_g_cm2,"
               "cos_theta,detector,component,bin_lo_keV,bin_hi_keV,"
               "sum_w,sum_w2,events,seed,analytic_direct"
            << (arrival ? ",incidence,mu_lo,mu_hi,mean_projected_area_cm2,"
                          "mean_projected_area_sigma_cm2"
                        : "")
            << '\n';
        std::ostringstream suffix_text;
        suffix_text << std::setprecision(17);
        if (arrival) {
            suffix_text << ',' << config.incidence << ',' << band.mu_lo << ','
                << band.mu_hi << ',' << mean_area_cm2 << ",0";
        }
        const std::string suffix = suffix_text.str();
        auto write_row = [&](SpectralComponent component, double lo,
                             double hi, double value) {
            output << (arrival ? "detector_arrival," : "detector,")
                   << config.energy_keV << ",,,,";
            if (!arrival)
                output << config.cos_theta;
            output << ',' << detector << ',' << component_name(component)
                   << ',' << lo << ',' << hi << ',' << value << ",0,1,0,"
                   << suffix << '\n';
        };
        for (const auto& atom : measure.atoms)
            write_row(atom.component, atom.energy_keV, atom.energy_keV,
                      atom.yield);
        for (const auto& bin : measure.continuous)
            write_row(bin.component, bin.lower_keV, bin.upper_keV,
                      bin.integral);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_detector_model_evaluate: "
                  << error.what() << '\n';
        return 2;
    }
}
