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

// Measure a detector payload's finite-distance depths (FiniteDistanceDepth)
// by Monte Carlo and write them into the payload. For each energy: a point
// source on the axis at a ladder of distances from the geometry's reference
// point (the crystal face or the endcap front, as for every user distance in
// CeeLo -- DetectorResponse::query_position -- and so for predict's), the
// full-energy and detected efficiencies of the payload's own geometry
// (EfficiencyCalculator), and for each the virtual disk -- depth delta and
// radius R -- that makes eps(d) s (s + x) constant across the ladder
// (x = d + delta, s = sqrt(x^2 + R^2): eps over the disk's solid angle, up to
// a constant), by grid search. A point whose efficiency has
// more than 0.3% MC noise is run again with up to 16x the events. A depth is
// kept only if every ladder point is detected with <= 2% noise and the fit's
// CV is <= 2%; otherwise (below ~30 keV on axis, where almost nothing gets
// through the front) the nearest kept energy's depth stands in, the higher
// one first.
//
//   ceelo_detector_finite_distance --payload IN --output OUT
//       [--energies E1,E2,...] [--distances D1,D2,...] [--events N]
//       [--threads T]

#include "efficiency/EfficiencyCalculator.h"
#include "full_spectrum/FullSpectrumModel.h"
#include "geometry/GeometryDescriptor.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = ceelo::full_spectrum;

namespace {

std::vector<double> number_list(const std::string& text)
{
    std::vector<double> out;
    std::stringstream list(text);
    std::string item;
    while (std::getline(list, item, ','))
        if (!item.empty())
            out.push_back(std::stod(item));
    return out;
}

struct DiskFit {
    double delta_cm = 0.0;
    double radius_cm = 0.0;
    double cv = 1e30;        ///< spread of eps / disk solid angle on the ladder
    double max_residual = 0.0;
};

DiskFit fit_disk(const std::vector<double>& distance, const std::vector<double>& eps,
                 double delta_max, double radius_max)
{
    DiskFit best;
    for (double radius = 0.0; radius <= radius_max + 1e-9; radius += 0.05) {
        for (double delta = 0.0; delta <= delta_max + 1e-9; delta += 0.01) {
            double s1 = 0.0, s2 = 0.0;
            for (std::size_t i = 0; i < distance.size(); ++i) {
                const double x = distance[i] + delta;
                const double s = std::sqrt(x * x + radius * radius);
                const double v = eps[i] * s * (s + x);
                s1 += v;
                s2 += v * v;
            }
            const double n = static_cast<double>(distance.size());
            const double mean = s1 / n;
            const double cv = mean > 0.0
                ? std::sqrt(std::max(0.0, s2 / n - mean * mean)) / mean : 1e30;
            if (cv < best.cv) {
                best.cv = cv;
                best.delta_cm = delta;
                best.radius_cm = radius;
                best.max_residual = mean;  // the plateau, until below
            }
        }
    }
    const double plateau = best.max_residual;
    best.max_residual = 0.0;
    for (std::size_t i = 0; i < distance.size(); ++i) {
        if (!(eps[i] > 0.0))
            continue;
        const double x = distance[i] + best.delta_cm;
        const double s = std::sqrt(x * x + best.radius_cm * best.radius_cm);
        best.max_residual = std::max(
            best.max_residual, std::abs(plateau / (s * (s + x)) / eps[i] - 1.0));
    }
    return best;
}

std::string read_text(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open " + path);
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

} // namespace

int main(int argc, char** argv)
{
    try {
        std::string input, output;
        std::vector<double> energies = {
            15, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 600, 800,
            1000, 1500, 2000, 3000, 5000, 8000, 12000};
        std::vector<double> distances = {5, 7.5, 10, 15, 20, 30, 50, 75, 100,
                                         150, 250, 500, 1000};
        std::uint64_t events = 1000000;
        unsigned threads = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            auto value = [&]() -> std::string {
                if (++i >= argc)
                    throw std::invalid_argument("missing value for " + option);
                return argv[i];
            };
            if (option == "--payload") input = value();
            else if (option == "--output") output = value();
            else if (option == "--energies") energies = number_list(value());
            else if (option == "--distances") distances = number_list(value());
            else if (option == "--events") events = std::stoull(value());
            else if (option == "--threads") threads = static_cast<unsigned>(std::stoul(value()));
            else throw std::invalid_argument("unknown option: " + option);
        }
        if (input.empty() || output.empty() || energies.empty()
            || distances.size() < 3 || events == 0)
            throw std::invalid_argument(
                "--payload and --output are required; at least three distances");
        std::sort(energies.begin(), energies.end());
        std::sort(distances.begin(), distances.end());

        fs::DetectorModelPayload payload =
            fs::DetectorPulseResponseModel::deserialize(read_text(input)).payload();
        if (!payload.geometry)
            throw std::runtime_error("the payload carries no detector geometry");
        const ceelo::GeometryDescriptor& gd = *payload.geometry;
        std::vector<std::unique_ptr<ceelo::Material>> materials;
        for (const auto& spec : gd.materials)
            materials.push_back(std::make_unique<ceelo::Material>(spec.to_material()));
        auto material = [&](int index) -> const ceelo::Material* {
            return materials.at(static_cast<std::size_t>(index)).get();
        };
        const double radius_max = 1.5 * payload.radius_cm;
        const double delta_max = 1.5 * payload.length_cm;
        // The reference point in the crystal frame, as query_position places it.
        const double reference_z = gd.reference_point == ceelo::ReferencePoint::EndcapFront
            ? -gd.endcap_front_offset_cm() : 0.0;

        constexpr double kTargetNoise = 0.003, kMaxNoise = 0.02, kMaxCv = 0.02;
        fs::FiniteDistanceDepth depth;
        std::vector<char> full_ok, detected_ok;
        std::cout << std::setprecision(5)
                  << "# " << payload.detector_name << ": crystal radius "
                  << payload.radius_cm << " cm; " << events << " events per point\n"
                  << "# E_keV  full: delta R (cv, max res)  detected: delta R (cv, max res)\n";
        for (double energy : energies) {
            std::vector<double> full, full_sigma, detected, detected_sigma;
            bool full_valid = true, detected_valid = true;
            for (double distance : distances) {
                auto run = [&](std::uint64_t n) {
                    ceelo::EfficiencyCalculator calc;
                    ceelo::apply_detector_side(calc, gd, material);
                    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, reference_z - distance));
                    return calc.compute(energy, n, threads);
                };
                auto result = run(events);
                const auto noise = [](double value, double sigma) {
                    return value > 0.0 ? sigma / value : 1.0;
                };
                const double worst = std::max(
                    noise(result.full_energy_peak_efficiency, result.fep_uncertainty),
                    noise(result.total_efficiency, result.total_uncertainty));
                if (worst > kTargetNoise && result.total_efficiency > 0.0) {
                    const double scale = std::min(
                        16.0, std::ceil(worst * worst / (kTargetNoise * kTargetNoise)));
                    result = run(static_cast<std::uint64_t>(scale * events));
                }
                full.push_back(result.full_energy_peak_efficiency);
                full_sigma.push_back(result.fep_uncertainty);
                detected.push_back(result.total_efficiency);
                detected_sigma.push_back(result.total_uncertainty);
                full_valid &= noise(result.full_energy_peak_efficiency,
                                    result.fep_uncertainty) <= kMaxNoise;
                detected_valid &= noise(result.total_efficiency,
                                        result.total_uncertainty) <= kMaxNoise;
            }
            const DiskFit f = fit_disk(distances, full, delta_max, radius_max);
            const DiskFit d = fit_disk(distances, detected, delta_max, radius_max);
            full_valid &= f.cv <= kMaxCv;
            detected_valid &= d.cv <= kMaxCv;
            depth.energies_keV.push_back(energy);
            depth.full_energy_cm.push_back(f.delta_cm);
            depth.full_energy_radius_cm.push_back(f.radius_cm);
            depth.detected_cm.push_back(d.delta_cm);
            depth.detected_radius_cm.push_back(d.radius_cm);
            full_ok.push_back(full_valid);
            detected_ok.push_back(detected_valid);
            std::cout << energy << "  " << f.delta_cm << ' ' << f.radius_cm << " ("
                      << f.cv << ", " << f.max_residual
                      << (full_valid ? "" : ", rejected") << ")  " << d.delta_cm
                      << ' ' << d.radius_cm << " (" << d.cv << ", " << d.max_residual
                      << (detected_valid ? "" : ", rejected") << ")\n" << std::flush;
        }
        // Rejected depths take the nearest kept energy's, the higher first.
        auto fill = [](std::vector<double>& value, const std::vector<char>& ok) {
            const std::size_t n = value.size();
            for (std::size_t i = 0; i < n; ++i) {
                if (ok[i])
                    continue;
                for (std::size_t step = 1; step < n; ++step) {
                    if (i + step < n && ok[i + step]) {
                        value[i] = value[i + step];
                        break;
                    }
                    if (step <= i && ok[i - step]) {
                        value[i] = value[i - step];
                        break;
                    }
                }
            }
            return std::any_of(ok.begin(), ok.end(), [](char c) { return c; });
        };
        if (!fill(depth.full_energy_cm, full_ok)
            || !fill(depth.full_energy_radius_cm, full_ok)
            || !fill(depth.detected_cm, detected_ok)
            || !fill(depth.detected_radius_cm, detected_ok))
            throw std::runtime_error("no energy gave a usable depth");
        payload.finite_distance = depth;
        const std::string text = fs::DetectorPulseResponseModel(payload).serialize();
        const std::string temporary = output + ".tmp";
        {
            std::ofstream out(temporary);
            if (!out)
                throw std::runtime_error("cannot write " + temporary);
            out << text;
        }
        if (std::rename(temporary.c_str(), output.c_str()) != 0)
            throw std::runtime_error("cannot replace " + output);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_detector_finite_distance: " << error.what() << '\n';
        return 2;
    }
}
