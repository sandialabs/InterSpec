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

// Synthetic held-out effective-Z fits for the FullSpecApprox forward model.
// This is an executable research harness, not a timing-sensitive unit test.

#include "full_spectrum/FullSpectrumModel.h"

#include <ceres/ceres.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = ceelo::full_spectrum;

namespace {

constexpr double kActivityScale = 1.0e8;

std::string read_text(const std::string& path)
{
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("cannot open detector payload: " + path);
    std::ostringstream text;
    text << stream.rdbuf();
    return text.str();
}

struct Case {
    std::string name;
    std::vector<fs::EmissionLine> lines;
    // activity / 1e8, areal density, effective Z, hydrogen mass fraction
    std::array<double, 4> truth;
    bool activity_only = false;
};

std::vector<double> channels() {
    std::vector<double> result(129);
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = 1600.0 * static_cast<double>(i) / 128.0;
    return result;
}

struct Residual {
    const fs::FastSpectrumModel* model;
    const std::vector<fs::EmissionLine>* lines;
    const std::vector<double>* observed;
    const std::vector<double>* channel_edges;

    template <typename Scalar>
    bool operator()(const Scalar* parameters, Scalar* residuals) const {
        const fs::ShieldParameters<Scalar> shield{
            parameters[1], parameters[2], parameters[3]};
        const fs::ResolutionParameters<Scalar> resolution{
            Scalar(0.21), Scalar(0.0023), Scalar(0.0)};
        const auto prediction = model->predict(
            *lines, parameters[0] * kActivityScale,
            shield, resolution, *channel_edges);
        for (std::size_t i = 0; i < prediction.size(); ++i) {
            residuals[i] = (prediction[i] - (*observed)[i])
                / std::sqrt((*observed)[i] + 1.0);
        }
        return true;
    }
};

struct SolveResult {
    std::array<double, 4> parameters{};
    std::array<double, 16> covariance{};
    double cost = std::numeric_limits<double>::infinity();
    bool converged = false;
    bool covariance_valid = false;
};

SolveResult solve(
    const fs::FastSpectrumModel& model,
    const Case& fit_case,
    const std::vector<double>& observed,
    const std::vector<double>& channel_edges,
    std::array<double, 4> initial,
    int fixed_parameter = -1,
    double fixed_value = 0.0)
{
    if (fixed_parameter >= 0)
        initial[static_cast<std::size_t>(fixed_parameter)] = fixed_value;
    SolveResult answer;
    answer.parameters = initial;
    auto* functor = new Residual{
        &model, &fit_case.lines, &observed, &channel_edges};
    auto* cost =
        new ceres::AutoDiffCostFunction<Residual, ceres::DYNAMIC, 4>(
            functor, static_cast<int>(observed.size()));
    ceres::Problem problem;
    problem.AddResidualBlock(cost, nullptr, answer.parameters.data());
    problem.SetParameterLowerBound(answer.parameters.data(), 0, 1e-6);
    problem.SetParameterUpperBound(answer.parameters.data(), 0, 100.0);
    problem.SetParameterLowerBound(answer.parameters.data(), 1, 0.0);
    problem.SetParameterUpperBound(answer.parameters.data(), 1, 100.0);
    problem.SetParameterLowerBound(answer.parameters.data(), 2, 1.0);
    problem.SetParameterUpperBound(answer.parameters.data(), 2, 92.0);
    problem.SetParameterLowerBound(answer.parameters.data(), 3, 0.0);
    problem.SetParameterUpperBound(answer.parameters.data(), 3, 0.5);

    std::vector<int> fixed;
    if (fit_case.activity_only) {
        fixed = {1, 2, 3};
        answer.parameters[1] = fit_case.truth[1];
        answer.parameters[2] = fit_case.truth[2];
        answer.parameters[3] = fit_case.truth[3];
    } else if (fixed_parameter >= 0) {
        fixed = {fixed_parameter};
    }
    if (!fixed.empty()) {
        problem.SetManifold(
            answer.parameters.data(), new ceres::SubsetManifold(4, fixed));
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 40;
    options.linear_solver_type = ceres::DENSE_QR;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-10;
    options.parameter_tolerance = 1e-10;
    options.num_threads = 1;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    answer.cost = summary.final_cost;
    answer.converged = summary.IsSolutionUsable();

    if (fixed_parameter < 0) {
        ceres::Covariance::Options covariance_options;
        covariance_options.algorithm_type = ceres::DENSE_SVD;
        ceres::Covariance covariance(covariance_options);
        const std::vector<std::pair<const double*, const double*>> blocks{
            {answer.parameters.data(), answer.parameters.data()}};
        if (covariance.Compute(blocks, &problem)) {
            covariance.GetCovarianceBlock(
                answer.parameters.data(), answer.parameters.data(),
                answer.covariance.data());
            answer.covariance_valid = true;
        }
    }
    return answer;
}

std::pair<double, double> profile_interval(
    const fs::FastSpectrumModel& model,
    const Case& fit_case,
    const std::vector<double>& observed,
    const std::vector<double>& channel_edges,
    const SolveResult& best,
    int parameter)
{
    if (fit_case.activity_only && parameter != 0)
        return {fit_case.truth[parameter], fit_case.truth[parameter]};
    double sigma = 0.0;
    if (best.covariance_valid) {
        sigma = std::sqrt(std::max(
            0.0, best.covariance[4 * parameter + parameter]));
    }
    if (!(sigma > 0.0) || !std::isfinite(sigma))
        sigma = parameter == 0 ? 0.25
            : (parameter == 1 ? 0.5
               : (parameter == 2 ? 5.0 : 0.05));
    if (fit_case.activity_only) {
        // With the shield fixed and fixed residual weights, activity is
        // exactly linear. Delta(cost)=0.5 therefore occurs at +/- one
        // covariance sigma and does not require a nested Ceres solve.
        return {
            std::max(1e-6, best.parameters[0] - sigma),
            std::min(100.0, best.parameters[0] + sigma)};
    }
    const std::array<double, 4> lower_bounds{1e-6, 0.0, 1.0, 0.0};
    const std::array<double, 4> upper_bounds{100.0, 100.0, 92.0, 0.5};
    const std::array<double, 4> initial_step_caps{0.01, 0.05, 1.0, 0.02};
    const double threshold = best.cost + 0.5;

    auto boundary = [&](double direction) {
        const double center = best.parameters[parameter];
        const double domain_bound = direction < 0.0
            ? lower_bounds[parameter] : upper_bounds[parameter];
        double inside = center;
        double step = std::min(0.1 * sigma, initial_step_caps[parameter]);
        step = std::max(step, 1e-8);
        double outside = center;
        bool bracketed = false;
        for (int expansion = 0; expansion < 14; ++expansion) {
            const double trial = direction < 0.0
                ? std::max(domain_bound, center - step)
                : std::min(domain_bound, center + step);
            if (trial == inside)
                return inside;
            const SolveResult point = solve(
                model, fit_case, observed, channel_edges,
                best.parameters, parameter, trial);
            // A fixed-coordinate solve can stop on a stringent tolerance
            // without Ceres classifying it as converged. Its finite minimized
            // cost is still the profile-likelihood quantity.
            if (std::isfinite(point.cost) && point.cost <= threshold) {
                inside = trial;
                if (trial == domain_bound)
                    return inside;
                step *= 2.0;
            } else {
                outside = trial;
                bracketed = true;
                break;
            }
        }
        if (!bracketed)
            return inside;
        for (int iteration = 0; iteration < 7; ++iteration) {
            const double trial = 0.5 * (inside + outside);
            const SolveResult point = solve(
                model, fit_case, observed, channel_edges,
                best.parameters, parameter, trial);
            if (std::isfinite(point.cost) && point.cost <= threshold)
                inside = trial;
            else
                outside = trial;
        }
        return inside;
    };
    return {boundary(-1.0), boundary(1.0)};
}

} // namespace

int main(int argc, char** argv) {
    std::string output_path;
    std::string selected_case;
    std::string detector_payload_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc)
            output_path = argv[++i];
        else if (arg == "--case" && i + 1 < argc)
            selected_case = argv[++i];
        else if (arg == "--detector-payload" && i + 1 < argc)
            detector_payload_path = argv[++i];
        else {
            std::cerr << "usage: full_spectrum_ceres_fit "
                         "[--output REPORT.csv] [--case NUCLIDE] "
                         "[--detector-payload PAYLOAD]\n";
            return 2;
        }
    }
    std::ofstream output_file;
    std::ostream* output = &std::cout;
    if (!output_path.empty()) {
        output_file.open(output_path);
        if (!output_file) {
            std::cerr << "cannot open output: " << output_path << '\n';
            return 2;
        }
        output = &output_file;
    }

    fs::FullSpectrumPayload payload;
    payload.source_quadrature_bins = 4;
    payload.detector = detector_payload_path.empty()
        ? fs::coaxial_hpge_payload()
        : fs::DetectorPulseResponseModel::deserialize(
              read_text(detector_payload_path)).payload();
    fs::FastSpectrumModel model(payload);
    const auto channel_edges = channels();
    const fs::ResolutionParameters<double> resolution{
        0.21, 0.0023, 0.0};
    const std::vector<Case> cases{
        {"Ba133",
         {{53.162, 0.021}, {79.614, 0.026}, {80.997, 0.340},
          {276.398, 0.0716}, {302.853, 0.183}, {356.013, 0.6205},
          {383.851, 0.0894}},
         {1.0, 1.7, 47.0, 0.0}, false},
        {"Cs137", {{661.657, 0.851}}, {1.0, 1.2, 28.0, 0.0}, true},
        {"Co60", {{1173.228, 0.9985}, {1332.492, 0.9998}},
         {1.0, 2.3, 31.0, 0.0}, false},
        {"Eu152",
         {{121.782, 0.285}, {244.697, 0.0755}, {344.279, 0.266},
          {778.904, 0.129}, {964.079, 0.146}, {1085.870, 0.102},
          {1112.076, 0.136}, {1408.013, 0.209}},
         {1.0, 3.1, 70.0, 0.0}, false},
        // Water-shielded Co-60: two high-energy Compton-dominated lines —
        // expected to REPORT the (T, Z, f_H) degeneracy, not resolve it.
        {"Co60_water", {{1173.228, 0.9985}, {1332.492, 0.9998}},
         {1.0, 8.0, 8.0, 0.1119}, false},
        // Water-shielded Eu-152: MEASURED degenerate — a light (Z~8)
        // hydrogenous shield attenuates almost purely by Compton even at
        // 122 keV (water PE fraction ~2%), so (T, Z, f_H) collapse to one
        // electron-areal-density parameter; profiles honestly cover truth.
        {"Eu152_water",
         {{121.782, 0.285}, {244.697, 0.0755}, {344.279, 0.266},
          {778.904, 0.129}, {964.079, 0.146}, {1085.870, 0.102},
          {1112.076, 0.136}, {1408.013, 0.209}},
         {1.0, 8.0, 8.0, 0.1119}, false},
        // Hydrogenous mid-Z composite (oil/PE in steel): the 122-344 keV
        // lines see real photoelectric attenuation from Z=26, which anchors
        // effective Z and lets Compton strength identify f_H.
        {"Eu152_fe_hydrogenous",
         {{121.782, 0.285}, {244.697, 0.0755}, {344.279, 0.266},
          {778.904, 0.129}, {964.079, 0.146}, {1085.870, 0.102},
          {1112.076, 0.136}, {1408.013, 0.209}},
         {1.0, 12.0, 26.0, 0.05}, false},
    };

    *output
        << "nuclide,converged,cost,activity,activity_bias,activity_sigma,"
           "activity_profile_low,activity_profile_high,areal_density,"
           "areal_density_bias,areal_density_sigma,areal_profile_low,"
           "areal_profile_high,effective_Z,effective_Z_bias,effective_Z_sigma,"
           "Z_profile_low,Z_profile_high,hydrogen_fraction,"
           "hydrogen_fraction_bias,hydrogen_fraction_sigma,"
           "fH_profile_low,fH_profile_high\n";
    for (const Case& fit_case : cases) {
        if (!selected_case.empty() && fit_case.name != selected_case)
            continue;
        std::cerr << "fitting " << fit_case.name << " with an HPGe response"
                  << std::endl;
        const fs::ShieldParameters<double> truth_shield{
            fit_case.truth[1], fit_case.truth[2], fit_case.truth[3]};
        const auto observed = model.predict(
            fit_case.lines, fit_case.truth[0] * kActivityScale,
            truth_shield, resolution, channel_edges);
        std::vector<std::array<double, 4>> starts{
            {0.4, 0.2, 10.0, 0.0}, {1.8, 5.0, 85.0, 0.3},
            {1.2, 1.0, 60.0, 0.1}};
        SolveResult best;
        int converged_starts = 0;
        for (auto start : starts) {
            const SolveResult candidate = solve(
                model, fit_case, observed, channel_edges, start);
            if (candidate.converged)
                ++converged_starts;
            if (candidate.cost < best.cost)
                best = candidate;
        }
        std::array<std::pair<double, double>, 4> profiles;
        for (int parameter = 0; parameter < 4; ++parameter) {
            profiles[parameter] = profile_interval(
                model, fit_case, observed, channel_edges, best, parameter);
        }
        auto sigma = [&](int parameter) {
            return best.covariance_valid
                ? std::sqrt(std::max(
                      0.0, best.covariance[4 * parameter + parameter]))
                : std::numeric_limits<double>::quiet_NaN();
        };
        *output << std::setprecision(10) << fit_case.name << ','
                  << converged_starts << '/' << starts.size() << ','
                  << best.cost << ','
                  << best.parameters[0] * kActivityScale << ','
                  << (best.parameters[0] - fit_case.truth[0])
                         * kActivityScale << ','
                  << sigma(0) * kActivityScale << ','
                  << profiles[0].first * kActivityScale << ','
                  << profiles[0].second * kActivityScale << ','
                  << best.parameters[1] << ','
                  << best.parameters[1] - fit_case.truth[1] << ','
                  << sigma(1) << ',' << profiles[1].first << ','
                  << profiles[1].second << ','
                  << best.parameters[2] << ','
                  << best.parameters[2] - fit_case.truth[2] << ','
                  << sigma(2) << ',' << profiles[2].first << ','
                  << profiles[2].second << ','
                  << best.parameters[3] << ','
                  << best.parameters[3] - fit_case.truth[3] << ','
                  << sigma(3) << ',' << profiles[3].first << ','
                  << profiles[3].second << '\n';
        output->flush();
    }
    if (!output_path.empty())
        std::cerr << "wrote " << output_path << '\n';
    return 0;
}
