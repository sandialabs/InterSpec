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

/// @file env_fit_payload_eval.cpp
/// @brief Dependency-free cross-language evaluator for CEEFIT runtime payloads.
///
/// The protocol is deliberately a small, strictly parsed token stream instead
/// of JSON so the production C++ evaluator acquires no serialization
/// dependency.  Validation tooling may load candidate payloads explicitly; the
/// ordinary runtime loader remains fail-closed.

#include "io/EnvFitPayload.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

constexpr std::size_t kMaximumQueries = 4096U;
constexpr std::size_t kMaximumEdges = 4097U;
constexpr std::size_t kMaximumCells = 1000000U;

void expect(std::istream& input, const char* expected)
{
    std::string token;
    if (!(input >> token) || token != expected)
        throw std::runtime_error(std::string("expected protocol token ") + expected);
}

std::size_t read_count(std::istream& input, const char* name,
                       std::size_t minimum, std::size_t maximum)
{
    unsigned long long value = 0;
    if (!(input >> value) || value < minimum || value > maximum)
        throw std::runtime_error(std::string("invalid ") + name);
    return static_cast<std::size_t>(value);
}

double read_finite(std::istream& input, const char* name)
{
    double value = 0.0;
    if (!(input >> value) || !std::isfinite(value))
        throw std::runtime_error(std::string("invalid finite value for ") + name);
    return value;
}

std::vector<double> read_edges(std::istream& input, const char* token,
                               std::size_t count)
{
    expect(input, token);
    std::vector<double> values(count);
    for (double& value : values) value = read_finite(input, token);
    for (std::size_t i = 1; i < values.size(); ++i)
        if (!(values[i] > values[i - 1U]))
            throw std::runtime_error(std::string(token) + " edges are not increasing");
    return values;
}

void write_vector(std::ostream& output, const char* token,
                  const std::vector<double>& values)
{
    output << token << ' ' << values.size();
    for (double value : values) {
        if (!std::isfinite(value))
            throw std::runtime_error(std::string("evaluator produced non-finite ") + token);
        output << ' ' << value;
    }
    output << '\n';
}

void write_named_vectors(std::ostream& output, const char* token,
                         const std::string& name,
                         const std::vector<double>& s_hat,
                         const std::vector<double>& incident)
{
    if (name.empty() || name.find_first_of(" \t\r\n") != std::string::npos)
        throw std::runtime_error("payload contains a protocol-unsafe model name");
    output << token << ' ' << name << ' ' << s_hat.size();
    for (double value : s_hat) {
        if (!std::isfinite(value)) throw std::runtime_error("non-finite fitted mass");
        output << ' ' << value;
    }
    output << ' ' << incident.size();
    for (double value : incident) {
        if (!std::isfinite(value)) throw std::runtime_error("non-finite incident mass");
        output << ' ' << value;
    }
    output << '\n';
}

void evaluate_requests(const std::string& payload_path,
                       std::istream& input, std::ostream& output)
{
    EnvFitPayloadEvaluator evaluator;
    std::string error;
    if (!evaluator.load(payload_path, &error, true))
        throw std::runtime_error("cannot load fitted candidate: " + error);

    expect(input, "CEEFIT-PARITY-REQUEST-v1");
    expect(input, "N");
    const std::size_t query_count = read_count(
        input, "query count", 1U, kMaximumQueries);

    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "CEEFIT-PARITY-RESULT-v1\nMODEL "
           << env_fit_payload_model_content_sha256(evaluator.payload())
           << "\nN " << query_count << '\n';

    for (std::size_t index = 0; index < query_count; ++index) {
        expect(input, "Q");
        const std::size_t supplied_index = read_count(
            input, "query index", 0U, kMaximumQueries - 1U);
        if (supplied_index != index)
            throw std::runtime_error("query indices are not contiguous and zero based");
        EnvFitRuntimeQuery query;
        query.e0_keV = read_finite(input, "source energy");
        query.distance_m = read_finite(input, "distance");
        query.source_height_m = read_finite(input, "source height");
        query.detector_height_m = read_finite(input, "detector height");
        query.direct_fluence_per_cm2 = read_finite(input, "direct fluence");
        const std::size_t energy_count = read_count(
            input, "energy-edge count", 2U, kMaximumEdges);
        const std::size_t mu_count = read_count(
            input, "mu-edge count", 2U, kMaximumEdges);
        if ((energy_count - 1U) > kMaximumCells / (mu_count - 1U))
            throw std::runtime_error("requested evaluation has too many cells");
        const std::vector<double> energy_edges = read_edges(input, "E", energy_count);
        const std::vector<double> mu_edges = read_edges(input, "M", mu_count);

        EnvFitEvaluation result;
        if (!evaluator.evaluate(query, energy_edges, mu_edges, result, &error))
            throw std::runtime_error("query " + std::to_string(index) + ": " + error);
        output << "Q " << index << ' ' << result.n_energy << ' ' << result.n_mu
               << ' ' << result.components.size() << ' ' << result.lines.size()
               << ' ' << result.direct_fluence_per_cm2 << '\n';
        write_vector(output, "TOTAL", result.continuum_s_hat_cell_masses);
        write_vector(output, "INCIDENT_TOTAL",
                     result.incident_continuum_cell_masses_per_cm2);
        for (const auto& component : result.components)
            write_named_vectors(output, "COMP", component.name,
                                component.s_hat_cell_masses,
                                component.incident_cell_masses_per_cm2);
        for (const auto& line : result.lines) {
            output << "LINE " << line.name << ' ' << line.energy_keV << ' '
                   << line.s_hat_mu_bin_masses.size();
            for (double value : line.s_hat_mu_bin_masses) output << ' ' << value;
            output << ' ' << line.incident_mu_bin_masses_per_cm2.size();
            for (double value : line.incident_mu_bin_masses_per_cm2)
                output << ' ' << value;
            output << '\n';
        }
        output << "ENDQ\n";
    }
    expect(input, "END");
    std::string trailing;
    if (input >> trailing) throw std::runtime_error("trailing request data");
    output << "END\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: env_fit_payload_eval PAYLOAD [REQUEST [RESULT]]\n";
        return 2;
    }
    try {
        if (argc == 4) {
            const auto result = std::filesystem::weakly_canonical(argv[3]);
            if (result == std::filesystem::weakly_canonical(argv[1]) ||
                result == std::filesystem::weakly_canonical(argv[2]))
                throw std::runtime_error(
                    "result file must not overwrite the payload or request");
        }
        std::unique_ptr<std::ifstream> input_file;
        std::unique_ptr<std::ofstream> output_file;
        std::istream* input = &std::cin;
        std::ostream* output = &std::cout;
        if (argc >= 3) {
            input_file = std::make_unique<std::ifstream>(argv[2]);
            if (!*input_file) throw std::runtime_error("cannot open request file");
            input = input_file.get();
        }
        if (argc == 4) {
            output_file = std::make_unique<std::ofstream>(argv[3]);
            if (!*output_file) throw std::runtime_error("cannot open result file");
            output = output_file.get();
        }
        evaluate_requests(argv[1], *input, *output);
        if (!*output) throw std::runtime_error("cannot write complete result");
    } catch (const std::exception& exception) {
        std::cerr << "env_fit_payload_eval: " << exception.what() << '\n';
        return 1;
    }
    return 0;
}
