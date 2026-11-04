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

#include <chrono>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

int main() {
    namespace fs = ceelo::full_spectrum;
    fs::FullSpectrumPayload payload;
    // The present research-speed configuration. Sixteen source quadrature
    // bins remains available when detector-folded continuum accuracy warrants
    // the extra cost.
    payload.source_quadrature_bins = 8;
    fs::FastSpectrumModel model(payload);
    std::vector<fs::EmissionLine> lines;
    for (int i = 0; i < 100; ++i)
        lines.push_back({10.0 * std::pow(1200.0, i / 99.0), 0.01});
    std::vector<double> channels(4097);
    for (std::size_t i = 0; i < channels.size(); ++i)
        channels[i] = 12000.0 * static_cast<double>(i) / 4096.0;
    const fs::ShieldParameters<double> shield{1.5, 26.5};
    const fs::ResolutionParameters<double> resolution{4.0, 1.5, 2e-4};

    for (int i = 0; i < 2; ++i)
        (void)model.predict(lines, 1e6, shield, resolution, channels);
    // CPU time, not wall: predict is single-threaded, so on an idle machine
    // the two agree, but a wall-clock gate on a shared machine measures the
    // machine rather than the code. This bench has been observed at 67 s of
    // wall for 4.9 s of CPU while a GEANT4 campaign held the cores.
    const std::clock_t start = std::clock();
    const auto spectrum =
        model.predict(lines, 1e6, shield, resolution, channels);
    const double elapsed = 1.0e3 * static_cast<double>(
        std::clock() - start) / CLOCKS_PER_SEC;
    const double counts =
        std::accumulate(spectrum.begin(), spectrum.end(), 0.0);
    std::cout << std::fixed << std::setprecision(3)
              << "lines=100 channels=4096 cpu_milliseconds=" << elapsed
              << " counts=" << counts
              << " payload_bytes=" << model.serialize().size() << '\n';
    return elapsed < 10000.0 && model.serialize().size() < 2 * 1024 * 1024
        ? 0 : 1;
}
