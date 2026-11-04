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

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ceelo::full_spectrum::electron_campaign {

struct Tally {
    double sum_w = 0.0;
    double sum_w2 = 0.0;
    void add(double weight) {
        sum_w += weight;
        sum_w2 += weight * weight;
    }
};

class Spectra {
public:
    Spectra(int bins, double maximum_energy_keV)
        : bins_(bins), maximum_energy_keV_(maximum_energy_keV) {
        if (bins < 1 || !(maximum_energy_keV > 0.0))
            throw std::invalid_argument("invalid spectrum binning");
    }

    void tally(const std::string& particle, const std::string& component,
               double energy_keV, double weight) {
        if (!(energy_keV >= 0.0) || energy_keV >= maximum_energy_keV_
            || !(weight > 0.0))
            return;
        auto& values = values_[{particle, component}];
        if (values.empty())
            values.resize(static_cast<std::size_t>(bins_));
        const int bin = std::min(
            bins_ - 1,
            static_cast<int>(
                energy_keV / maximum_energy_keV_ * bins_));
        values[static_cast<std::size_t>(bin)].add(weight);
    }

    void write(const std::string& path, const std::string& campaign,
               const std::string& source, double primary_energy_keV,
               int Z, double areal_density_g_cm2,
               std::uint64_t normalization_events, std::uint64_t seed,
               const std::string& physics) const {
        std::ofstream stream(path);
        if (!stream)
            throw std::runtime_error("cannot open " + path);
        stream << std::setprecision(12)
            << "campaign,source,primary_energy_keV,Z,"
               "areal_density_g_cm2,particle,component,bin_lo_keV,"
               "bin_hi_keV,sum_w,sum_w2,normalization_events,seed,physics\n";
        for (const auto& entry : values_) {
            for (int bin = 0; bin < bins_; ++bin) {
                const double lo =
                    maximum_energy_keV_ * bin / bins_;
                const double hi =
                    maximum_energy_keV_ * (bin + 1) / bins_;
                const Tally& tally =
                    entry.second[static_cast<std::size_t>(bin)];
                stream << campaign << ',' << source << ','
                       << primary_energy_keV << ',' << Z << ','
                       << areal_density_g_cm2 << ','
                       << entry.first.first << ',' << entry.first.second << ','
                       << lo << ',' << hi << ',' << tally.sum_w << ','
                       << tally.sum_w2 << ',' << normalization_events << ','
                       << seed << ',' << physics << '\n';
            }
        }
    }

private:
    int bins_;
    double maximum_energy_keV_;
    std::map<std::pair<std::string, std::string>, std::vector<Tally>> values_;
};

} // namespace ceelo::full_spectrum::electron_campaign
