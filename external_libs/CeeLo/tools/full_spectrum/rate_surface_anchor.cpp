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

// Print the detected-rate surface anchor log odds for (energy, cos theta)
// pairs read from stdin, one pair per line, for the detector a geometry file
// describes (GeometryDescriptor's "CEELO_GEOMETRY 1" text).
//
// The fitter subtracts exactly what the evaluator adds back, so both sides
// MUST use one implementation. Rather than porting the anchor to numpy, the
// Python fitter shells out to this tool, which calls the same
// `rate_surface_anchor_log_odds` the C++ evaluator uses and therefore the
// same cross-section tables. Parity is exact by construction.
//
// --info prints, instead, what the fitter derives from the geometry as one
// JSON object: the passive-layer landmark, the crystal's K edge (the K-line
// element is its heaviest constituent) and the beam envelope.

#include "cross_sections/CrossSectionData.h"
#include "full_spectrum/FullSpectrumModel.h"
#include "geometry/GeometryDescriptor.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    try {
        ceelo::full_spectrum::DetectorRateSurfaceAnchor anchor;
        anchor.enabled = true;
        std::string geometry_file;
        bool info = false;
        for (int index = 1; index < argc; ++index) {
            const std::string flag = argv[index];
            auto value = [&]() -> std::string {
                if (++index >= argc)
                    throw std::invalid_argument("missing value for " + flag);
                return argv[index];
            };
            if (flag == "--geometry")
                geometry_file = value();
            else if (flag == "--blend-kev")
                anchor.blend_keV = std::stod(value());
            else if (flag == "--blend-power")
                anchor.blend_power = std::stod(value());
            else if (flag == "--azimuth-points")
                anchor.azimuth_quadrature_points =
                    static_cast<std::uint32_t>(std::stoul(value()));
            else if (flag == "--info")
                info = true;
            else
                throw std::invalid_argument("unknown option: " + flag);
        }
        if (geometry_file.empty())
            throw std::invalid_argument("--geometry FILE is required");
        std::ifstream in(geometry_file);
        if (!in)
            throw std::runtime_error("cannot read " + geometry_file);
        std::ostringstream text;
        text << in.rdbuf();
        const ceelo::GeometryDescriptor geometry =
            ceelo::GeometryDescriptor::from_text(text.str());
        ceelo::full_spectrum::set_rate_surface_anchor_geometry(anchor, geometry);
        std::cout << std::setprecision(17);
        if (info) {
            int k_element = 0;
            for (const auto& c : anchor.crystal_composition)
                k_element = std::max(k_element, static_cast<int>(c.Z));
            const ceelo::FluorescenceData* fluorescence =
                ceelo::CrossSectionData::instance().fluorescence(k_element);
            std::cout << "{\"passive_landmark_keV\": "
                      << ceelo::full_spectrum::passive_landmark_keV(anchor)
                      << ", \"k_line_element\": " << k_element
                      << ", \"k_edge_keV\": "
                      << (fluorescence ? fluorescence->k_edge_keV : 0.0)
                      << ", \"envelope_radius_cm\": " << anchor.envelope_radius_cm
                      << ", \"envelope_length_cm\": " << anchor.envelope_length_cm
                      << "}\n";
            return 0;
        }
        double energy_keV = 0.0;
        double cos_theta = 0.0;
        while (std::cin >> energy_keV >> cos_theta) {
            std::cout << ceelo::full_spectrum::rate_surface_anchor_log_odds(
                             anchor, energy_keV, cos_theta)
                      << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ceelo_rate_surface_anchor: " << error.what() << '\n';
        return 2;
    }
}
