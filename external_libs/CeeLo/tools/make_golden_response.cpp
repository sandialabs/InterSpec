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

// make_golden_response -- generate a parameterized DetectorResponse XML plus
// a fresh MC probe bank CSV, for use as committed regression fixtures (e.g.
// InterSpec's target/testing/test_data/ceelo_drf/).
//
// All MC is deterministically seeded: rerunning with identical arguments
// reproduces the same fixtures (in distribution across threads; the probe
// truth values regenerate within their recorded MC sigma).
//
// Usage:
//   make_golden_response <preset> <out_dir> [profile] [precision] [n_probes] [seed]
//     preset     nai3x3 | hpge_coax | czt_box | detective_x
//     profile    far-field | general | contact       (default general)
//     precision  per-node MC FEP target              (default 0.003)
//     n_probes   probe-bank size                     (default 60)
//     seed       base seed                           (default 1)
//
// Probe CSV columns match the campaign truth-bank layout (subset):
//   point,E_keV,d_cm,cos_theta,phi_deg,eps_fep,fep_unc,eps_tot,tot_unc,seed

#include "io/DetectorResponse.h"
#include "io/ResponseGenerator.h"
#include "materials/Material.h"

#include "response_presets.h"

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

using namespace ceelo;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::printf(
            "Usage: make_golden_response <preset> <out_dir> [profile] "
            "[precision] [n_probes] [seed] [max_events] [max_cpu_s] "
            "[precision_profile] [certify]\n"
            "  presets: nai3x3 | hpge_coax | czt_box | detective_x\n"
            "  max_events / max_cpu_s raise the per-node caps so a tight "
            "precision is actually reached (high-statistics fixtures);\n"
            "  max_cpu_s is CPU time summed across worker threads (~10x the "
            "old wall-seconds argument).\n"
            "  precision_profile: uniform (default) | relax_mild (D0 graded "
            "map). Committed goldens stay uniform.\n"
            "  certify: pass the literal 'certify' to attach an accuracy "
            "certificate (default off; goldens ship without one).\n");
        return 1;
    }
    const std::string preset = argv[1];
    const std::string out_dir = argv[2];

    GenerationOptions opts;
    opts.detector_name = preset;
    if (argc > 3) {
        if (std::strcmp(argv[3], "far-field") == 0)
            opts.profile = ResponseProfile::FarField;
        else if (std::strcmp(argv[3], "contact") == 0)
            opts.profile = ResponseProfile::Contact;
    }
    if (argc > 4) opts.node_fep_precision = std::atof(argv[4]);
    if (argc > 5) { /* n_probes handled below */ }
    if (argc > 6) opts.base_seed = std::strtoull(argv[6], nullptr, 10);
    if (argc > 7) opts.max_events_per_node = std::strtoull(argv[7], nullptr, 10);
    if (argc > 8) opts.max_cpu_seconds_per_node = std::atof(argv[8]);
    if (argc > 9 && std::strcmp(argv[9], "relax_mild") == 0)
        opts.precision_profile = ceelo::PrecisionProfile::RelaxMild;
    const bool do_certify = (argc > 10 && std::strcmp(argv[10], "certify") == 0);
    const int n_probes = argc > 5 ? std::atoi(argv[5]) : 60;

    opts.progress = [](double frac, const std::string& stage) {
        static int last = -1;
        const int pct = static_cast<int>(100.0 * frac);
        if (pct / 5 != last / 5) {
            std::printf("  [%3d%%] %s\n", pct, stage.c_str());
            std::fflush(stdout);
        }
        last = pct;
    };

    const GeometryDescriptor gd = preset_descriptor(preset);

    std::printf("=== %s: generating response (profile %s, precision %.4f, "
                "~%d MC nodes) ===\n",
                preset.c_str(), argc > 3 ? argv[3] : "general",
                opts.node_fep_precision,
                ResponseGenerator::estimated_node_count(gd, opts));
    std::shared_ptr<DetectorResponse> resp =
        ResponseGenerator::generate(gd, opts);

    if (do_certify) {
        std::printf("=== %s: accuracy certificate (48 probes) ===\n",
                    preset.c_str());
        ResponseGenerator::certify(*resp, gd, opts);
        std::printf("  certificate: fep p95 %.3f%% max %.3f%% (%.1f CPU-s)\n",
                    100.0 * resp->certificate.fep_p95,
                    100.0 * resp->certificate.fep_max,
                    resp->certificate.cpu_seconds);
    }

    const std::string xml_path = out_dir + "/" + preset + "_response.xml";
    {
        std::ofstream f(xml_path);
        if (!f) {
            std::printf("cannot write %s\n", xml_path.c_str());
            return 1;
        }
        f << resp->to_xml_string();
    }
    std::printf("wrote %s (content hash %llu)\n", xml_path.c_str(),
                static_cast<unsigned long long>(resp->content_hash()));

    // Fresh probe bank: Halton offset 5000 -> disjoint from every generation
    // scan (and from any parity split -- campaign pitfall).
    std::printf("=== %s: probe bank (%d points) ===\n", preset.c_str(), n_probes);
    GenerationOptions probe_opts = opts;
    probe_opts.node_fep_precision = std::min(0.002, opts.node_fep_precision);
    probe_opts.progress = opts.progress;
    const std::vector<ProbePoint> bank =
        ResponseGenerator::probe_bank(gd, probe_opts, n_probes, 5000);

    const std::string csv_path = out_dir + "/" + preset + "_probe.csv";
    FILE* f = std::fopen(csv_path.c_str(), "w");
    if (!f) {
        std::printf("cannot write %s\n", csv_path.c_str());
        return 1;
    }
    std::fprintf(f, "# CeeLo golden probe bank: %s (Halton offset 5000)\n",
                 preset.c_str());
    std::fprintf(f, "# response: %s_response.xml (hash %llu)\n", preset.c_str(),
                 static_cast<unsigned long long>(resp->content_hash()));
    std::fprintf(f,
        "point,E_keV,d_cm,cos_theta,phi_deg,eps_fep,fep_unc,eps_tot,tot_unc,seed\n");
    for (size_t i = 0; i < bank.size(); ++i) {
        const ProbePoint& p = bank[i];
        std::fprintf(f, "%zu,%.4f,%.5f,%.6f,%.3f,%.8e,%.3e,%.8e,%.3e,%llu\n",
                     i, p.energy_keV, p.d_cm, p.cos_theta, p.phi_deg, p.eps_fep,
                     p.fep_unc, p.eps_tot, p.tot_unc,
                     static_cast<unsigned long long>(p.seed));
    }
    std::fclose(f);
    std::printf("wrote %s\n", csv_path.c_str());
    return 0;
}
