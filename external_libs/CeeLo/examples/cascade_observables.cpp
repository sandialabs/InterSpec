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

// Cascade-summing MC observable producer for profiling/compare_validation.py.
//
// Reads the GEANT4 reference table (tests/data/geant4_reference/cascade_summing_multi.csv),
// runs compute_cascade on the SAME "alcyl" geometry the reference was generated in, and
// writes the engine's per-decay observable areas to our_cascade_multi.csv. The dashboard
// then gates engine/G4 ratios against the reference's bands, exactly as it does the
// per-energy efficiency configs.
//
//   full rows: FullRealization summed-spectrum area over [lo,hi].
//   cond rows: Conditional per-decay area of the peak line = A_FR(peak) * k_cond/k_full,
//              anchoring the Conditional estimator to FullRealization's per-decay
//              normalization while injecting only its summing-physics ratio.
//
// Requires -DCEELO_WITH_SANDIADECAY=ON. SANDIA_DECAY_XML_PATH and CASCADE_SUMMING_REF are
// provided by the build.
//
// Usage: cascade_observables [--ref <csv>] [--out <csv>]
//                            [--events-full N] [--events-cond N]

#include "cascade_ref_common.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_ref;

int main(int argc, char** argv) {
    std::string ref_path = CASCADE_SUMMING_REF;
    std::string out_path = "our_cascade_multi.csv";
    uint64_t events_full = 8000000;
    uint64_t events_cond = 4000000;

    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--ref") && i + 1 < argc) ref_path = argv[++i];
        else if (!std::strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
        else if (!std::strcmp(argv[i], "--events-full") && i + 1 < argc)
            events_full = std::strtoull(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--events-cond") && i + 1 < argc)
            events_cond = std::strtoull(argv[++i], nullptr, 10);
        else {
            std::fprintf(stderr, "unknown/again arg: %s\n", argv[i]);
            return 2;
        }
    }

    const std::vector<Obs> rows = load_reference(ref_path);
    if (rows.empty()) {
        std::fprintf(stderr, "no reference rows read from %s\n", ref_path.c_str());
        return 1;
    }

    // Preserve the reference file's nuclide order (nicer diffs than alphabetical).
    std::vector<std::string> order;
    std::map<std::string, std::vector<Obs>> by_nuc;
    for (const Obs& o : rows) {
        if (by_nuc.find(o.nuclide) == by_nuc.end()) order.push_back(o.nuclide);
        by_nuc[o.nuclide].push_back(o);
    }

    std::FILE* out = std::fopen(out_path.c_str(), "w");
    if (!out) {
        std::fprintf(stderr, "cannot open %s for writing\n", out_path.c_str());
        return 1;
    }
    std::fprintf(out, "# CeeLo cascade-summing MC results (per-decay areas).\n");
    std::fprintf(out, "# Producer: cascade_observables (FullRealization %llu, Conditional %llu events).\n",
                 (unsigned long long)events_full, (unsigned long long)events_cond);
    std::fprintf(out, "nuclide,estimator,name,area,area_unc\n");

    for (const std::string& nuc : order) {
        const std::vector<Obs>& these = by_nuc[nuc];
        const double source_cm = these.front().source_cm;
        const int emax = these.front().emax;

        std::vector<double> cond_peaks;
        for (const Obs& o : these)
            if (o.estimator == "cond") cond_peaks.push_back(o.peak_keV);

        std::fprintf(stderr, "[%s] FullRealization (emax=%d, src=%.2f cm)...\n",
                     nuc.c_str(), emax, source_cm);
        const FullResult fr =
            run_fullrealization(nuc, emax, events_full, source_cm, cond_peaks);

        std::vector<double> kc, kcu;
        if (!cond_peaks.empty()) {
            std::fprintf(stderr, "[%s] Conditional (%zu peaks)...\n",
                         nuc.c_str(), cond_peaks.size());
            run_conditional(nuc, events_cond, source_cm, cond_peaks, kc, kcu);
        }

        std::size_t ci = 0;
        for (const Obs& o : these) {
            double val = 0.0, unc = 0.0;
            if (o.estimator == "full") {
                val = area(fr.spectrum, o.lo, o.hi);
                unc = area_unc(fr.spectrum_unc, o.lo, o.hi);
            } else {  // cond: A_FR(window) * k_cond / k_full
                const double afr = area(fr.spectrum, o.lo, o.hi);
                const double afru = area_unc(fr.spectrum_unc, o.lo, o.hi);
                const double kf = fr.k[ci], kfu = fr.k_unc[ci];
                const double kcc = kc[ci], kccu = kcu[ci];
                ++ci;
                if (kf > 0.0 && kcc > 0.0 && afr > 0.0) {
                    val = afr * kcc / kf;
                    double rel = (afru / afr) * (afru / afr)
                               + (kccu / kcc) * (kccu / kcc)
                               + (kfu / kf) * (kfu / kf);
                    unc = val * std::sqrt(rel);
                }
            }
            std::fprintf(out, "%s,%s,%s,%.6e,%.6e\n",
                         o.nuclide.c_str(), o.estimator.c_str(), o.name.c_str(),
                         val, unc);
        }
    }

    std::fclose(out);
    std::fprintf(stderr, "wrote %s\n", out_path.c_str());
    return 0;
}
