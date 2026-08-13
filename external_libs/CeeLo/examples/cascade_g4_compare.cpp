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

// Engine side of the W(theta) vs GEANT4 cross-check.
//
// Reproduces the beta-absorbing shielded Co-60 geometry of export_cascade_g4
// (NaI 3"x3", point source at -2 cm, 0.2 cm Fe shell) and runs compute_cascade
// (FullRealization) twice -- with the gamma-gamma angular correlation enabled
// (a2/a4 from the enriched SandiaDecay data) and forced isotropic -- producing
// the per-decay summed deposit spectrum each way. Extracts the 1173 / 1332 keV
// photopeak and 2505 keV sum-peak areas and the within-code W-on/W-off ratios,
// which cancel geometry/efficiency and isolate the angular-correlation effect for
// comparison against GEANT4 run with/without --correlated-gamma.
//
// Writes eng_corr.csv / eng_iso.csv (per-decay deposit spectra) for direct
// spectrum comparison with the G4 *_histogram.csv files.
//
// Usage: cascade_g4_compare [events]   (default 40,000,000)

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_adapter;

namespace {
const SandiaDecay::SandiaDecayDataBase& db() {
  static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
  return database;
}

struct Area { double a = 0, u = 0; };

// Sum counts/decay over [e-w, e+w], combining per-bin uncertainties in quadrature.
Area peak_area(const CascadeResult& r, double e, double w) {
  Area out;
  const std::size_t n = r.summed_spectrum.size();
  for (std::size_t i = 0; i < n; ++i) {
    const double elo = i;  // 1 keV bins from 0
    if (elo >= e - w && elo <= e + w) {
      out.a += r.summed_spectrum[i];
      out.u += static_cast<double>(r.summed_spectrum_uncertainty[i]) *
               r.summed_spectrum_uncertainty[i];
    }
  }
  out.u = std::sqrt(out.u);
  return out;
}

void write_spectrum(const char* path, const CascadeResult& r) {
  FILE* f = std::fopen(path, "w");
  std::fprintf(f, "energy_keV,counts_per_decay,uncertainty\n");
  for (std::size_t i = 0; i < r.summed_spectrum.size(); ++i)
    std::fprintf(f, "%zu,%.8g,%.3g\n", i, r.summed_spectrum[i],
                 r.summed_spectrum_uncertainty[i]);
  std::fclose(f);
}

CascadeResult run(const std::vector<DecayCascade>& casc, uint64_t events) {
  Material nai = make_NaI();
  Material fe = make_Iron();
  EfficiencyCalculator calc;
  calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
  calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -2.0));
  calc.add_source_shield(&fe, 0.2);
  CascadeConfig cfg;
  cfg.cascades = casc;
  cfg.method = CascadeMethod::FullRealization;
  cfg.num_events = events;
  cfg.peaks = {PeakWindow{1173.2, 1.5}, PeakWindow{1332.5, 1.5}};
  for (int e = 0; e <= 3000; ++e) cfg.spectrum_bin_edges.push_back((float)e);
  return calc.compute_cascade(cfg);
}

void report(const char* label, double e, double w,
            const CascadeResult& on, const CascadeResult& off) {
  const Area a = peak_area(on, e, w), b = peak_area(off, e, w);
  const double ratio = b.a > 0 ? a.a / b.a : 0;
  const double ru = (a.a > 0 && b.a > 0)
      ? ratio * std::sqrt(std::pow(a.u / a.a, 2) + std::pow(b.u / b.a, 2)) : 0;
  std::printf("%-12s  W=%.6e+-%.1e  iso=%.6e+-%.1e  W/iso=%.5f+-%.5f\n",
              label, a.a, a.u, b.a, b.u, ratio, ru);
}
} // namespace

int main(int argc, char** argv) {
  const uint64_t events = (argc > 1) ? std::strtoull(argv[1], nullptr, 10) : 40000000ULL;
  CascadeOptions on_opts;                          // angular_correlations=true
  CascadeOptions off_opts; off_opts.angular_correlations = false;
  const auto casc_on = build_cascades(db(), "Co60", on_opts);
  const auto casc_off = build_cascades(db(), "Co60", off_opts);

  std::printf("Co-60 shielded (NaI 3x3, pt -2cm, 0.2cm Fe), %llu events, FullRealization\n",
              (unsigned long long)events);
  const CascadeResult on = run(casc_on, events);
  const CascadeResult off = run(casc_off, events);

  write_spectrum("eng_corr.csv", on);
  write_spectrum("eng_iso.csv", off);

  std::printf("%-12s  %-26s %-26s %s\n", "peak", "W(theta) [cts/decay]",
              "isotropic [cts/decay]", "ratio");
  report("1173 phot", 1173.2, 2.5, on, off);
  report("1332 phot", 1332.5, 2.5, on, off);
  report("2505 sum",  2505.0, 3.0, on, off);
  return 0;
}
