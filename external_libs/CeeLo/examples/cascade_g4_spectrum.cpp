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

// General engine-vs-GEANT4 cascade-summing spectrum driver.
//
// Sets up a 3"x3" NaI with a point source at -d cm behind a thin Al shell, then:
//   - exports the matched GDML (<nuclide>.gdml) for the G4 harness,
//   - runs compute_cascade (FullRealization) with the angular correlation ON and
//     (optionally) OFF, writing the per-decay summed deposit spectra.
//
// The thin Al shell stops the electron-capture conversion + Auger electrons
// (which the photon-only engine does not emit) so the comparison against G4's
// full atomic-deexcitation decay is a clean gamma + x-ray summing test; the
// shell's small x-ray attenuation is modelled identically in both codes.
//
// Usage: cascade_g4_spectrum <nuclide> [events] [dist_cm] [shield_mm_Al] [emax_keV]

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace ceelo;
using namespace ceelo::cascade_adapter;

namespace {
const SandiaDecay::SandiaDecayDataBase& db() {
  static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
  return database;
}
void write_spectrum(const std::string& path, const CascadeResult& r) {
  FILE* f = std::fopen(path.c_str(), "w");
  std::fprintf(f, "energy_keV,counts_per_decay,uncertainty\n");
  for (std::size_t i = 0; i < r.summed_spectrum.size(); ++i)
    std::fprintf(f, "%zu,%.8g,%.3g\n", i, r.summed_spectrum[i],
                 r.summed_spectrum_uncertainty[i]);
  std::fclose(f);
}
} // namespace

int main(int argc, char** argv) {
  if (argc < 2) { std::fprintf(stderr, "Usage: %s <nuclide> [events] [dist_cm] [shield_mm] [emax]\n", argv[0]); return 1; }
  const std::string nuc = argv[1];
  const uint64_t events = (argc > 2) ? std::strtoull(argv[2], nullptr, 10) : 80000000ULL;
  const double dist = (argc > 3) ? std::atof(argv[3]) : 2.0;
  const double shield_cm = (argc > 4) ? std::atof(argv[4]) / 10.0 : 0.03;  // mm -> cm
  const int emax = (argc > 5) ? std::atoi(argv[5]) : 700;
  // Source geometry: "shell" = point source behind a thin Al spherical shell
  // (default); "alcyl" = source distributed in a small solid Al cylinder of
  // radius/half-height shield_cm. The latter is deadlock-free in G4 (no
  // near-degenerate inner sphere surface) for high-multiplicity low-energy
  // electron emitters like Am-241, while still stopping the EC/IC conversion +
  // Auger electrons and self-attenuating the x-rays identically in both codes.
  const std::string srcmode = (argc > 6) ? argv[6] : "shell";

  static Material nai = make_NaI();
  static Material al = make_Aluminum();

  auto run = [&](bool corr) {
    CascadeOptions opt; opt.angular_correlations = corr;
    const auto casc = build_cascades(db(), nuc, opt);
    EfficiencyCalculator calc;
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    if (srcmode == "alcyl") {
      calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -dist), shield_cm, shield_cm);
      calc.set_source_material(&al);
    } else {
      calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -dist));
      if (shield_cm > 0) calc.add_source_shield(&al, shield_cm);
    }
    if (corr)  // export GDML once (geometry is the same either way)
      calc.export_geant4_gdml(nuc + ".gdml", /*vacuum_world=*/true);
    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.method = CascadeMethod::FullRealization;
    cfg.num_events = events;
    for (int e = 0; e <= emax; ++e) cfg.spectrum_bin_edges.push_back((float)e);
    return calc.compute_cascade(cfg);
  };

  std::printf("%s: NaI 3x3, -%.1f cm, src=%s (%.2f mm Al), %llu events\n",
              nuc.c_str(), dist, srcmode.c_str(), shield_cm * 10.0,
              (unsigned long long)events);
  const CascadeResult on = run(true);
  write_spectrum(nuc + "_eng.csv", on);
  const CascadeResult off = run(false);
  write_spectrum(nuc + "_eng_iso.csv", off);
  std::printf("wrote %s.gdml, %s_eng.csv (W on), %s_eng_iso.csv (W off)\n",
              nuc.c_str(), nuc.c_str(), nuc.c_str());
  return 0;
}
