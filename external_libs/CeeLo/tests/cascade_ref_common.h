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

// Shared helpers for the cascade-summing GEANT4 gate: the "alcyl" compute_cascade
// run, spectrum-window integration, and the reference-CSV reader. Included by BOTH
// tests/test_cascade_summing.cpp (the ctest regression gate) and
// examples/cascade_observables.cpp (the compare_validation.py MC producer) so the
// geometry, the observable definitions, and the parsed reference share ONE
// implementation and cannot drift.
//
// The includer must provide two compile definitions:
//   SANDIA_DECAY_XML_PATH  - path to sandia.decay.xml (for the adapter)
//   CASCADE_SUMMING_REF    - path to tests/data/geant4_reference/cascade_summing_multi.csv
#ifndef CEELO_CASCADE_REF_COMMON_H
#define CEELO_CASCADE_REF_COMMON_H

#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include "cascade/SandiaDecayCascade.h"
#include "SandiaDecay.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace ceelo {
namespace cascade_ref {

// The adapter types (the SandiaDecay -> DecayCascade builder) live in cascade_adapter.
using cascade_adapter::CascadeOptions;
using cascade_adapter::build_cascades;

// One summing observable parsed from the reference CSV.
//   estimator = "full" -> FullRealization summed-spectrum area over [lo,hi]
//   estimator = "cond" -> Conditional per-decay area of the line at peak_keV
struct Obs {
    std::string nuclide;
    std::string estimator;  // "full" | "cond"
    std::string name;
    int lo = 0, hi = 0;
    double peak_keV = 0.0;   // cond rows only
    double g4 = 0.0, g4_sig = 0.0;
    double rlo = 0.0, rhi = 0.0;
    double source_cm = 0.03;
    int emax = 0;
};

// All observables for one nuclide, split by estimator (they share source_cm/emax).
struct NuclideRefs {
    std::string nuclide;
    double source_cm = 0.03;
    int emax = 0;
    std::vector<Obs> full;
    std::vector<Obs> cond;
};

// FullRealization (spectrum) + optional per-peak summing factors, in ONE run.
struct FullResult {
    std::vector<float> spectrum;       // per-decay counts, 1 keV bins (lower edge = index)
    std::vector<float> spectrum_unc;   // per-bin 1sigma
    std::vector<double> k;             // summing factor per requested peak
    std::vector<double> k_unc;
};

inline const SandiaDecay::SandiaDecayDataBase& db() {
    static SandiaDecay::SandiaDecayDataBase database(SANDIA_DECAY_XML_PATH);
    return database;
}

// Sum per-decay counts (or uncertainties in quadrature) over the window [lo,hi] keV.
inline double area(const std::vector<float>& y, int lo, int hi) {
    double s = 0.0;
    for (int e = lo; e <= hi && e < static_cast<int>(y.size()); ++e) s += y[e];
    return s;
}
inline double area_unc(const std::vector<float>& u, int lo, int hi) {
    double s = 0.0;
    for (int e = lo; e <= hi && e < static_cast<int>(u.size()); ++e) s += u[e] * u[e];
    return std::sqrt(s);
}

// The alcyl geometry, shared by every run below: 3"x3" NaI, source distributed in a
// solid Al cylinder (r = halfz = source_cm) 2 cm behind the front face.
inline void configure_alcyl(EfficiencyCalculator& calc, double source_cm) {
    static Material nai = make_NaI();
    static Material al = make_Aluminum();
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
    calc.set_cylindrical_source(Eigen::Vector3d(0.0, 0.0, -2.0), source_cm, source_cm);
    calc.set_source_material(&al);
}

// FullRealization run: summed spectrum over [0,emax] plus the summing factor at each
// requested peak (angular correlations + vacancy model on by default).
inline FullResult run_fullrealization(const std::string& nuc, int emax_keV,
                                      uint64_t num_events, double source_cm,
                                      const std::vector<double>& peaks_keV = {}) {
    CascadeOptions opt;
    const auto casc = build_cascades(db(), nuc, opt);

    EfficiencyCalculator calc;
    configure_alcyl(calc, source_cm);

    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.method = CascadeMethod::FullRealization;
    cfg.num_events = num_events;
    cfg.num_threads = 4;  // fixed => reproducible event split / seeds
    for (int e = 0; e <= emax_keV; ++e)
        cfg.spectrum_bin_edges.push_back(static_cast<float>(e));
    for (double p : peaks_keV) cfg.peaks.push_back({p, 1.5});

    const CascadeResult res = calc.compute_cascade(cfg);
    FullResult out;
    out.spectrum = res.summed_spectrum;
    out.spectrum_unc = res.summed_spectrum_uncertainty;
    for (std::size_t i = 0; i < peaks_keV.size(); ++i) {
        const double k = (i < res.peaks.size() && res.peaks[i].found)
                             ? res.peaks[i].summing_factor : 0.0;
        const double ku = (i < res.peaks.size() && res.peaks[i].found)
                              ? res.peaks[i].summing_factor_unc : 0.0;
        out.k.push_back(k);
        out.k_unc.push_back(ku);
    }
    return out;
}

// Conditional run in the SAME alcyl geometry: per-peak summing factor + 1sigma.
inline void run_conditional(const std::string& nuc, uint64_t num_events,
                            double source_cm, const std::vector<double>& peaks_keV,
                            std::vector<double>& k_out, std::vector<double>& k_unc_out) {
    CascadeOptions opt;
    const auto casc = build_cascades(db(), nuc, opt);

    EfficiencyCalculator calc;
    configure_alcyl(calc, source_cm);

    CascadeConfig cfg;
    cfg.cascades = casc;
    cfg.method = CascadeMethod::Conditional;
    cfg.num_events = num_events;
    cfg.num_threads = 4;
    for (double p : peaks_keV) cfg.peaks.push_back({p, 1.5});

    const CascadeResult res = calc.compute_cascade(cfg);
    k_out.clear();
    k_unc_out.clear();
    for (std::size_t i = 0; i < peaks_keV.size(); ++i) {
        k_out.push_back((i < res.peaks.size() && res.peaks[i].found)
                            ? res.peaks[i].summing_factor : 0.0);
        k_unc_out.push_back((i < res.peaks.size() && res.peaks[i].found)
                                ? res.peaks[i].summing_factor_unc : 0.0);
    }
}

// Parse the reference CSV into flat rows (skips '#' comment and header lines).
inline std::vector<Obs> load_reference(const std::string& path) {
    std::vector<Obs> rows;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("nuclide", 0) == 0) continue;  // header
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> c;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (c.size() < 12) continue;
        Obs o;
        o.nuclide = c[0];
        o.estimator = c[1];
        o.name = c[2];
        o.lo = std::atoi(c[3].c_str());
        o.hi = std::atoi(c[4].c_str());
        o.peak_keV = std::atof(c[5].c_str());
        o.g4 = std::atof(c[6].c_str());
        o.g4_sig = std::atof(c[7].c_str());
        o.rlo = std::atof(c[8].c_str());
        o.rhi = std::atof(c[9].c_str());
        o.source_cm = std::atof(c[10].c_str());
        o.emax = std::atoi(c[11].c_str());
        rows.push_back(o);
    }
    return rows;
}

// Group flat rows by nuclide (preserving source_cm/emax and the full/cond split).
inline std::map<std::string, NuclideRefs> group_reference(const std::vector<Obs>& rows) {
    std::map<std::string, NuclideRefs> by_nuc;
    for (const Obs& o : rows) {
        NuclideRefs& n = by_nuc[o.nuclide];
        n.nuclide = o.nuclide;
        n.source_cm = o.source_cm;
        n.emax = o.emax;
        if (o.estimator == "cond") n.cond.push_back(o);
        else n.full.push_back(o);
    }
    return by_nuc;
}

}  // namespace cascade_ref
}  // namespace ceelo

#endif  // CEELO_CASCADE_REF_COMMON_H
