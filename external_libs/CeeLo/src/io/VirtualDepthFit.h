#ifndef CEELO_IO_VIRTUAL_DEPTH_FIT_H
#define CEELO_IO_VIRTUAL_DEPTH_FIT_H
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

// Virtual interaction depth delta(E) fitter + exporter.
//
// Problem: a face-referenced "intrinsic" efficiency drifts with distance because
// the point from which the crystal effectively radiates moves deeper into the
// crystal as energy rises. InterSpec models this with a single energy-INDEPENDENT
// `setback`; this utility measures the energy-DEPENDENT virtual depth delta(E) so
// the source integral can use solid-angle-to-(d + setback + delta(E)).
//
// Method (generalized from the efficiency-grid study's fit_vpd): for each
// energy, find the depth delta that makes eps(d) / f_Omega(d + delta, R) constant
// across distance, where f_Omega is the disk solid-angle fraction of the crystal
// face. delta(E) is then fit to an exp-of-log power series so the coefficients
// slot into InterSpec's efficiency-coefficient convention (energy in MeV).
//
// This is a STANDALONE orchestration utility: fit_virtual_depth() drives
// EfficiencyCalculator::compute(); nothing here changes MC transport.

#include "geometry/Geometry.h"   // DetectorShape
#include "materials/Material.h"  // Material

#include <cstdint>
#include <string>
#include <vector>

namespace ceelo {

/// Self-describing detector identity carried in the exported delta(E) file.
struct DetectorDescriptor {
    std::string name;            ///< e.g. "NaI 3x3"
    std::string material_name;   ///< Material::name(), e.g. "NaI"
    DetectorShape shape = DetectorShape::Cylinder;
    double crystal_radius_cm = 0.0;    ///< R used in f_Omega (cylinder radius, or R_eff for a box)
    double crystal_diameter_cm = 0.0;  ///< 2 * crystal_radius_cm (convenience for InterSpec)
    double crystal_length_cm = 0.0;    ///< full crystal length along the axis
    double half_x_cm = 0.0;            ///< box only (0 for cylinder)
    double half_y_cm = 0.0;            ///< box only (0 for cylinder)
};

/// Build a descriptor from a detector spec. `dimensions` follow set_detector():
/// Cylinder = {radius, half_length}; Box = {half_x, half_y, length}. For a box,
/// crystal_radius_cm is the area-equivalent disk radius R_eff = 2*sqrt(hx*hy/pi)
/// (an approximation; the CV fit absorbs residual face-shape error).
DetectorDescriptor make_descriptor(std::string name, DetectorShape shape,
                                   const Material& material,
                                   const std::vector<double>& dimensions);

/// Configuration for a full delta(E) fit.
struct VpdFitConfig {
    std::vector<double> energies_keV;        ///< energies at which to fit delta
    std::vector<double> distance_ladder_cm;  ///< distance nodes for the per-energy CV fit
    bool include_near_field = true;          ///< if false, the caller's ladder is used as-is
                                             ///<   (Phase-4 fit-far stress); informational flag
    double delta_max_cm = 6.0;               ///< CV grid-search upper bound
    double delta_step_cm = 0.02;             ///< CV grid-search step
    uint64_t events_per_point = 400000;      ///< MC events per (energy, distance)
    unsigned num_threads = 0;                ///< 0 = auto
    int poly_order = 3;                      ///< # exp-of-log terms (p0..p_{order-1})
};

/// Per-energy fit result + the raw distance scan that produced it.
struct VpdPoint {
    double energy_keV = 0.0;
    double delta_cm = 0.0;        ///< fitted virtual depth
    double intrinsic_mean = 0.0;  ///< mean of eps / f_Omega(d + delta) (the plateau)
    double cv_residual = 0.0;     ///< coefficient of variation at best delta (fit quality)
    double max_abs_res = 0.0;     ///< max |pred - eps| / eps over the ladder
    std::vector<double> dist_cm;  ///< raw distances
    std::vector<double> eps;      ///< raw FEP efficiencies
    std::vector<double> eps_unc;  ///< raw MC 1-sigma uncertainties
};

/// Full delta(E) fit for one detector.
struct VpdFit {
    DetectorDescriptor detector;
    std::vector<double> coeffs;             ///< p0..p_{order-1}; delta = exp(sum p_k u^k), u = ln(E_MeV)
    std::vector<VpdPoint> points;           ///< raw per-energy fits (sorted by energy)
    std::vector<double> log_fit_residuals;  ///< ln(delta_k) - model(E_k), per energy
    double valid_e_min_keV = 0.0;
    double valid_e_max_keV = 0.0;
    double rms_log_residual = 0.0;          ///< overall exp-of-log fit quality

    /// Evaluate the fitted model at an arbitrary energy (energy clamped to the
    /// [valid_e_min, valid_e_max] range before evaluation).
    double delta_at(double energy_keV) const;
};

// ---- Pure (no-MC) building blocks ----

/// Per-energy delta fit by CV minimization over a delta grid. Mirrors
/// the efficiency-grid study's fit_vpd, generalized to arbitrary R and
/// distance ladder, using disk_solid_angle_fraction(d + delta, R). `eps_unc` is
/// carried into the returned VpdPoint (not used by the objective).
VpdPoint fit_delta_one_energy(double R_cm,
                              const std::vector<double>& dist_cm,
                              const std::vector<double>& eps,
                              const std::vector<double>& eps_unc,
                              double delta_max_cm = 6.0,
                              double delta_step_cm = 0.02);

/// Least-squares fit of ln(delta) on powers of u = ln(E_MeV). Returns coeffs
/// p0..p_{order-1}. delta is floored at `delta_floor_cm` before taking the log
/// to guard degenerate / near-zero low-energy fits.
std::vector<double> fit_exp_log_coeffs(const std::vector<double>& energies_keV,
                                       const std::vector<double>& delta_cm,
                                       int order,
                                       double delta_floor_cm = 1e-3);

/// Evaluate delta(E) = exp(sum_k coeffs[k] * u^k), u = ln(E_MeV). E in keV.
double eval_exp_log(const std::vector<double>& coeffs, double energy_keV);

// ---- Full pipeline (drives MC) ----

/// Configure an EfficiencyCalculator for `descriptor`/`material`/`dimensions`,
/// run compute() over the (energy x distance) grid, fit delta per energy, and
/// fit the exp-of-log coefficients. This is the only function here that runs MC.
VpdFit fit_virtual_depth(const DetectorDescriptor& descriptor,
                         const Material& material,
                         const std::vector<double>& dimensions,
                         const VpdFitConfig& cfg);

// ---- Self-describing text I/O (CSV-with-#-comments, no JSON) ----

/// Write the fit to a self-describing text file. Returns false on open failure.
bool save_vpd(const VpdFit& fit, const std::string& filename);

/// Read a fit written by save_vpd. Throws std::runtime_error on parse failure.
VpdFit load_vpd(const std::string& filename);

} // namespace ceelo

#endif // CEELO_IO_VIRTUAL_DEPTH_FIT_H
