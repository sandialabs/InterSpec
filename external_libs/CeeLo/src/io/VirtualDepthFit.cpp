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

#include "io/VirtualDepthFit.h"

#include "io/SolidAngle.h"
#include "efficiency/EfficiencyCalculator.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ceelo {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// u = ln(E in MeV); E supplied in keV (InterSpec coefficient convention).
inline double log_energy_mev(double energy_keV) {
    return std::log(energy_keV / 1000.0);
}

} // namespace

DetectorDescriptor make_descriptor(std::string name, DetectorShape shape,
                                   const Material& material,
                                   const std::vector<double>& dimensions) {
    DetectorDescriptor d;
    d.name = std::move(name);
    d.material_name = material.name();
    d.shape = shape;
    if (shape == DetectorShape::Cylinder) {
        // dimensions = {radius, half_length}
        d.crystal_radius_cm = dimensions.at(0);
        d.crystal_length_cm = 2.0 * dimensions.at(1);
    } else {
        // dimensions = {half_x, half_y, length}
        d.half_x_cm = dimensions.at(0);
        d.half_y_cm = dimensions.at(1);
        d.crystal_length_cm = dimensions.at(2);
        // Area-equivalent disk radius for the rectangular face (area 2hx * 2hy):
        //   pi R_eff^2 = 4 hx hy  =>  R_eff = 2 sqrt(hx hy / pi)
        d.crystal_radius_cm = 2.0 * std::sqrt(d.half_x_cm * d.half_y_cm / kPi);
    }
    d.crystal_diameter_cm = 2.0 * d.crystal_radius_cm;
    return d;
}

VpdPoint fit_delta_one_energy(double R_cm,
                              const std::vector<double>& dist_cm,
                              const std::vector<double>& eps,
                              const std::vector<double>& eps_unc,
                              double delta_max_cm,
                              double delta_step_cm) {
    VpdPoint pt;
    pt.dist_cm = dist_cm;
    pt.eps = eps;
    pt.eps_unc = eps_unc;

    double best_delta = 0.0, best_cv = 1e30, best_mean = 0.0;
    for (double delta = 0.0; delta <= delta_max_cm + 1e-9; delta += delta_step_cm) {
        double s = 0.0, s2 = 0.0;
        int n = 0;
        for (size_t i = 0; i < dist_cm.size(); ++i) {
            const double f = disk_solid_angle_fraction(dist_cm[i] + delta, R_cm);
            if (f <= 0.0) continue;
            const double v = eps[i] / f;
            s += v;
            s2 += v * v;
            ++n;
        }
        if (n == 0) continue;
        const double mean = s / n;
        const double var = s2 / n - mean * mean;
        const double cv = (mean > 0.0) ? std::sqrt(std::max(0.0, var)) / mean : 1e30;
        if (cv < best_cv) {
            best_cv = cv;
            best_delta = delta;
            best_mean = mean;
        }
    }

    pt.delta_cm = best_delta;
    pt.intrinsic_mean = best_mean;
    pt.cv_residual = best_cv;

    double maxres = 0.0;
    for (size_t i = 0; i < dist_cm.size(); ++i) {
        if (eps[i] <= 0.0) continue;
        const double pred = best_mean * disk_solid_angle_fraction(dist_cm[i] + best_delta, R_cm);
        maxres = std::max(maxres, std::abs(pred - eps[i]) / eps[i]);
    }
    pt.max_abs_res = maxres;
    return pt;
}

std::vector<double> fit_exp_log_coeffs(const std::vector<double>& energies_keV,
                                       const std::vector<double>& delta_cm,
                                       int order,
                                       double delta_floor_cm) {
    const size_t n = energies_keV.size();
    if (n == 0 || order < 1) return {};
    // Cannot fit more coefficients than data points.
    const int eff_order = std::min<int>(order, static_cast<int>(n));

    Eigen::MatrixXd A(n, eff_order);
    Eigen::VectorXd b(n);
    for (size_t k = 0; k < n; ++k) {
        const double u = log_energy_mev(energies_keV[k]);
        double pow_u = 1.0;
        for (int j = 0; j < eff_order; ++j) {
            A(k, j) = pow_u;
            pow_u *= u;
        }
        b(k) = std::log(std::max(delta_cm[k], delta_floor_cm));
    }

    const Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);
    std::vector<double> coeffs(static_cast<size_t>(eff_order));
    for (int j = 0; j < eff_order; ++j) coeffs[j] = c(j);
    return coeffs;
}

double eval_exp_log(const std::vector<double>& coeffs, double energy_keV) {
    if (coeffs.empty()) return 0.0;
    const double u = log_energy_mev(energy_keV);
    double sum = 0.0, pow_u = 1.0;
    for (double c : coeffs) {
        sum += c * pow_u;
        pow_u *= u;
    }
    return std::exp(sum);
}

double VpdFit::delta_at(double energy_keV) const {
    double e = energy_keV;
    if (valid_e_min_keV > 0.0 && e < valid_e_min_keV) e = valid_e_min_keV;
    if (valid_e_max_keV > 0.0 && e > valid_e_max_keV) e = valid_e_max_keV;
    return eval_exp_log(coeffs, e);
}

VpdFit fit_virtual_depth(const DetectorDescriptor& descriptor,
                         const Material& material,
                         const std::vector<double>& dimensions,
                         const VpdFitConfig& cfg) {
    VpdFit fit;
    fit.detector = descriptor;

    EfficiencyCalculator calc;
    calc.set_detector(descriptor.shape, &material, dimensions);

    // Sort energies ascending for a well-conditioned validity range / readout.
    std::vector<double> energies = cfg.energies_keV;
    std::sort(energies.begin(), energies.end());

    std::vector<double> fit_e, fit_delta;
    for (double E : energies) {
        std::vector<double> eps, unc;
        eps.reserve(cfg.distance_ladder_cm.size());
        unc.reserve(cfg.distance_ladder_cm.size());
        for (double d : cfg.distance_ladder_cm) {
            calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -d));
            const EfficiencyResult r =
                calc.compute(E, cfg.events_per_point, cfg.num_threads);
            eps.push_back(r.full_energy_peak_efficiency);
            unc.push_back(r.fep_uncertainty);
        }
        VpdPoint pt = fit_delta_one_energy(descriptor.crystal_radius_cm,
                                           cfg.distance_ladder_cm, eps, unc,
                                           cfg.delta_max_cm, cfg.delta_step_cm);
        pt.energy_keV = E;
        fit.points.push_back(pt);
        fit_e.push_back(E);
        fit_delta.push_back(pt.delta_cm);
    }

    fit.coeffs = fit_exp_log_coeffs(fit_e, fit_delta, cfg.poly_order);
    fit.valid_e_min_keV = energies.empty() ? 0.0 : energies.front();
    fit.valid_e_max_keV = energies.empty() ? 0.0 : energies.back();

    // Residuals of the exp-of-log fit (in ln delta).
    double ss = 0.0;
    for (size_t k = 0; k < fit_e.size(); ++k) {
        const double model = eval_exp_log(fit.coeffs, fit_e[k]);
        const double res = std::log(std::max(fit_delta[k], 1e-3)) -
                           std::log(std::max(model, 1e-3));
        fit.log_fit_residuals.push_back(res);
        ss += res * res;
    }
    fit.rms_log_residual = fit_e.empty() ? 0.0 : std::sqrt(ss / fit_e.size());
    return fit;
}

// ---------------------------------------------------------------------------
// Text I/O
// ---------------------------------------------------------------------------

namespace {

std::string shape_to_string(DetectorShape s) {
    return s == DetectorShape::Box ? "Box" : "Cylinder";
}

DetectorShape shape_from_string(const std::string& s) {
    return s == "Box" ? DetectorShape::Box : DetectorShape::Cylinder;
}

/// If `line` is "# key: value", return true and fill key/value (trimmed).
bool parse_meta(const std::string& line, std::string& key, std::string& value) {
    if (line.empty() || line[0] != '#') return false;
    const size_t colon = line.find(':');
    if (colon == std::string::npos) return false;
    size_t ks = 1;
    while (ks < colon && std::isspace(static_cast<unsigned char>(line[ks]))) ++ks;
    key = line.substr(ks, colon - ks);
    size_t vs = colon + 1;
    while (vs < line.size() && std::isspace(static_cast<unsigned char>(line[vs]))) ++vs;
    value = line.substr(vs);
    // trim trailing whitespace
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back())))
        value.pop_back();
    return true;
}

std::vector<double> parse_csv_doubles(const std::string& line) {
    std::vector<double> out;
    std::stringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        try {
            out.push_back(std::stod(tok));
        } catch (...) {
            return {}; // not a numeric row
        }
    }
    return out;
}

} // namespace

bool save_vpd(const VpdFit& fit, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) return false;
    out << std::setprecision(8);

    const DetectorDescriptor& d = fit.detector;
    out << "# CeeLo virtual interaction depth delta(E) export\n";
    out << "# format: 1\n";
    out << "# detector: " << d.name << "\n";
    out << "# material: " << d.material_name << "\n";
    out << "# shape: " << shape_to_string(d.shape) << "\n";
    out << "# crystal_radius_cm: " << d.crystal_radius_cm << "\n";
    out << "# crystal_diameter_cm: " << d.crystal_diameter_cm << "\n";
    out << "# crystal_length_cm: " << d.crystal_length_cm << "\n";
    out << "# half_x_cm: " << d.half_x_cm << "\n";
    out << "# half_y_cm: " << d.half_y_cm << "\n";
    out << "# valid_energy_keV: " << fit.valid_e_min_keV << " " << fit.valid_e_max_keV << "\n";
    out << "# model: delta_cm = exp(p0 + p1*u + p2*u^2 + ...), u = ln(E_MeV)\n";
    out << "# poly_order: " << fit.coeffs.size() << "\n";
    out << "# rms_log_residual: " << fit.rms_log_residual << "\n";
    out << "# solid_angle: disk, Omega/4pi = 0.5*(1 - x/sqrt(x^2+R^2)), x = d + delta\n";

    out << "COEFFS\n";
    for (size_t j = 0; j < fit.coeffs.size(); ++j)
        out << (j ? "," : "") << "p" << j;
    out << "\n";
    for (size_t j = 0; j < fit.coeffs.size(); ++j)
        out << (j ? "," : "") << fit.coeffs[j];
    out << "\n";

    out << "POINTS\n";
    out << "# E_keV,delta_cm,intrinsic_mean,cv_residual,max_abs_res_pct,log_fit_residual\n";
    for (size_t k = 0; k < fit.points.size(); ++k) {
        const VpdPoint& p = fit.points[k];
        const double res = (k < fit.log_fit_residuals.size()) ? fit.log_fit_residuals[k] : 0.0;
        out << p.energy_keV << "," << p.delta_cm << "," << p.intrinsic_mean << ","
            << p.cv_residual << "," << (100.0 * p.max_abs_res) << "," << res << "\n";
    }

    out << "RAWGRID\n";
    out << "# E_keV,dist_cm,eps,eps_unc\n";
    for (const VpdPoint& p : fit.points) {
        for (size_t i = 0; i < p.dist_cm.size(); ++i) {
            out << p.energy_keV << "," << p.dist_cm[i] << "," << p.eps[i] << ","
                << (i < p.eps_unc.size() ? p.eps_unc[i] : 0.0) << "\n";
        }
    }
    return static_cast<bool>(out);
}

VpdFit load_vpd(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("load_vpd: cannot open " + filename);

    VpdFit fit;
    DetectorDescriptor& d = fit.detector;

    enum Section { NONE, COEFFS, POINTS, RAWGRID } section = NONE;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::string key, value;
        if (line[0] == '#') {
            if (!parse_meta(line, key, value)) continue;
            if (key == "detector") d.name = value;
            else if (key == "material") d.material_name = value;
            else if (key == "shape") d.shape = shape_from_string(value);
            else if (key == "crystal_radius_cm") d.crystal_radius_cm = std::stod(value);
            else if (key == "crystal_diameter_cm") d.crystal_diameter_cm = std::stod(value);
            else if (key == "crystal_length_cm") d.crystal_length_cm = std::stod(value);
            else if (key == "half_x_cm") d.half_x_cm = std::stod(value);
            else if (key == "half_y_cm") d.half_y_cm = std::stod(value);
            else if (key == "rms_log_residual") fit.rms_log_residual = std::stod(value);
            else if (key == "valid_energy_keV") {
                std::stringstream ss(value);
                ss >> fit.valid_e_min_keV >> fit.valid_e_max_keV;
            }
            continue;
        }

        if (line == "COEFFS") { section = COEFFS; continue; }
        if (line == "POINTS") { section = POINTS; continue; }
        if (line == "RAWGRID") { section = RAWGRID; continue; }

        if (section == COEFFS) {
            // Skip the "p0,p1,.." header; take the first numeric row.
            std::vector<double> v = parse_csv_doubles(line);
            if (!v.empty()) fit.coeffs = v;
        } else if (section == POINTS) {
            std::vector<double> v = parse_csv_doubles(line);
            if (v.size() >= 6) {
                VpdPoint p;
                p.energy_keV = v[0];
                p.delta_cm = v[1];
                p.intrinsic_mean = v[2];
                p.cv_residual = v[3];
                p.max_abs_res = v[4] / 100.0;
                fit.points.push_back(p);
                fit.log_fit_residuals.push_back(v[5]);
            }
        } else if (section == RAWGRID) {
            std::vector<double> v = parse_csv_doubles(line);
            if (v.size() >= 4) {
                // Append to the matching point (by energy).
                for (VpdPoint& p : fit.points) {
                    if (std::abs(p.energy_keV - v[0]) < 1e-6) {
                        p.dist_cm.push_back(v[1]);
                        p.eps.push_back(v[2]);
                        p.eps_unc.push_back(v[3]);
                        break;
                    }
                }
            }
        }
    }
    return fit;
}

} // namespace ceelo
