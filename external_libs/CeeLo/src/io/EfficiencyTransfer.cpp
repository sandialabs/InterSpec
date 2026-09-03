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

/// @file EfficiencyTransfer.cpp
/// @brief EFFTRAN-style efficiency transfer (see EfficiencyTransfer.h).

#include "io/EfficiencyTransfer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ceelo {

namespace {

// Segment index of E among ascending edges: number of edges strictly below E.
size_t segment_of(double e_keV, const std::vector<double>& edges) {
    size_t s = 0;
    for (const double x : edges)
        if (e_keV > x) ++s;
    return s;
}

}  // namespace

EfficiencyTransfer::EfficiencyTransfer(const Geometry& geom,
                                       const Eigen::Vector3d& ref_pos_cm,
                                       AnchorCurve anchor, MuChoice mu,
                                       int n_rays,
                                       std::vector<double> crystal_edges_keV,
                                       double passive_compton_recapture)
    : geom_(geom), ref_pos_(ref_pos_cm), mu_(mu), n_rays_(n_rays),
      recap_(passive_compton_recapture) {
    const size_t n = anchor.energies_keV.size();
    if (n == 0 || anchor.eff.size() != n)
        throw std::runtime_error("EfficiencyTransfer: empty/mismatched anchor");

    // Sort the anchor points by energy (the ln-eta interpolation needs an
    // ascending grid; callers usually pass one already).
    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return anchor.energies_keV[a] < anchor.energies_keV[b];
    });

    q_ref_ = make_aperture_quadrature(geom_, ref_pos_, n_rays_);

    // ln eta(E) = ln( eps_ref(E) / K(E, ref) ) at every anchor energy. The
    // reference kernel MUST fold the same recapture as eps_at's target kernel,
    // else the transfer ratio K(target)/K(ref) is inconsistent (and exactness
    // at the anchor distance is lost).
    std::vector<double> lnE, ln_eta;
    lnE.reserve(n);
    ln_eta.reserve(n);
    for (const size_t i : order) {
        const double E = anchor.energies_keV[i];
        const double eff = anchor.eff[i];
        const double K = q_ref_.interaction_omega(E, mu_, recap_);
        if (E <= 0.0 || eff <= 0.0 || K <= 0.0)
            throw std::runtime_error(
                "EfficiencyTransfer: non-positive anchor E/eff or zero kernel "
                "(is the reference position inside the detector, or the anchor "
                "efficiency zero?)");
        lnE.push_back(std::log(E));
        ln_eta.push_back(std::log(eff / K));
    }

    // Group into K-edge segments (both flanks become segment endpoints); build
    // a monotone-cubic curve per segment so no interpolant bridges an edge.
    std::sort(crystal_edges_keV.begin(), crystal_edges_keV.end());
    size_t i0 = 0;
    while (i0 < n) {
        const size_t seg = segment_of(std::exp(lnE[i0]), crystal_edges_keV);
        size_t i1 = i0;
        while (i1 + 1 < n &&
               segment_of(std::exp(lnE[i1 + 1]), crystal_edges_keV) == seg)
            ++i1;
        EtaSeg es;
        es.lo = lnE[i0];
        es.hi = lnE[i1];
        if (i1 > i0) {
            es.curve = Pchip(std::vector<double>(lnE.begin() + i0,
                                                 lnE.begin() + i1 + 1),
                             std::vector<double>(ln_eta.begin() + i0,
                                                 ln_eta.begin() + i1 + 1));
        } else {
            es.constant = ln_eta[i0];  // single node in this segment
        }
        eta_segs_.push_back(std::move(es));
        i0 = i1 + 1;
    }
}

double EfficiencyTransfer::eta_ref(double energy_keV) const {
    if (eta_segs_.empty()) return 0.0;
    const double lnE = std::log(energy_keV);
    // Pick the segment whose lnE range contains the query; else the nearest.
    const EtaSeg* best = &eta_segs_.front();
    double best_dist = std::numeric_limits<double>::max();
    for (const EtaSeg& s : eta_segs_) {
        if (lnE >= s.lo - 1e-12 && lnE <= s.hi + 1e-12) { best = &s; best_dist = 0.0; break; }
        const double d = std::min(std::fabs(lnE - s.lo), std::fabs(lnE - s.hi));
        if (d < best_dist) { best_dist = d; best = &s; }
    }
    const double q = std::min(std::max(lnE, best->lo), best->hi);  // clamp
    const double ln_eta = best->curve.valid() ? best->curve(q) : best->constant;
    return std::exp(ln_eta);
}

double EfficiencyTransfer::k_ref(double energy_keV) const {
    return q_ref_.interaction_omega(energy_keV, mu_, recap_);
}

double EfficiencyTransfer::eps_at(double energy_keV,
                                  const ApertureQuadrature& q_target) const {
    return eta_ref(energy_keV) *
           q_target.interaction_omega(energy_keV, mu_, recap_);
}

double EfficiencyTransfer::eps_at(double energy_keV,
                                  const Eigen::Vector3d& target_pos_cm) const {
    // Behind the crystal-face plane the removal-only K ratio is not valid
    // (side/back entry); refuse rather than return a biased guess.
    const double d = target_pos_cm.norm();
    const double cos_theta = (d > 0.0) ? (-target_pos_cm.z() / d) : 1.0;
    if (cos_theta < 0.0) return std::numeric_limits<double>::quiet_NaN();
    return eps_at(energy_keV, make_target_quadrature(target_pos_cm));
}

// ---------------------------------------------------------------------------
// make_transfer_response
// ---------------------------------------------------------------------------

std::shared_ptr<DetectorResponse> make_transfer_response(
    const GeometryDescriptor& descriptor, const AnchorCurve& fep_anchor,
    const Eigen::Vector3d& ref_pos_cm, const AnchorCurve* tot_anchor,
    const TransferResponseOptions& opts) {
    if (fep_anchor.energies_keV.size() < 2 ||
        fep_anchor.eff.size() != fep_anchor.energies_keV.size())
        throw std::runtime_error(
            "make_transfer_response: need >= 2 ascending anchor points");

    auto resp = std::make_shared<DetectorResponse>();
    resp->descriptor = descriptor;

    // Stored mu tables (self-contained; eta*K cannot drift with the XS data).
    for (size_t i = 0; i < descriptor.materials.size(); ++i) {
        const Material m = descriptor.materials[i].to_material();
        resp->mu_tables.push_back(MuTable::sample(m, static_cast<int>(i)));
    }

    const double e_min = opts.e_min_keV > 0.0 ? opts.e_min_keV
                                              : fep_anchor.energies_keV.front();
    const double e_max = opts.e_max_keV > 0.0 ? opts.e_max_keV
                                              : fep_anchor.energies_keV.back();
    resp->provenance.method = ProductionMethod::CurveTransfer;  //no Monte Carlo at all
    resp->provenance.profile = ResponseProfile::FarField;
    resp->provenance.kernel_n_rays = opts.kernel_n_rays;
    resp->provenance.detector_name = opts.detector_name;
    resp->provenance.valid_e_min_keV = e_min;
    resp->provenance.valid_e_max_keV = e_max;

    // Build geometry + mu lookups ONCE (finalize appends owned materials).
    resp->finalize();

    const double a = descriptor.transverse_half_extent();
    resp->provenance.min_distance_cm = 2.0 * a;  // far-field validity floor

    // Dead-layer/endcap scatter-in credit for the total kernel. Set BEFORE any
    // kernel_K(NoRayleigh) call so the anchor b(E) and every query fold it with
    // the same coefficient (the ratio is what corrects the near-field total).
    resp->scatter_in_recapture = opts.scatter_in_recapture;

    const ApertureQuadrature q_ref = resp->make_quadrature(ref_pos_cm);
    const std::vector<double> edges = descriptor.crystal_k_edges(e_min, e_max);

    // Angle-flat eta table: two identical sentinel cos-theta nodes.
    EtaTable& t = resp->eta_fep;
    t.energies_keV = fep_anchor.energies_keV;
    t.cos_thetas = {opts.cos_theta_lo, 1.0};
    t.edges_keV = edges;
    const size_t ne = t.energies_keV.size();
    t.ln_eta.assign(ne * 2, 0.0);
    t.frac_sigma.assign(ne * 2, 0.0);
    for (size_t e = 0; e < ne; ++e) {
        const double E = t.energies_keV[e];
        const double K = resp->kernel_K(E, q_ref, MuChoice::Total);
        if (E <= 0.0 || fep_anchor.eff[e] <= 0.0 || K <= 0.0)
            throw std::runtime_error("make_transfer_response: zero kernel/eff");
        const double ln_eta = std::log(fep_anchor.eff[e] / K);
        const double fs = e < fep_anchor.frac_sigma.size()
                              ? fep_anchor.frac_sigma[e] : 0.0;
        for (size_t c = 0; c < 2; ++c) {
            t.ln_eta[t.index(e, c, 0)] = ln_eta;
            t.frac_sigma[t.index(e, c, 0)] = fs;
        }
    }
    t.finalize();

    // Total-efficiency tier: b(E) from the total anchor (angle-flat), else the
    // bare-crystal kernel-exact tier.
    if (tot_anchor && tot_anchor->energies_keV.size() >= 2 &&
        tot_anchor->eff.size() == tot_anchor->energies_keV.size()) {
        resp->tot_eff.tier = TotEffTier::BCurve;
        for (size_t i = 0; i < tot_anchor->energies_keV.size(); ++i) {
            const double E = tot_anchor->energies_keV[i];
            const double K_nors = resp->kernel_K(E, q_ref, MuChoice::NoRayleigh);
            if (E <= 0.0 || tot_anchor->eff[i] <= 0.0 || K_nors <= 0.0)
                continue;
            resp->tot_eff.b_energies_keV.push_back(E);
            resp->tot_eff.ln_b.push_back(std::log(tot_anchor->eff[i] / K_nors));
        }
        if (resp->tot_eff.b_energies_keV.size() < 2)
            resp->tot_eff.tier = TotEffTier::KernelExact;  // fall back
    } else {
        resp->tot_eff.tier = TotEffTier::KernelExact;
    }
    resp->tot_eff.finalize();

    // Honest off-axis/near sigma (the angle-flat eta has no theta residual).
    resp->model_transfer = opts.model_transfer;

    return resp;
}

} // namespace ceelo
