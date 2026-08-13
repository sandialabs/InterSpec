#ifndef CEELO_IO_EFFICIENCY_TRANSFER_H
#define CEELO_IO_EFFICIENCY_TRANSFER_H
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

/// @file EfficiencyTransfer.h
/// @brief EFFTRAN-style efficiency transfer: anchor an efficiency curve at ONE
/// reference source position, then transfer it to any other position by the
/// ratio of the interaction-weighted effective solid angle.
///
/// EFFTRAN transfers a known efficiency curve between geometries via
///     eps_target(E) = eps_ref(E) * T_target(E) / T_ref(E)
/// where T(E) is the attenuation-weighted effective solid angle with chord
/// saturation. In CeeLo terms T(E) is exactly
/// `ApertureQuadrature::interaction_omega(E, MuChoice)` (io/ResponseKernel.h) --
/// the interaction-weighted effective solid angle over the EXACT ray-traced
/// detector geometry (`MuChoice::Total` = FEP kernel; `MuChoice::NoRayleigh` =
/// total-efficiency kernel). So the transfer is precisely holding the response
/// factorization's residual eta(E) = eps_ref(E)/K(E, ref) constant across
/// positions (DetectorResponse.h Eq. 3, with eta angle-flat). Unlike EFFTRAN
/// this kernel is fully 3-D (off-axis native), Rayleigh-aware, and chord-exact.
///
/// This header provides two things:
///   1. `EfficiencyTransfer` -- the standalone, MC-free transfer primitive on a
///      live `Geometry` (live cross sections). Anchor from CeeLo MC results OR
///      user-measured (E, eff) points (the classic EFFTRAN lab workflow).
///   2. `make_transfer_response` -- the same idea packaged as a normal
///      `DetectorResponse` (stored mu tables, angle-flat eta), so XML,
///      extended-source, cascade, grounding and the sigma machinery all work.
///
/// Validity: the K ratio transfers geometry to ~1% down to d ~ 1.1-1.8a
/// (campaign D1) FOR A FIXED ANGLE (distance transfer). Off-axis the residual
/// eta(E, theta) varies (campaign D2: ~5-12 cos-theta nodes for tight 1%), so
/// pure angle-flat transfer degrades off-axis -- `make_transfer_response`
/// therefore attaches a SigmaTransferModel that inflates sigma off-axis/near.

#include "geometry/Geometry.h"
#include "io/DetectorResponse.h"
#include "io/Pchip.h"
#include "io/ResponseKernel.h"

#include <Eigen/Core>

#include <functional>
#include <memory>
#include <vector>

namespace ceelo {

/// Anchor efficiency curve: (E, eff[, frac_sigma]) at ONE reference position.
/// Supplied from CeeLo MC results OR user-measured points (true EFFTRAN).
struct AnchorCurve {
    std::vector<double> energies_keV;   ///< strictly ascending
    std::vector<double> eff;            ///< absolute efficiency at the ref pos
    std::vector<double> frac_sigma;     ///< optional per-point (empty => 0)
};

/// Standalone EFFTRAN transfer primitive over a live detector `Geometry`.
/// Cheap: one aperture quadrature is built at construction; each transfer is a
/// pair of `interaction_omega` evaluations (sub-ms). Not thread-safe to build,
/// but the const eval methods are safe once constructed.
class EfficiencyTransfer {
public:
    /// @param geom              live detector geometry (crystal + layers).
    /// @param ref_pos_cm        reference source point, crystal-face frame.
    /// @param anchor            efficiency curve known at ref_pos_cm.
    /// @param mu                Total = FEP transfer; NoRayleigh = total-eff.
    /// @param n_rays            aperture-quadrature ray count.
    /// @param crystal_edges_keV crystal K-edges to segment the ln-eta E-interp
    ///                          (no interpolant bridges an edge); empty = one
    ///                          global segment.
    /// @param passive_compton_recapture dead-layer/endcap scatter-in credit for
    ///                          the total kernel (mu == NoRayleigh); ignored for
    ///                          FEP. Pass kTotalScatterInRecapture for total
    ///                          transfer on thick-dead-layer detectors.
    EfficiencyTransfer(const Geometry& geom, const Eigen::Vector3d& ref_pos_cm,
                       AnchorCurve anchor, MuChoice mu = MuChoice::Total,
                       int n_rays = 8192,
                       std::vector<double> crystal_edges_keV = {},
                       double passive_compton_recapture = 0.0);

    /// eps(E, target) = eta_ref(E) * K(E, target). The position overload
    /// refuses behind-plane targets (cos_theta < 0): returns NaN.
    double eps_at(double energy_keV, const Eigen::Vector3d& target_pos_cm) const;
    double eps_at(double energy_keV, const ApertureQuadrature& q_target) const;

    /// MEFFTRAN: extended/shielded target via a per-ray source-transmission
    /// factor folded into the target kernel (build t_src from
    /// fep_survival_removal_mu for FEP -- with the material-aware
    /// kn_in_window_fraction(E, win, mat) f_win when the source/shield
    /// material is known -- mu_total for total; see
    /// DetectorResponse::eps_fep_element).
    template <typename TsrcFn>
    double eps_at_with_tsrc(double energy_keV, const ApertureQuadrature& q_target,
                            TsrcFn&& t_src) const {
        return eta_ref(energy_keV) *
               q_target.interaction_omega_with_tsrc(
                   energy_keV, mu_, std::forward<TsrcFn>(t_src), recap_);
    }

    /// The transferable residual eta(E) = eps_ref(E)/K(E, ref) (E-interpolated,
    /// K-edge-segmented), and the reference kernel K(E, ref).
    double eta_ref(double energy_keV) const;
    double k_ref(double energy_keV) const;

    ApertureQuadrature make_target_quadrature(const Eigen::Vector3d& t) const {
        return make_aperture_quadrature(geom_, t, n_rays_);
    }

    /// Collimator validity gate at a target: s = transmitted/geometric active
    /// solid angle. s < 0.3 => scatter/leak the removal-only ratio can't
    /// represent (warn); s < 0.05 => refuse. See ResponseKernel hole_fraction.
    double hole_fraction(double energy_keV, const ApertureQuadrature& q) const {
        return q.hole_fraction(energy_keV);
    }

    MuChoice mu() const { return mu_; }

private:
    // One K-edge segment of the ln-eta(lnE) curve.
    struct EtaSeg {
        double lo = 0.0, hi = 0.0;   // lnE range (tiles the anchor span)
        Pchip curve;                 // valid() when >= 2 nodes
        double constant = 0.0;       // used when the segment has 1 node
    };

    const Geometry& geom_;
    Eigen::Vector3d ref_pos_;
    MuChoice mu_;
    int n_rays_;
    double recap_ = 0.0;             // passive Compton scatter-in credit (total)
    ApertureQuadrature q_ref_;       // built once at ref_pos_
    std::vector<EtaSeg> eta_segs_;   // segmented ln-eta(lnE)
};

/// Default dead-layer/endcap Compton scatter-in recapture fraction for TOTAL
/// efficiency transfer. Calibrated against HPGe near-field total MC
/// (efficiency-transfer study): the coax-HPGe total near-field under-prediction drops
/// from ~2.6% at contact (r=0) to ~1.2% at r=0.8, monotonically and without
/// over-correcting (r=1.0 still leaves ~0.9%), while thin-can NaI (small passive
/// path) is unaffected. Physically the near-crystal passive layers recapture
/// most Compton interactions (forward photon + recoil electron into the
/// contiguous active volume), so r is high but held sub-unity for safety on
/// thicker passive stacks. FEP is unaffected (peak photons cannot be scatter-in).
inline constexpr double kTotalScatterInRecapture = 0.8;

// ---------------------------------------------------------------------------
// DetectorResponse producer (the recommended stored artifact)
// ---------------------------------------------------------------------------

struct TransferResponseOptions {
    double cos_theta_lo = 0.02;      ///< low sentinel angular node
    double e_min_keV = 0.0;          ///< 0 => span the anchor
    double e_max_keV = 0.0;
    int kernel_n_rays = 2048;        ///< evaluation-time quadrature rays
    SigmaTransferModel model_transfer{};  ///< off-axis/near sigma envelope
    /// Total-kernel scatter-in recapture (see kTotalScatterInRecapture); only
    /// affects eps_total, and only for detectors with passive layers.
    double scatter_in_recapture = kTotalScatterInRecapture;
    std::string detector_name;
};

/// Build a self-contained `DetectorResponse` from a single anchor curve: an
/// angle-flat eta table (two identical sentinel cos-theta nodes) plus the
/// query-time kernel over the stored mu tables. `tot_anchor` (optional) adds an
/// angle-flat total-efficiency b(E) tier; without it the total efficiency is
/// the bare-crystal kernel-exact tier. The result flags off-axis/near queries
/// via the attached SigmaTransferModel and provenance min-distance.
std::shared_ptr<DetectorResponse> make_transfer_response(
    const GeometryDescriptor& descriptor, const AnchorCurve& fep_anchor,
    const Eigen::Vector3d& ref_pos_cm, const AnchorCurve* tot_anchor = nullptr,
    const TransferResponseOptions& opts = {});

} // namespace ceelo

#endif // CEELO_IO_EFFICIENCY_TRANSFER_H
