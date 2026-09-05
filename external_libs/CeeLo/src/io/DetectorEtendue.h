#ifndef CEELO_IO_DETECTOR_ETENDUE_H
#define CEELO_IO_DETECTOR_ETENDUE_H
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

/// @file DetectorEtendue.h
/// @brief Detector-side line set for extended-source efficiency integrals.
///
/// The per-element extended-source kernel (spec Eq. 5) integrates over SOURCE
/// points and, for each, over an aperture fan of rays toward the crystal:
///
///     eps = Int_V dV rho(r) P(r) (1/4pi) Int dOmega k(r, w) T(r, w)
///
/// Reversing the order of integration turns it into an integral over LINES that
/// hit the active crystal.  The phase-space measure seen from a source point,
/// dV dOmega, equals the etendue dA |w.n| dOmega ds seen from the crystal hull
/// (dA on the hull, n its outward normal, w the photon direction, s the
/// position along the line), so
///
///     eps = (1/4pi) Int_S dA Int dOmega |w.n| k(line) Int_chord ds rho P T
///
/// with k the same per-line detector kernel as ApertureQuadrature's
/// (sum over scoring segments of exp(-tau_before)(1 - exp(-mu l))).  The
/// detector side of every line - hull point, direction, material segments -
/// depends on neither the source nor the energy, so a fixed line set can be
/// built once and reused for every energy and every source configuration; the
/// source chord integral is analytic for a uniform source.  A host (InterSpec's
/// volumetric-source fit) does that reuse; this file only builds the lines and
/// evaluates the detector kernel along them.
///
/// Sampling (deterministic, Halton): a hull point on the crystal's convex hull
/// (front face + side wall; box: front + four sides), allocated among the faces
/// by projected area toward the target, then a direction in the cone that
/// bounds a target sphere as seen from that point (the source assembly's
/// bounding sphere, so a small far source is not starved of lines).  Directions
/// that leave the hull (w.n <= 0 for the OUTWARD direction) or point behind the
/// face plane get zero weight and are dropped WITHOUT renormalising - they are
/// counted in n, which is what keeps the estimator unbiased.  Lines with no
/// active-crystal chord (through the bore, a bulletized corner, a dead layer
/// only) are dropped too: they can never contribute.
///
/// Conventions match ResponseKernel.h: crystal front face at z = 0, crystal
/// along +z, source side at z < 0.  `KernelRay::dir` is the PHOTON direction
/// (into the crystal); `omega_w` is the line's etendue weight divided by 4 pi,
/// in cm^2 (not the dimensionless solid-angle fraction of an aperture
/// quadrature) - the host multiplies by a chord integral in cm to get an
/// efficiency.

#include "geometry/Geometry.h"
#include "io/ResponseKernel.h"

#include <Eigen/Core>

#include <cstdint>
#include <vector>

namespace ceelo {

/// A fixed set of lines through the active crystal; see the file docs.
struct EtendueLineSet {
    /// One KernelRay per kept line: `omega_w` = etendue weight / 4pi (cm^2),
    /// `dir` = photon direction (into the crystal), `segs` = the full trace
    /// from outside the outermost shell, sorted along the direction, `active_len`
    /// > 0.  `cone_omega_frac` is 0 (no single cone) and `n_rays_total` is the
    /// number of lines SAMPLED, dropped ones included.
    ApertureQuadrature q;

    /// Hull point and photon direction per kept line, double precision,
    /// parallel to `q.rays`.  A host intersects these with its own source
    /// geometry, where float directions would not do.
    std::vector<Eigen::Vector3d> origin;
    std::vector<Eigen::Vector3d> dir;

    /// The target sphere the directions were confined to (radius <= 0: the full
    /// sphere of directions was sampled at every hull point).
    Eigen::Vector3d target_centre{0.0, 0.0, 0.0};
    double target_radius = 0.0;

    /// Sum of the estimator's etendue weights over every sampled line that got
    /// a non-zero weight, BEFORE the active-chord drop (cm^2 sr).  In
    /// expectation this is Int_S dA Int_cone dOmega |w.n|; a measure check.
    double total_etendue = 0.0;
};

/// One sampled point on the crystal's convex hull.
struct HullPoint {
    Eigen::Vector3d point;    ///< cm, crystal frame
    Eigen::Vector3d normal;   ///< outward unit normal
    double area_weight = 0.0; ///< A_face / p_face (cm^2): this point's share of the hull area
};

/// `n` hull points (Halton bases 2, 3, 5 at index_offset + i).  The faces (front + side wall; box:
/// front + four sides) are allocated by projected area toward `toward` (a unit vector from the hull
/// centre toward the source region), floored at 5% of each face's true area so no face is starved;
/// `area_weight` carries the allocation, so the estimator is unbiased whatever it is.  With
/// `have_direction == false` the faces are allocated by plain area.
void sample_hull_points(const Geometry& geom, int n, const Eigen::Vector3d& toward,
                        bool have_direction, uint64_t index_offset,
                        std::vector<HullPoint>& out);

/// Trace the line through hull point `x` with OUTWARD unit direction `w_out` (toward the source
/// side; the photon travels along -w_out) and append it to `set` with the given etendue weight
/// (cm^2 sr: the caller's `area_weight * |w.n| / p(w) / n`, where p is the density its direction
/// proposal used).  Returns false - and appends nothing - when the line has no active-crystal
/// chord.  Requires w_out.n > 0 and w_out.z < 0 (asserted): a direction that does not leave the
/// hull toward the source side must be dropped by the caller WITHOUT renormalising.
bool append_etendue_line(EtendueLineSet& set, const Geometry& geom,
                         const Eigen::Vector3d& x, const Eigen::Vector3d& w_out,
                         double etendue_cm2_sr);

/// Build `n_lines` lines for `geom` aimed at the sphere (`target_centre_cm`,
/// `target_radius_cm`) in the crystal frame: hull points from sample_hull_points, directions
/// equal-solid-angle in the cone bounding the sphere as seen from each point (the full sphere of
/// directions when the point is inside it, or when `target_radius_cm <= 0`).  `index_offset`
/// shifts the Halton index range (a disjoint set for an independent estimate).  A host with a
/// better direction proposal (e.g. aimed at points of its own source solid) composes the two
/// functions above instead.
///
/// Cost: n_lines traces through the detector geometry (comparable to
/// make_aperture_quadrature at the same count).
EtendueLineSet build_etendue_lines(const Geometry& geom,
                                   const Eigen::Vector3d& target_centre_cm,
                                   double target_radius_cm,
                                   int n_lines,
                                   uint64_t index_offset = 0);

/// Per-line interaction probability with LIVE cross sections
/// (Material::macroscopic_xs): the chord-model kernel
///   p_i = sum_{scoring segs} exp(-tau_before(E)) (1 - exp(-mu*(E) l_seg))
/// exactly as ApertureQuadrature::interaction_omega accumulates it, but
/// returned per ray, PARALLEL to `q.rays` (0 for a ray without an active
/// chord).  So  sum_i q.rays[i].omega_w * p_out[i] == q.interaction_omega(...).
/// DetectorResponse::fep_line_probabilities is the stored-mu-table counterpart.
void line_interaction_probabilities(const ApertureQuadrature& q,
                                    double energy_keV, MuChoice mu,
                                    std::vector<double>& p_out,
                                    double passive_compton_recapture = 0.0);

}  // namespace ceelo

#endif  // CEELO_IO_DETECTOR_ETENDUE_H
