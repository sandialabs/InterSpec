#pragma once
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

/// @file ConvexSeparation.h
/// @brief How far apart two convex solids are, and where a ray enters one - the test of whether a
///        source (with its shielding) stays clear of a detector's outside
///        (GeometryDescriptor::Envelope::solid()), wherever it sits around the detector.

#include <Eigen/Core>

#include <cstdint>

namespace ceelo {

/// A point, sphere, finite cylinder or box, in any position and orientation.  Units are whatever
/// the caller uses consistently (cm throughout CeeLo).
struct ConvexSolid {
    enum class Kind : uint8_t { Point, Sphere, Cylinder, Box };

    Kind kind = Kind::Point;
    Eigen::Vector3d center = Eigen::Vector3d::Zero();

    double radius = 0.0;                                     ///< sphere and cylinder
    Eigen::Vector3d axis = Eigen::Vector3d::UnitZ();         ///< cylinder: unit axis
    double half_length = 0.0;                                ///< cylinder: half-length along `axis`
    Eigen::Vector3d half_extents = Eigen::Vector3d::Zero();  ///< box: along `rotation`'s columns
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();  ///< box: local -> world

    static ConvexSolid point(const Eigen::Vector3d& c);
    static ConvexSolid sphere(const Eigen::Vector3d& c, double r);
    static ConvexSolid cylinder(const Eigen::Vector3d& c, const Eigen::Vector3d& axis, double r,
                                double half_length);
    static ConvexSolid box(const Eigen::Vector3d& c, const Eigen::Vector3d& half_extents,
                           const Eigen::Matrix3d& rotation = Eigen::Matrix3d::Identity());

    /// The point of this solid nearest `p` (`p` itself when inside).
    Eigen::Vector3d project(const Eigen::Vector3d& p) const;

    /// This solid moved rigidly: a point x of it goes to `R * x + t`.  `R` must be a rotation.
    ConvexSolid transformed(const Eigen::Matrix3d& R, const Eigen::Vector3d& t) const;
};

/// Distance between two convex solids - 0 when they overlap - and optionally the closest point of
/// each.  Alternating projections: for disjoint compact convex sets they converge to a closest
/// pair, and for overlapping ones to a common point, stopping once the pair moves by less than `tol`
/// (or after `max_iter` passes).  Near tangency - a curved face against a flat one - convergence is
/// slow: solids that touch, or overlap by up to about 1e-4 of their size, can read a small positive
/// separation (micrometres for detector-sized solids).  So a separation is a clearance, not an exact
/// touching test; it is run from each solid and the closer answer kept.
double convex_separation(const ConvexSolid& a, const ConvexSolid& b,
                         Eigen::Vector3d* closest_a = nullptr, Eigen::Vector3d* closest_b = nullptr,
                         double tol = 1.0e-10, int max_iter = 5000);

/// Distance from `origin`, along the unit direction `dir`, to where the ray enters `solid` - 0 when
/// `origin` is inside it, +infinity when the ray misses it.  A point solid is never entered.
double convex_entry_distance(const ConvexSolid& solid, const Eigen::Vector3d& origin,
                             const Eigen::Vector3d& dir);

}  // namespace ceelo
