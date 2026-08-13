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

#include "geometry/SourceGeometry.h"
#include "geometry/Cylinder.h"
#include "geometry/Box.h"
#include "geometry/Sphere.h"
#include "physics/ElectronCsda.h"
#include "transport/ComptonScatter.h"
#include "transport/TransportUtils.h"

#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>

namespace ceelo {

namespace {

/// Path length of a ray through a spherical shell [r_in, r_out] centered at
/// `center`. Sums the near-wall and far-wall crossings (a ray born in the shell
/// can cross the inner void and re-enter the far wall). r_in = 0 → solid ball.
/// Mirrors the inner+outer sphere logic the point-source shield loop has always
/// used, factored so the spherical *source material* and shields share it.
double spherical_shell_path(const Eigen::Vector3d& center,
                            const Eigen::Vector3d& position,
                            const Eigen::Vector3d& direction,
                            double r_in, double r_out) {
    auto outer_hit = intersect_sphere(position, direction, center, r_out);
    if (!outer_hit || !outer_hit->valid()) return 0.0;

    double t0_out = std::max(outer_hit->t_enter, 0.0);
    double t1_out = outer_hit->t_exit;

    if (r_in > 1e-10) {
        auto inner_hit = intersect_sphere(position, direction, center, r_in);
        if (inner_hit && inner_hit->valid()) {
            double path = 0.0;
            double t1_in = inner_hit->t_exit;
            if (t1_in < t1_out) path = t1_out - t1_in;
            double t0_in = std::max(inner_hit->t_enter, 0.0);
            if (t0_out < t0_in) path += t0_in - t0_out;
            return path;
        }
    }
    return t1_out - t0_out;
}

/// Path length of a ray (in the cylinder-local frame) through the region between
/// an inner and outer coaxial cylinder. Inner cylinder = (r_in, half-length
/// hz_in); outer = (r_out, hz_out). r_in = 0 → solid. Mirrors the inner+outer
/// tube logic the cylindrical shield loop uses; reused for the annular
/// (hollow-bore) source material where hz_in == hz_out.
double cyl_shell_path(const Eigen::Vector3d& local_pos,
                      const Eigen::Vector3d& local_dir,
                      double r_in, double hz_in,
                      double r_out, double hz_out) {
    auto outer_hit = intersect_cylinder(local_pos, local_dir, r_out, -hz_out, hz_out);
    if (!outer_hit || !outer_hit->valid()) return 0.0;

    double t0_out = std::max(outer_hit->t_enter, 0.0);
    double t1_out = outer_hit->t_exit;

    if (r_in > 1e-10) {
        auto inner_hit = intersect_cylinder(local_pos, local_dir, r_in, -hz_in, hz_in);
        if (inner_hit && inner_hit->valid()) {
            double path = 0.0;
            double t1_in = inner_hit->t_exit;
            if (t1_in < t1_out) path = t1_out - t1_in;
            double t0_in = std::max(inner_hit->t_enter, 0.0);
            if (t0_out < t0_in) path += t0_in - t0_out;
            return path;
        }
    }
    return t1_out - t0_out;
}

/// Path length of a ray (in the box-local frame) through the shell between an
/// inner and outer centered box. Inner all-zero → solid. Mirrors cyl_shell_path
/// (and the inner+outer box logic the rectangular shield loop uses); used for
/// the hollow-box (crate/container-wall) source material.
double box_shell_path(const Eigen::Vector3d& local_pos,
                      const Eigen::Vector3d& local_dir,
                      const Eigen::Vector3d& inner_hd,
                      const Eigen::Vector3d& outer_hd) {
    auto outer_hit = intersect_box(local_pos, local_dir,
                                   outer_hd.x(), outer_hd.y(),
                                   -outer_hd.z(), outer_hd.z());
    if (!outer_hit || !outer_hit->valid()) return 0.0;

    double t0_out = std::max(outer_hit->t_enter, 0.0);
    double t1_out = outer_hit->t_exit;

    if (inner_hd.minCoeff() > 1e-10) {
        auto inner_hit = intersect_box(local_pos, local_dir,
                                       inner_hd.x(), inner_hd.y(),
                                       -inner_hd.z(), inner_hd.z());
        if (inner_hit && inner_hit->valid()) {
            double path = 0.0;
            double t1_in = inner_hit->t_exit;
            if (t1_in < t1_out) path = t1_out - t1_in;
            double t0_in = std::max(inner_hit->t_enter, 0.0);
            if (t0_out < t0_in) path += t0_in - t0_out;
            return path;
        }
    }
    return t1_out - t0_out;
}

} // namespace

void SourceGeometry::add_shield(const Material* mat, double thickness) {
    assert(mat && thickness > 0.0);
    shields_.push_back({mat, thickness, thickness, thickness});
}

void SourceGeometry::add_shield(const Material* mat, double t_radial, double t_end) {
    assert(mat && t_radial >= 0.0 && t_end >= 0.0 && (t_radial > 0.0 || t_end > 0.0));
    assert(!configured_ || shape_ == Shape::Cylindrical);
    shields_.push_back({mat, t_radial, t_radial, t_end});
}

void SourceGeometry::add_shield(const Material* mat, double t_x, double t_y, double t_z) {
    assert(mat && t_x >= 0.0 && t_y >= 0.0 && t_z >= 0.0
           && (t_x > 0.0 || t_y > 0.0 || t_z > 0.0));
    assert(!configured_ || shape_ == Shape::Rectangular);
    shields_.push_back({mat, t_x, t_y, t_z});
}

void SourceGeometry::configure_point(const Eigen::Vector3d& position) {
    shape_ = Shape::Point;
    point_position_ = position;
    configured_ = true;
}

void SourceGeometry::configure_cylindrical(const Eigen::Vector3d& center, double radius,
                                            double half_length, const Eigen::Matrix3d& rotation,
                                            double inner_radius) {
    assert(inner_radius >= 0.0 && inner_radius < radius);
    shape_ = Shape::Cylindrical;
    cyl_center_ = center;
    cyl_radius_ = radius;
    cyl_inner_r_ = inner_radius;
    cyl_half_length_ = half_length;
    cyl_rotation_ = rotation;
    configured_ = true;
}

void SourceGeometry::configure_spherical(const Eigen::Vector3d& center, double inner_radius,
                                          double outer_radius, const Eigen::Matrix3d& rotation) {
    assert(inner_radius >= 0.0 && inner_radius < outer_radius);
    shape_ = Shape::Sphere;
    sphere_center_ = center;
    sphere_inner_r_ = inner_radius;
    sphere_radius_ = outer_radius;
    sphere_rotation_ = rotation;
    configured_ = true;
}

void SourceGeometry::configure_rectangular(const Eigen::Vector3d& center,
                                            const Eigen::Vector3d& half_dims,
                                            const Eigen::Matrix3d& rotation,
                                            const Eigen::Vector3d& inner_half_dims) {
    // Inner void must be all-zero (solid) or strictly inside the outer box.
    assert((inner_half_dims.array() == 0.0).all()
           || ((inner_half_dims.array() >= 0.0).all()
               && (inner_half_dims.array() < half_dims.array()).all()));
    shape_ = Shape::Rectangular;
    rect_center_ = center;
    rect_half_dims_ = half_dims;
    rect_inner_half_dims_ = inner_half_dims;
    rect_rotation_ = rotation;
    configured_ = true;
}

void SourceGeometry::configure_marinelli(double well_inner_radius, double outer_radius,
                                          double z_bk, double z_we, double z_bot) {
    shape_ = Shape::Marinelli;
    marinelli_well_r_ = well_inner_radius;
    marinelli_outer_r_ = outer_radius;
    marinelli_z_bk_ = z_bk;
    marinelli_z_we_ = z_we;
    marinelli_z_bot_ = z_bot;
    configured_ = true;
}

void SourceGeometry::set_detector_bounds(double det_radius, double det_z_min, double det_z_max) {
    det_bound_radius_ = det_radius;
    det_bound_z_min_ = det_z_min;
    det_bound_z_max_ = det_z_max;
}

bool SourceGeometry::is_inside_marinelli(const Eigen::Vector3d& pos) const {
    double r2 = pos.x() * pos.x() + pos.y() * pos.y();
    double z = pos.z();
    double or2 = marinelli_outer_r_ * marinelli_outer_r_;
    double wr2 = marinelli_well_r_ * marinelli_well_r_;

    if (r2 > or2) return false;
    if (z < marinelli_z_bot_ || z > marinelli_z_we_) return false;

    // In the disk region (z < z_bk): all r ∈ [0, outer_r] are valid
    if (z <= marinelli_z_bk_) return true;

    // In the ring region (z > z_bk): only r ∈ [well_r, outer_r]
    return r2 >= wr2;
}

SourceGeometry::MarinelliExitInfo SourceGeometry::marinelli_exit_info(
    const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) const
{
    // Find nearest exit distance from the L-shaped Marinelli volume,
    // tracking which boundary surface is hit.
    // Boundaries:
    //   1. Outer cylinder: r = outer_r, z ∈ [z_bot, z_we]
    //   2. Well inner wall: r = well_r, z ∈ [z_bk, z_we]
    //   3. Bottom plane: z = z_bot
    //   4. Top of ring: z = z_we (r >= well_r)
    //   5. Well opening: z = z_bk, r < well_r (open — no beaker wall)

    constexpr double kEps = 1e-10;
    double t_min = 1e30;
    MarinelliExitSurface surface = MarinelliExitSurface::None;

    double px = pos.x(), py = pos.y(), pz = pos.z();
    double dx = dir.x(), dy = dir.y(), dz = dir.z();
    double r2 = px * px + py * py;

    // --- Outer cylinder (exit from inside) ---
    {
        double a = dx * dx + dy * dy;
        double b = 2.0 * (px * dx + py * dy);
        double c = r2 - marinelli_outer_r_ * marinelli_outer_r_;
        double disc = b * b - 4.0 * a * c;
        if (disc > 0.0 && a > kEps) {
            double sqrt_disc = std::sqrt(disc);
            double t1 = (-b + sqrt_disc) / (2.0 * a);
            if (t1 > kEps) {
                double z_hit = pz + t1 * dz;
                if (z_hit >= marinelli_z_bot_ - kEps && z_hit <= marinelli_z_we_ + kEps) {
                    if (t1 < t_min) {
                        t_min = t1;
                        surface = MarinelliExitSurface::OuterCylinder;
                    }
                }
            }
        }
    }

    // --- Well inner wall (r = well_r, z ∈ [z_bk, z_we]) ---
    {
        double a = dx * dx + dy * dy;
        double b = 2.0 * (px * dx + py * dy);
        double c = r2 - marinelli_well_r_ * marinelli_well_r_;
        double disc = b * b - 4.0 * a * c;
        if (disc > 0.0 && a > kEps) {
            double sqrt_disc = std::sqrt(disc);
            double t1 = (-b - sqrt_disc) / (2.0 * a);
            double t2 = (-b + sqrt_disc) / (2.0 * a);

            for (double t : {t1, t2}) {
                if (t > kEps) {
                    double z_hit = pz + t * dz;
                    if (z_hit >= marinelli_z_bk_ - kEps && z_hit <= marinelli_z_we_ + kEps) {
                        if (t < t_min) {
                            t_min = t;
                            surface = MarinelliExitSurface::WellWall;
                        }
                        break;
                    }
                }
            }
        }
    }

    // --- Bottom plane: z = z_bot ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_bot_ - pz) / dz;
        if (t > kEps) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            if (r2_hit <= marinelli_outer_r_ * marinelli_outer_r_ + kEps) {
                if (t < t_min) {
                    t_min = t;
                    surface = MarinelliExitSurface::Bottom;
                }
            }
        }
    }

    // --- Top plane: z = z_we (ring region, r >= well_r) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_we_ - pz) / dz;
        if (t > kEps) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            double wr2 = marinelli_well_r_ * marinelli_well_r_;
            if (r2_hit >= wr2 - kEps && r2_hit <= marinelli_outer_r_ * marinelli_outer_r_ + kEps) {
                if (t < t_min) {
                    t_min = t;
                    surface = MarinelliExitSurface::Top;
                }
            }
        }
    }

    // --- Well opening: z = z_bk, r < well_r (open face toward detector) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_bk_ - pz) / dz;
        if (t > kEps) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            double wr2 = marinelli_well_r_ * marinelli_well_r_;
            if (r2_hit < wr2 + kEps) {
                if (t < t_min) {
                    t_min = t;
                    surface = MarinelliExitSurface::WellOpening;
                }
            }
        }
    }

    if (t_min < 1e29) {
        return {t_min, surface};
    }
    return {0.0, MarinelliExitSurface::None};
}

double SourceGeometry::marinelli_boundary_distance(const Eigen::Vector3d& pos,
                                                    const Eigen::Vector3d& dir) const {
    return marinelli_exit_info(pos, dir).distance;
}

double SourceGeometry::marinelli_wall_factor(const Eigen::Vector3d& position,
                                              const Eigen::Vector3d& direction) const {
    if (!is_inside_marinelli(position)) return 0.0;

    auto exit_info = marinelli_exit_info(position, direction);

    // Well opening has no beaker wall
    if (exit_info.surface == MarinelliExitSurface::WellOpening ||
        exit_info.surface == MarinelliExitSurface::None) {
        return 0.0;
    }

    // Compute cosine of incidence angle at exit surface
    double cos_inc = 0.0;
    Eigen::Vector3d exit_pos = position + direction * exit_info.distance;

    switch (exit_info.surface) {
    case MarinelliExitSurface::OuterCylinder:
    case MarinelliExitSurface::WellWall: {
        // Cylindrical surface: normal is radial in xy-plane
        double ex = exit_pos.x(), ey = exit_pos.y();
        double r_exit = std::sqrt(ex * ex + ey * ey);
        if (r_exit > 1e-10) {
            cos_inc = std::abs(direction.x() * ex / r_exit +
                               direction.y() * ey / r_exit);
        }
        break;
    }
    case MarinelliExitSurface::Bottom:
    case MarinelliExitSurface::Top:
        // Flat surface: normal is along z
        cos_inc = std::abs(direction.z());
        break;
    default:
        break;
    }

    // Cap to prevent extreme paths at grazing incidence (max 20x thickness)
    cos_inc = std::max(cos_inc, 0.05);

    return 1.0 / cos_inc;
}

SourceGeometry::MarinelliReentryInfo SourceGeometry::compute_marinelli_reentry(
    const Eigen::Vector3d& pos, const Eigen::Vector3d& dir) const
{
    // Check if a photon from OUTSIDE the Marinelli water can re-enter.
    // Three entry paths:
    //   1. Through well wall → ring water (outward from well cavity)
    //   2. Through well opening → disk water (downward, no wall)
    //   3. Through ring top → ring water (downward from above z_we, no wall)

    constexpr double kEps = 1e-9;
    double t_best = 1e30;
    double wall_path_best = 0.0;

    double px = pos.x(), py = pos.y(), pz = pos.z();
    double dx = dir.x(), dy = dir.y(), dz = dir.z();
    double or2 = marinelli_outer_r_ * marinelli_outer_r_;

    // Beaker wall thickness (from shields — first layer)
    double wall_t = 0.0;
    if (!shields_.empty()) wall_t = shields_[0].scalar_thickness();

    // --- 1. Ring entry through well wall (r = marinelli_well_r_) ---
    // Photon crosses the well outer surface (= water inner boundary) at r = well_r
    {
        double a = dx * dx + dy * dy;
        double b = 2.0 * (px * dx + py * dy);
        double r2_cur = px * px + py * py;
        double c = r2_cur - marinelli_well_r_ * marinelli_well_r_;
        double disc = b * b - 4.0 * a * c;
        if (disc > 0.0 && a > kEps) {
            double sqrt_disc = std::sqrt(disc);
            // We want the FIRST positive intersection (entering the cylinder from inside)
            double t1 = (-b - sqrt_disc) / (2.0 * a);
            double t2 = (-b + sqrt_disc) / (2.0 * a);
            // From inside the well (r < well_r), the exit is the + root
            for (double t : {t1, t2}) {
                if (t > kEps) {
                    double z_hit = pz + t * dz;
                    if (z_hit >= marinelli_z_bk_ - kEps && z_hit <= marinelli_z_we_ + kEps) {
                        // Check that the hit point is within the ring (r_hit ≈ well_r, which is ≥ well_r)
                        if (t < t_best) {
                            t_best = t;
                            // Compute wall path: distance through the cylindrical shell
                            // from r = well_r - wall_t to r = well_r
                            if (wall_t > 0.0) {
                                double r_inner = marinelli_well_r_ - wall_t;
                                double c_inner = r2_cur - r_inner * r_inner;
                                double disc_inner = b * b - 4.0 * a * c_inner;
                                if (disc_inner > 0.0) {
                                    double t_inner = (-b + std::sqrt(disc_inner)) / (2.0 * a);
                                    if (t_inner > kEps && t_inner < t) {
                                        wall_path_best = t - t_inner;
                                    } else {
                                        // Already inside wall
                                        wall_path_best = t;
                                    }
                                } else {
                                    wall_path_best = wall_t; // fallback
                                }
                            } else {
                                wall_path_best = 0.0;
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    // --- 2. Disk entry through well opening (z = z_bk, r < well_r, no wall) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_bk_ - pz) / dz;
        if (t > kEps && t < t_best) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            // Well opening: r < well_r enters disk water (no beaker wall)
            if (r2_hit < marinelli_well_r_ * marinelli_well_r_) {
                // Also verify it's within the outer radius
                if (r2_hit <= or2 + kEps) {
                    t_best = t;
                    wall_path_best = 0.0; // No wall at well opening
                }
            }
        }
    }

    // --- 3. Ring top entry (z = z_we, r ∈ [well_r, outer_r], no wall) ---
    if (std::abs(dz) > kEps) {
        double t = (marinelli_z_we_ - pz) / dz;
        if (t > kEps && t < t_best) {
            double rx = px + t * dx, ry = py + t * dy;
            double r2_hit = rx * rx + ry * ry;
            double wr2 = marinelli_well_r_ * marinelli_well_r_;
            if (r2_hit >= wr2 - kEps && r2_hit <= or2 + kEps) {
                t_best = t;
                wall_path_best = 0.0; // No wall at ring top
            }
        }
    }

    if (t_best < 1e29) {
        // Advance slightly past the boundary into the water
        Eigen::Vector3d entry = pos + dir * (t_best + 1e-8);
        return {true, t_best, wall_path_best, entry};
    }

    return {false, 0.0, 0.0, {0, 0, 0}};
}

double SourceGeometry::source_material_path(const Eigen::Vector3d& position,
                                             const Eigen::Vector3d& direction) const {
    if (!source_material_) return 0.0;

    if (shape_ == Shape::Point) {
        return 0.0; // Point sources have no volume
    }

    if (shape_ == Shape::Cylindrical) {
        // Transform to source local frame
        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);
        Eigen::Vector3d local_dir = cyl_rotation_ * direction;

        // Annulus when a bore is present (subtract the non-attenuating core);
        // cyl_inner_r_ == 0 reduces to the solid-cylinder exit path.
        return cyl_shell_path(local_pos, local_dir,
                              cyl_inner_r_, cyl_half_length_,
                              cyl_radius_, cyl_half_length_);
    }

    if (shape_ == Shape::Sphere) {
        // Shell when a void center is present; inner_r == 0 → solid ball.
        return spherical_shell_path(sphere_center_, position, direction,
                                    sphere_inner_r_, sphere_radius_);
    }

    if (shape_ == Shape::Rectangular) {
        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);
        Eigen::Vector3d local_dir = rect_rotation_ * direction;

        // Hollow box shell: subtract ray overlap with the non-attenuating
        // inner void (a ray born in the shell can cross the void and re-enter
        // the far wall). Solid boxes keep the original exit-path expression
        // below, unchanged.
        if (rect_inner_half_dims_.minCoeff() > 1e-10) {
            return box_shell_path(local_pos, local_dir,
                                  rect_inner_half_dims_, rect_half_dims_);
        }

        auto hit = intersect_box(local_pos, local_dir,
                                 rect_half_dims_.x(), rect_half_dims_.y(),
                                 -rect_half_dims_.z(), rect_half_dims_.z());
        if (!hit || !hit->valid()) return 0.0;

        return std::max(hit->t_exit, 0.0);
    }

    if (shape_ == Shape::Marinelli) {
        if (!is_inside_marinelli(position)) return 0.0;
        return marinelli_boundary_distance(position, direction);
    }

    return 0.0;
}

double SourceGeometry::min_distance_to_boundary(
    const Eigen::Vector3d& position, bool include_shields) const {
    // Envelope = source-material extent, plus (when include_shields) the
    // summed shield thicknesses per axis. For a convex envelope the nearest
    // boundary point from an interior point is the perpendicular foot on the
    // closest face, so the per-axis min below is exact (never an over-estimate).
    switch (shape_) {
    case Shape::Point: {
        // Point source: no material volume — electrons originate in the
        // spherical shield shells, which fill outward from the point.
        if (!include_shields) return 0.0;
        double r_out = 0.0;
        for (const auto& layer : shields_) r_out += layer.scalar_thickness();
        double r = (position - point_position_).norm();
        return std::max(r_out - r, 0.0);
    }
    case Shape::Cylindrical: {
        double r_out = cyl_radius_;
        double l_out = cyl_half_length_;
        if (include_shields)
            for (const auto& layer : shields_) { r_out += layer.tx; l_out += layer.tz; }
        Eigen::Vector3d local = cyl_rotation_ * (position - cyl_center_);
        double rho = std::hypot(local.x(), local.y());
        double d = std::min(r_out - rho, l_out - std::abs(local.z()));
        // Hollow bore: the inner void wall is also a material boundary when we
        // measure to the source-material surface (include_shields == false).
        if (!include_shields && cyl_inner_r_ > 1e-10)
            d = std::min(d, rho - cyl_inner_r_);
        return std::max(d, 0.0);
    }
    case Shape::Sphere: {
        double r_out = sphere_radius_;
        if (include_shields)
            for (const auto& layer : shields_) r_out += layer.scalar_thickness();
        double r = (position - sphere_center_).norm();
        double d = r_out - r;
        if (!include_shields && sphere_inner_r_ > 1e-10)
            d = std::min(d, r - sphere_inner_r_);
        return std::max(d, 0.0);
    }
    case Shape::Rectangular: {
        Eigen::Vector3d h = rect_half_dims_;
        if (include_shields)
            for (const auto& layer : shields_) {
                h.x() += layer.tx; h.y() += layer.ty; h.z() += layer.tz;
            }
        Eigen::Vector3d local = rect_rotation_ * (position - rect_center_);
        double dx = h.x() - std::abs(local.x());
        double dy = h.y() - std::abs(local.y());
        double dz = h.z() - std::abs(local.z());
        double d = std::min({dx, dy, dz});
        // Hollow shell: the inner-void wall is also a material boundary when
        // we measure to the source-material surface (include_shields == false).
        // For each axis where the point lies beyond the inner box, the
        // perpendicular foot on that inner face is a candidate. Conservative
        // (<= true distance) since the shell is non-convex — same caveat as
        // the cylinder bore.
        if (!include_shields && rect_inner_half_dims_.minCoeff() > 1e-10) {
            for (int i = 0; i < 3; ++i) {
                const double beyond = std::abs(local[i]) - rect_inner_half_dims_[i];
                if (beyond > 0.0)
                    d = std::min(d, beyond);
            }
        }
        return std::max(d, 0.0);
    }
    case Shape::Marinelli: {
        // Conservative (≤ true) distance from an interior water point to the
        // nearest water boundary (beaker wall or re-entrant well void) of the
        // L-shaped fill.  Under-estimate is safe for the containment test; the
        // beaker walls ARE the boundary, so include_shields does not extend it.
        // Water = disk (z∈[z_bot,z_bk], r∈[0,outer]) ∪ ring (z∈[z_bk,z_we],
        // r∈[well,outer]).
        const double rho = std::hypot(position.x(), position.y());
        const double z   = position.z();
        double d = marinelli_outer_r_ - rho;            // outer beaker wall
        d = std::min(d, z - marinelli_z_bot_);          // beaker bottom
        if (z >= marinelli_z_bk_) {                      // ring region
            d = std::min(d, rho - marinelli_well_r_);   // well inner wall (void)
            d = std::min(d, marinelli_z_we_ - z);       // ring top
        } else if (rho < marinelli_well_r_) {           // disk under the well
            d = std::min(d, marinelli_z_bk_ - z);       // well floor (void above)
        }
        return std::max(d, 0.0);
    }
    default:
        return 0.0;
    }
}

double SourceGeometry::compute_transmission(const Eigen::Vector3d& position,
                                             const Eigen::Vector3d& direction,
                                             double energy_MeV) const {
    double total_mu_L = 0.0;

    // Source material attenuation
    if (source_material_) {
        double path = source_material_path(position, direction);
        if (path > 0.0) {
            total_mu_L += source_material_->mu_total(energy_MeV) * path;
        }
    }

    // Shielding layers
    if (shape_ == Shape::Point) {
        // Spherical shells: path = thickness (direction-independent)
        for (const auto& layer : shields_) {
            assert(layer.is_uniform());
            total_mu_L += layer.material->mu_total(energy_MeV) * layer.scalar_thickness();
        }
    } else if (shape_ == Shape::Cylindrical) {
        // Cylindrical shells: compute path through each shell
        double inner_r = cyl_radius_;
        double inner_half_z = cyl_half_length_;

        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);
        Eigen::Vector3d local_dir = cyl_rotation_ * direction;

        for (const auto& layer : shields_) {
            assert(layer.tx == layer.ty);
            double outer_r = inner_r + layer.tx;
            double outer_half_z = inner_half_z + layer.tz;

            auto outer_hit = intersect_cylinder(local_pos, local_dir,
                                                outer_r, -outer_half_z, outer_half_z);
            auto inner_hit = intersect_cylinder(local_pos, local_dir,
                                                inner_r, -inner_half_z, inner_half_z);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;

                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    // Path from inner exit to outer exit
                    if (t1_in < t1_out) {
                        path = t1_out - t1_in;
                    }
                    // Also front segment if applicable
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) {
                        path += t0_in - t0_out;
                    }
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 0.0) {
                total_mu_L += layer.material->mu_total(energy_MeV) * path;
            }

            inner_r = outer_r;
            inner_half_z = outer_half_z;
        }
    } else if (shape_ == Shape::Rectangular) {
        double inner_hx = rect_half_dims_.x();
        double inner_hy = rect_half_dims_.y();
        double inner_hz = rect_half_dims_.z();

        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);
        Eigen::Vector3d local_dir = rect_rotation_ * direction;

        for (const auto& layer : shields_) {
            double outer_hx = inner_hx + layer.tx;
            double outer_hy = inner_hy + layer.ty;
            double outer_hz = inner_hz + layer.tz;

            auto outer_hit = intersect_box(local_pos, local_dir,
                                           outer_hx, outer_hy, -outer_hz, outer_hz);
            auto inner_hit = intersect_box(local_pos, local_dir,
                                           inner_hx, inner_hy, -inner_hz, inner_hz);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;

                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) {
                        path = t1_out - t1_in;
                    }
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) {
                        path += t0_in - t0_out;
                    }
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 0.0) {
                total_mu_L += layer.material->mu_total(energy_MeV) * path;
            }

            inner_hx = outer_hx;
            inner_hy = outer_hy;
            inner_hz = outer_hz;
        }
    } else if (shape_ == Shape::Sphere) {
        // Spherical shells: ray-trace each shell (the emission point is offset
        // from the center, so the path is direction-dependent — unlike a point
        // source emitting from the exact center).
        double inner_r = sphere_radius_;
        for (const auto& layer : shields_) {
            assert(layer.is_uniform());
            double outer_r = inner_r + layer.scalar_thickness();
            double path = spherical_shell_path(sphere_center_, position, direction,
                                               inner_r, outer_r);
            if (path > 0.0)
                total_mu_L += layer.material->mu_total(energy_MeV) * path;
            inner_r = outer_r;
        }
    } else if (shape_ == Shape::Marinelli) {
        // Beaker wall: path depends on exit surface. Well opening has no wall.
        double wall_factor = marinelli_wall_factor(position, direction);
        if (wall_factor > 0.0) {
            for (const auto& layer : shields_) {
                assert(layer.is_uniform());
                total_mu_L += layer.material->mu_total(energy_MeV) * layer.scalar_thickness() * wall_factor;
            }
        }
    }

    return std::exp(-total_mu_L);
}

double SourceGeometry::point_source_transmission(double energy_MeV) const {
    double total_mu_L = 0.0;
    for (const auto& layer : shields_) {
        assert(layer.is_uniform());
        total_mu_L += layer.material->mu_total(energy_MeV) * layer.scalar_thickness();
    }
    return std::exp(-total_mu_L);
}

SourceGeometry::SourceTransmissionResult SourceGeometry::compute_transmission_fep_only(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_MeV,
    std::mt19937_64& rng) const
{
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double energy_keV = energy_MeV * 1e3;
    double weight = 1.0;
    Eigen::Vector3d dir = direction;
    double total_path = 0.0;  // accumulate path for exit_position

    // Helper: process a single material layer (path length known)
    auto process_layer = [&](const Material* mat, double path) {
        if (!mat || path <= 0.0) return;
        total_path += path;
        MacroscopicXS xs = mat->macroscopic_xs(energy_MeV);
        double mu_no_rs = xs.mu_pe + xs.mu_cs + xs.mu_pp;
        double mu_rs = xs.mu_rs;

        weight *= std::exp(-mu_no_rs * path);

        // Stochastic Rayleigh
        double lambda_rs = mu_rs * path;
        int n_scatter = 0;
        if (lambda_rs > 0.01) {
            std::poisson_distribution<int> poisson(lambda_rs);
            n_scatter = poisson(rng);
        } else if (lambda_rs > 1e-10) {
            n_scatter = (uniform(rng) < lambda_rs) ? 1 : 0;
        }

        if (n_scatter > 0) {
            int Z_rs = mat->select_element(energy_MeV, 2 /*RS*/, uniform(rng));
            for (int i = 0; i < n_scatter; ++i) {
                double cos_theta = sample_rayleigh_cos_theta(Z_rs, energy_keV, rng);
                double phi = kTwoPi * uniform(rng);
                dir = rotate_direction(dir, cos_theta, phi);
            }
        }
    };

    // Source material attenuation
    if (source_material_) {
        double path = source_material_path(position, direction);
        process_layer(source_material_, path);
    }

    // Shielding layers
    if (shape_ == Shape::Point) {
        for (const auto& layer : shields_) {
            assert(layer.is_uniform());
            process_layer(layer.material, layer.scalar_thickness());
        }
    } else if (shape_ == Shape::Cylindrical) {
        double inner_r = cyl_radius_;
        double inner_half_z = cyl_half_length_;
        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);

        for (const auto& layer : shields_) {
            // Recompute local_dir each iteration: process_layer may change dir
            // via Rayleigh scattering in an earlier layer.
            Eigen::Vector3d local_dir = cyl_rotation_ * dir;

            assert(layer.tx == layer.ty);
            double outer_r = inner_r + layer.tx;
            double outer_half_z = inner_half_z + layer.tz;

            auto outer_hit = intersect_cylinder(local_pos, local_dir,
                                                outer_r, -outer_half_z, outer_half_z);
            auto inner_hit = intersect_cylinder(local_pos, local_dir,
                                                inner_r, -inner_half_z, inner_half_z);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            process_layer(layer.material, path);

            inner_r = outer_r;
            inner_half_z = outer_half_z;
        }
    } else if (shape_ == Shape::Rectangular) {
        double inner_hx = rect_half_dims_.x();
        double inner_hy = rect_half_dims_.y();
        double inner_hz = rect_half_dims_.z();
        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);

        for (const auto& layer : shields_) {
            // Recompute local_dir each iteration: process_layer may change dir
            // via Rayleigh scattering in an earlier layer.
            Eigen::Vector3d local_dir = rect_rotation_ * dir;

            double outer_hx = inner_hx + layer.tx;
            double outer_hy = inner_hy + layer.ty;
            double outer_hz = inner_hz + layer.tz;

            auto outer_hit = intersect_box(local_pos, local_dir,
                                           outer_hx, outer_hy, -outer_hz, outer_hz);
            auto inner_hit = intersect_box(local_pos, local_dir,
                                           inner_hx, inner_hy, -inner_hz, inner_hz);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            process_layer(layer.material, path);

            inner_hx = outer_hx;
            inner_hy = outer_hy;
            inner_hz = outer_hz;
        }
    } else if (shape_ == Shape::Sphere) {
        // Spherical shells: ray-trace each shell with the current direction
        // (Rayleigh in an earlier layer may have changed `dir`).
        double inner_r = sphere_radius_;
        for (const auto& layer : shields_) {
            assert(layer.is_uniform());
            double outer_r = inner_r + layer.scalar_thickness();
            double path = spherical_shell_path(sphere_center_, position, dir,
                                               inner_r, outer_r);
            process_layer(layer.material, path);
            inner_r = outer_r;
        }
    } else if (shape_ == Shape::Marinelli) {
        // Beaker wall: path depends on exit surface. Well opening has no wall.
        double wall_factor = marinelli_wall_factor(position, direction);
        if (wall_factor > 0.0) {
            for (const auto& layer : shields_) {
                assert(layer.is_uniform());
                process_layer(layer.material, layer.scalar_thickness() * wall_factor);
            }
        }
    }

    return {weight, dir, position + dir * total_path};
}

std::vector<SourceGeometry::SourcePathSegment> SourceGeometry::trace_source_segments(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV) const
{
    std::vector<SourcePathSegment> segments;
    trace_source_segments(position, direction, energy_keV, segments);
    return segments;
}

void SourceGeometry::trace_source_segments(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double /*energy_keV*/,
    std::vector<SourcePathSegment>& segments,
    std::size_t max_segments) const
{
    segments.clear();
    segments.reserve(1 + shields_.size());

    // Source material (extended sources only)
    if (source_material_ && shape_ != Shape::Point) {
        double path = source_material_path(position, direction);
        if (path > 1e-10) {
            segments.push_back({source_material_, path});
            if (segments.size() >= max_segments) return;
        }
    }

    // Shielding layers
    if (shape_ == Shape::Point || shape_ == Shape::Sphere) {
        // Concentric spherical shells. Point: shells grow from the point
        // (inner_r starts at 0). Sphere: shells wrap the source ball (inner_r
        // starts at the source outer radius). Both use the shared shell-path
        // helper (the emission point may be offset from the center for a
        // spherical volume source, so the path is direction-dependent).
        const Eigen::Vector3d& center =
            (shape_ == Shape::Point) ? point_position_ : sphere_center_;
        double inner_r = (shape_ == Shape::Point) ? 0.0 : sphere_radius_;
        for (const auto& layer : shields_) {
            assert(layer.is_uniform());
            double outer_r = inner_r + layer.scalar_thickness();
            double path = spherical_shell_path(center, position, direction,
                                               inner_r, outer_r);
            if (path > 1e-10) {
                segments.push_back({layer.material, path});
                if (segments.size() >= max_segments) return;
            }
            inner_r = outer_r;
        }
    } else if (shape_ == Shape::Cylindrical) {
        double inner_r = cyl_radius_;
        double inner_half_z = cyl_half_length_;
        Eigen::Vector3d local_pos = cyl_rotation_ * (position - cyl_center_);
        Eigen::Vector3d local_dir = cyl_rotation_ * direction;

        for (const auto& layer : shields_) {
            assert(layer.tx == layer.ty);
            double outer_r = inner_r + layer.tx;
            double outer_half_z = inner_half_z + layer.tz;

            auto outer_hit = intersect_cylinder(local_pos, local_dir,
                                                outer_r, -outer_half_z, outer_half_z);
            auto inner_hit = intersect_cylinder(local_pos, local_dir,
                                                inner_r, -inner_half_z, inner_half_z);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 1e-10) {
                segments.push_back({layer.material, path});
                if (segments.size() >= max_segments) return;
            }

            inner_r = outer_r;
            inner_half_z = outer_half_z;
        }
    } else if (shape_ == Shape::Rectangular) {
        double inner_hx = rect_half_dims_.x();
        double inner_hy = rect_half_dims_.y();
        double inner_hz = rect_half_dims_.z();
        Eigen::Vector3d local_pos = rect_rotation_ * (position - rect_center_);
        Eigen::Vector3d local_dir = rect_rotation_ * direction;

        for (const auto& layer : shields_) {
            double outer_hx = inner_hx + layer.tx;
            double outer_hy = inner_hy + layer.ty;
            double outer_hz = inner_hz + layer.tz;

            auto outer_hit = intersect_box(local_pos, local_dir,
                                           outer_hx, outer_hy, -outer_hz, outer_hz);
            auto inner_hit = intersect_box(local_pos, local_dir,
                                           inner_hx, inner_hy, -inner_hz, inner_hz);

            double path = 0.0;
            if (outer_hit && outer_hit->valid()) {
                double t0_out = std::max(outer_hit->t_enter, 0.0);
                double t1_out = outer_hit->t_exit;
                if (inner_hit && inner_hit->valid()) {
                    double t1_in = inner_hit->t_exit;
                    if (t1_in < t1_out) path = t1_out - t1_in;
                    double t0_in = std::max(inner_hit->t_enter, 0.0);
                    if (t0_out < t0_in) path += t0_in - t0_out;
                } else {
                    path = t1_out - t0_out;
                }
            }

            if (path > 1e-10) {
                segments.push_back({layer.material, path});
                if (segments.size() >= max_segments) return;
            }

            inner_hx = outer_hx;
            inner_hy = outer_hy;
            inner_hz = outer_hz;
        }
    } else if (shape_ == Shape::Marinelli) {
        // Beaker wall: path depends on exit surface. Well opening has no wall.
        // Only add wall if photon is inside sample volume.
        double wall_factor = marinelli_wall_factor(position, direction);
        if (wall_factor > 0.0) {
            for (const auto& layer : shields_) {
                assert(layer.is_uniform());
                double path = layer.scalar_thickness() * wall_factor;
                if (path > 1e-10) {
                    segments.push_back({layer.material, path});
                    if (segments.size() >= max_segments) return;
                }
            }
        }
    }
}

double SourceGeometry::no_interaction_probability(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV,
    double* path_cm) const
{
    const auto segments = trace_source_segments(position, direction, energy_keV);
    const double energy_MeV = energy_keV * 1e-3;
    double tau = 0.0;
    double path = 0.0;
    for (const auto& seg : segments) {
        if (!seg.material || seg.length <= 0.0) continue;
        path += seg.length;
        tau += seg.material->macroscopic_xs(energy_MeV).mu_total() * seg.length;
    }
    if (path_cm) *path_cm = path;
    return std::exp(-tau);
}

SourceGeometry::SourceFullTransportResult SourceGeometry::transport_source_photon(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV,
    std::mt19937_64& rng) const
{
    // Compton-angle biasing applies only to the primary photon (this public
    // entry point); annihilation/brems recursion always runs analog. Never
    // biased for Marinelli (re-entry transports share this entry point and
    // the auto-policy never enables it there).
    const int budget = (compton_bias_.cone_fraction > 0.0 &&
                        shape_ != Shape::Marinelli)
                           ? compton_bias_.max_biased_vertices
                           : 0;
    return transport_source_photon_impl(position, direction, energy_keV, rng,
                                        budget);
}

SourceGeometry::SourceFullTransportResult SourceGeometry::transport_source_photon_impl(
    const Eigen::Vector3d& position,
    const Eigen::Vector3d& direction,
    double energy_keV,
    std::mt19937_64& rng,
    int bias_budget) const
{
    constexpr int kMaxSourceScatters = 50;
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    Eigen::Vector3d pos = position;
    Eigen::Vector3d dir = direction;
    double energy = energy_keV;

    // Accumulate source electrons and secondary photons (annihilation gammas,
    // bremsstrahlung) across multiple Compton scatters
    std::vector<SourceElectron> result_electrons;
    std::vector<SourceSecondaryPhoton> result_secondaries;

    // True once any interaction (Compton/PE/PP/Rayleigh) has occurred,
    // including scatters in earlier re-trace iterations.
    bool any_interaction = false;

    // Accumulated Compton-angle importance-sampling weight (1.0 = analog).
    double w_bias = 1.0;

    // Transport a brems photon through the remaining source geometry and keep
    // the survivor (and any nested secondaries) as event secondaries.
    auto transport_brems_photon = [&](const BremsPhoton& bp) {
        auto br = transport_source_photon_impl(
            bp.position, bp.direction, bp.energy_keV, rng, /*bias_budget=*/0);
        if (br.survived)
            result_secondaries.push_back(
                {br.position, br.direction, br.energy_keV});
        for (auto& sec : br.secondaries)
            result_secondaries.push_back(std::move(sec));
        for (auto& e : br.source_electrons)
            result_electrons.push_back(std::move(e));
    };

    // Detect whether the source geometry is a single material (no internal
    // interfaces). Computed once per primary.
    const Material* first_mat = source_material_;
    bool single_material = true;
    for (const auto& layer : shields_) {
        if (!layer.material) continue;
        if (!first_mat) first_mat = layer.material;
        else if (layer.material != first_mat) { single_material = false; break; }
    }

    // Handle an electron (Compton recoil, photoelectron, pair e+/e-) born in
    // the source geometry at e_pos in material `mat`.
    //
    // ALL shapes (incl. Marinelli as of the regime-aware escape work): Moliere
    // condensed-history walk through the source segments — emits bremsstrahlung
    // along the way and detects escape through the outer boundary (the direction
    // diffusion is essential in dense shields; straight-line CSDA overestimates
    // escape several-fold).  Marinelli previously used the legacy straight-line
    // try_source_electron for escape; it now shares this walk + the regime-aware
    // exit-state skin-escape gate (light water → fully analog), so cfg 8 joins
    // the gated/validated set.
    //
    // Perf: the geometry-tracing Moliere walk dominates the cost of the
    // secondary channels — it re-traces the source segments at every substep
    // (profiling: trace_source_segments is the second-hottest leaf). Most
    // electrons are born deep enough that they stop without leaving their
    // birth material. For those we (a) cannot escape and (b) can sample brems
    // with the geometry-free ElectronCsda::sample_brems_in_material (the same
    // SB physics, creation-point approximation — exact when the range is small
    // vs the medium, which is precisely the contained case; this is the
    // sampler already used for Marinelli). The walk is reserved for electrons
    // that can reach a material interface or escape.
    //
    // Containment test: the electron stays within material `mat` if its CSDA
    // range in `mat` is below the distance to the nearest surface of the `mat`
    // region. `d_stay` always measures to the boundary of the birth material's
    // own region, and containment means the electron never leaves `mat`, so the
    // matching reach bound is `mat`'s own range/density — exact, not merely
    // conservative (using any other material's range/density would be a
    // mismatch, since the g/cm^2 CSDA range is itself Z/A-dependent). The
    // per-config escape-energy floor (set_source_electron_threshold, default
    // 50 keV) is the other knob.
    auto process_electron = [&](const Eigen::Vector3d& e_pos,
                                const Eigen::Vector3d& e_dir,
                                double T_keV, const Material& mat) {
        const bool want_brems = source_brems_enabled_ &&
            T_keV > kMoliereBremsThreshold_keV;
        const bool want_escape = source_electron_enabled_ &&
            T_keV > source_electron_threshold_keV_;
        if (!want_brems && !want_escape) return;

        // The electron's reach (cm) in its birth material `mat`: the straight-
        // line displacement cannot exceed the CSDA path length. Exact for the
        // containment test below, which only ever asks whether the electron
        // stays inside the `mat` region.
        double range_cm = 0.0;
        if (mat.density() > 0.0) {
            range_cm =
                ElectronCsda::instance().range_gcm2_material(mat, T_keV)
                / mat.density();
        }

        // Distance the electron is guaranteed to remain inside material `mat`:
        // the outer boundary for a single-material geometry, otherwise the
        // source-material surface for a source-material electron. Shield-born
        // electrons in a multi-material geometry take the full walk (0 here).
        double d_stay = 0.0;
        if (single_material)
            d_stay = min_distance_to_boundary(e_pos, /*include_shields=*/true);
        else if (&mat == source_material_)
            d_stay = min_distance_to_boundary(e_pos, /*include_shields=*/false);

        if (range_cm < d_stay) {
            // Contained in `mat`: cannot escape, cannot cross an interface.
            if (want_brems) {
                auto brems = ElectronCsda::instance().sample_brems_in_material(
                    mat, e_pos, e_dir, T_keV, rng);
                for (const auto& bp : brems) transport_brems_photon(bp);
            }
            return;  // no escaped electron possible
        }

        auto walk = ElectronCsda::instance().walk_in_source_geometry(
            *this, mat, e_pos, e_dir, T_keV, rng);

        if (source_brems_enabled_)
            for (const auto& bp : walk.brems_photons) transport_brems_photon(bp);

        if (source_electron_enabled_ && walk.escaped &&
            walk.exit_KE_keV > source_electron_threshold_keV_) {
            result_electrons.push_back(
                {walk.exit_position, walk.exit_direction, walk.exit_KE_keV});
        }
    };

    std::vector<SourcePathSegment> segments;  // reused across retrace iters (E2)
    for (int retrace = 0; retrace < kMaxSourceScatters; ++retrace) {
        trace_source_segments(pos, dir, energy, segments);

        if (segments.empty()) {
            return {true, pos, dir, energy, std::move(result_secondaries),
                    std::move(result_electrons), any_interaction, w_bias};
        }

        bool scattered = false;
        double distance_along_ray = 0.0;

        for (const auto& seg : segments) {
            double energy_MeV = energy * 1e-3;
            MacroscopicXS xs = seg.material->macroscopic_xs(energy_MeV);
            double mu_total = xs.mu_total();

            if (mu_total <= 0.0) {
                distance_along_ray += seg.length;
                continue;
            }

            double s = sample_interaction_distance(mu_total, rng);

            if (s >= seg.length) {
                // No interaction in this segment — advance past it
                distance_along_ray += seg.length;
                continue;
            }

            // Interaction occurs at distance s within this segment
            pos += dir * (distance_along_ray + s);

            InteractionType type = select_interaction(xs, uniform(rng));
            any_interaction = true;

            switch (type) {
            case InteractionType::Photoelectric: {
                // Photoelectron carries (nearly) the full photon energy;
                // binding energy is negligible against the 200 keV brems floor.
                if (energy > kMoliereBremsThreshold_keV) {
                    Eigen::Vector3d e_dir =
                        sample_photoelectron_direction(dir, energy, rng);
                    process_electron(pos, e_dir, energy, *seg.material);
                }
                return {false, pos, dir, energy, std::move(result_secondaries),
                        std::move(result_electrons), true, w_bias};
            }

            case InteractionType::PairProduction: {
                // Generate two back-to-back annihilation gammas (electron rest mass energy)
                // and transport each through the remaining source material.
                constexpr double kElectronMassKeV = 510.998950;

                // The e+/e- pair stopping in the shield: bremsstrahlung and
                // possible escape. Pair kinetic energy split sampled
                // uniformly; both particles start (approximately) along the
                // photon direction.
                double pair_T = energy - 2.0 * kElectronMassKeV;
                if (pair_T > 0.0) {
                    double u = uniform(rng);
                    process_electron(pos, dir, u * pair_T, *seg.material);
                    process_electron(pos, dir, (1.0 - u) * pair_T, *seg.material);
                }

                SourceFullTransportResult pp_result{false, pos, dir, energy,
                                                    std::move(result_secondaries),
                                                    std::move(result_electrons),
                                                    true, w_bias};
                Eigen::Vector3d ann_dir1 = sample_isotropic_direction(rng);
                Eigen::Vector3d ann_dir2 = -ann_dir1;
                for (const auto& ann_dir : {ann_dir1, ann_dir2}) {
                    auto ann = transport_source_photon_impl(
                        pos, ann_dir, kElectronMassKeV, rng, /*bias_budget=*/0);
                    if (ann.survived)
                        pp_result.secondaries.push_back(
                            {ann.position, ann.direction, ann.energy_keV});
                    // Propagate secondaries and source electrons from recursive calls
                    for (auto& sec : ann.secondaries)
                        pp_result.secondaries.push_back(std::move(sec));
                    for (auto& e : ann.source_electrons)
                        pp_result.source_electrons.push_back(std::move(e));
                }
                return pp_result;
            }

            case InteractionType::Compton: {
                int Z_cs = seg.material->select_element(energy_MeV, 1 /*CS*/, uniform(rng));

                // Compton-angle mixture biasing (primary's first K vertices):
                // sample the outgoing direction from
                //   q(w) = (1-gamma) p(w) + gamma * 1[w in cone]/Omega_d,
                // p(w) = f(mu)/(2pi*N) the analog KN x S density on the
                // sphere, cone subtending the detector bounding sphere from
                // this vertex; multiply the event weight by p/q <= 1/(1-gamma).
                // Energy/recoil at the sampled angle go through the exact
                // analog kinematics (finish_compton_at_angle).
                ComptonScatterResult scatter;
                bool biased_vertex = false;
                if (bias_budget > 0 && det_bound_radius_ > 0.0) {
                    const double gamma = compton_bias_.cone_fraction;
                    const Eigen::Vector3d det_c(
                        0.0, 0.0, 0.5 * (det_bound_z_min_ + det_bound_z_max_));
                    const double half_h =
                        0.5 * (det_bound_z_max_ - det_bound_z_min_);
                    const double Rb = 1.05 * std::sqrt(
                        det_bound_radius_ * det_bound_radius_ + half_h * half_h);
                    Eigen::Vector3d axis = det_c - pos;
                    const double dist = axis.norm();
                    if (dist > Rb) {
                        axis /= dist;
                        const double sin_half = Rb / dist;
                        const double cos_half =
                            std::sqrt(1.0 - sin_half * sin_half);
                        const double omega_d = kTwoPi * (1.0 - cos_half);

                        Eigen::Vector3d w_dir;
                        if (uniform(rng) < gamma) {
                            double ct = 1.0 - uniform(rng) * (1.0 - cos_half);
                            double ph = kTwoPi * uniform(rng);
                            w_dir = rotate_direction(axis, ct, ph);
                        } else {
                            double mu_a =
                                sample_compton_cos_theta(energy, Z_cs, rng);
                            double ph = kTwoPi * uniform(rng);
                            w_dir = rotate_direction(dir, mu_a, ph);
                        }
                        const double mu = std::max(
                            -1.0, std::min(1.0, w_dir.dot(dir)));
                        const double p_sphere =
                            compton_angular_pdf_unnorm(energy, Z_cs, mu) /
                            (kTwoPi * compton_angular_norm(energy, Z_cs));
                        const double c_sphere =
                            (w_dir.dot(axis) >= cos_half) ? 1.0 / omega_d : 0.0;
                        const double q =
                            (1.0 - gamma) * p_sphere + gamma * c_sphere;
                        if (q > 0.0) {
                            w_bias *= p_sphere / q;
                            scatter = finish_compton_at_angle(
                                energy, dir, mu, rng, Z_cs,
                                /*doppler=*/false, &w_dir);
                            biased_vertex = true;
                            --bias_budget;
                        }
                    }
                }
                if (!biased_vertex) {
                    scatter = sample_compton_scatter(energy, dir, rng, Z_cs);
                }

                // Recoil electron: bremsstrahlung (all shapes) and Moliere
                // escape (non-Marinelli) via process_electron.
                double recoil_T = energy - scatter.scattered_energy_keV;
                if (recoil_T > std::min(kMoliereBremsThreshold_keV,
                                        source_electron_threshold_keV_)) {
                    Eigen::Vector3d e_dir = compton_recoil_direction(
                        dir, scatter.new_direction,
                        energy, scatter.scattered_energy_keV);
                    process_electron(pos, e_dir, recoil_T, *seg.material);
                }

                energy = scatter.scattered_energy_keV;
                dir = scatter.new_direction;
                if (energy < 1.0) {
                    return {false, pos, dir, energy, std::move(result_secondaries),
                            std::move(result_electrons), true, w_bias};
                }
                scattered = true;
                break;
            }

            case InteractionType::Rayleigh: {
                int Z_rs = seg.material->select_element(energy_MeV, 2 /*RS*/, uniform(rng));
                double cos_theta = sample_rayleigh_cos_theta(Z_rs, energy, rng);
                double phi = kTwoPi * uniform(rng);
                dir = rotate_direction(dir, cos_theta, phi);
                scattered = true;
                break;
            }
            }

            if (scattered) break;  // Re-trace from new position/direction
        }

        if (!scattered) {
            // Walked all segments without interaction — photon exited
            pos += dir * distance_along_ray;
            return {true, pos, dir, energy, std::move(result_secondaries),
                    std::move(result_electrons), any_interaction, w_bias};
        }
    }

    // Hit max scatters — kill the photon
    return {false, pos, dir, energy, std::move(result_secondaries),
            std::move(result_electrons), true, w_bias};
}

double SourceGeometry::outermost_extent_radius() const {
    double r = 0.0;

    if (shape_ == Shape::Point) {
        // Total shield thickness
        for (const auto& layer : shields_) {
            r += layer.scalar_thickness();
        }
    } else if (shape_ == Shape::Cylindrical) {
        // Bounding-sphere radius: shields grow radially by tx and axially by tz.
        double rr = cyl_radius_;
        double hz = cyl_half_length_;
        for (const auto& layer : shields_) {
            rr += layer.tx;
            hz += layer.tz;
        }
        r = std::sqrt(rr * rr + hz * hz);
    } else if (shape_ == Shape::Rectangular) {
        // Bounding-sphere radius: half-diagonal of the outermost box.
        Eigen::Vector3d hd = rect_half_dims_;
        for (const auto& layer : shields_) {
            hd += Eigen::Vector3d(layer.tx, layer.ty, layer.tz);
        }
        r = hd.norm();
    } else if (shape_ == Shape::Sphere) {
        r = sphere_radius_;
        for (const auto& layer : shields_) {
            r += layer.scalar_thickness();
        }
    } else if (shape_ == Shape::Marinelli) {
        r = marinelli_outer_r_;
        for (const auto& layer : shields_) {
            r += layer.scalar_thickness();
        }
    }

    return r;
}

} // namespace ceelo
