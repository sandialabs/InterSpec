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

#include "geometry/RayTrace.h"
#include "geometry/Cylinder.h"
#include "geometry/Box.h"
#include "materials/Material.h"
#include "transport/ComptonScatter.h"

#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits>
#include <random>

namespace ceelo {

namespace {
constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
} // anonymous namespace

// ---- Geometry member implementations for ray tracing ----

Geometry::Geometry() = default;

void Geometry::set_detector(DetectorShape shape, const Material* material,
                            const std::vector<double>& dimensions) {
    shape_ = shape;
    detector_material_ = material;

    // Redefining the crystal invalidates everything sized against the old one:
    // a fillet radius, bore or dead layer that no longer fits would otherwise
    // survive silently (a stale fillet leaves rho_c <= 0 with no diagnostic in
    // any build). Clear them all and make the caller re-declare what it wants.
    bullet_radius_ = 0.0;
    bore_hole_.reset();
    dead_layer_.reset();

    if (shape == DetectorShape::Cylinder) {
        assert(dimensions.size() >= 2);
        radius_ = dimensions[0];
        length_ = dimensions[1];
    } else {
        assert(dimensions.size() >= 3);
        half_x_ = dimensions[0];
        half_y_ = dimensions[1];
        length_ = dimensions[2];
    }
}

/// The bore must stay inside the crystal wall over its whole depth.  With a
/// bulletized front edge the crystal narrows to `radius_ - bullet_radius` near
/// the front face, so a deep enough bore can be admissible against the full
/// radius and still be wider than the crystal where it ends.
static bool bore_fits(double bore_radius, double bore_depth,
                      double radius, double length, double bullet_radius) {
    if (bore_radius >= radius) return false;
    if (bullet_radius <= 0.0) return true;
    const double z_bore_start = length - bore_depth;   // apex, from the front
    if (z_bore_start >= bullet_radius) return true;    // ends behind the fillet
    // Crystal radius at the bore's apex, on the fillet arc.
    const double rho_c = radius - bullet_radius;
    const double dz = bullet_radius - z_bore_start;
    const double r_here = rho_c + std::sqrt(std::max(
        0.0, bullet_radius * bullet_radius - dz * dz));
    return bore_radius < r_here;
}

void Geometry::set_bore_hole(double bore_radius, double bore_depth,
                             bool rounded_tip) {
    assert(shape_ == DetectorShape::Cylinder && "Bore hole only for cylindrical detectors");
    assert(bore_radius > 0.0 && bore_radius < radius_);
    assert(bore_depth > 0.0 && bore_depth < length_);
    // The tip hemisphere has to fit inside the bore it caps.
    assert(!rounded_tip || bore_radius <= bore_depth);
    assert(bore_fits(bore_radius, bore_depth, radius_, length_, bullet_radius_)
           && "Bore is wider than the crystal where it ends");
    bore_hole_ = BoreHoleConfig{bore_radius, bore_depth, rounded_tip};
}

void Geometry::set_bullet_radius(double bullet_radius) {
    assert(shape_ == DetectorShape::Cylinder &&
           "Bulletization only for cylindrical detectors");
    assert(bullet_radius >= 0.0);
    // rho_c = radius_ - bullet_radius must stay positive: a fully rounded nose
    // (bullet_radius == radius_) is a different solid and is not supported.
    assert(bullet_radius < radius_);
    assert(bullet_radius < length_);
    // Checked in both setters so the order they are called in does not matter.
    assert((!bore_hole_ || bore_fits(bore_hole_->radius, bore_hole_->depth,
                                     radius_, length_, bullet_radius))
           && "Fillet would narrow the crystal past the existing bore");
    bullet_radius_ = bullet_radius;
}

void Geometry::set_dead_layer(double front, double side, double back) {
    dead_layer_ = DeadLayerConfig{front, side, back};
}

void Geometry::add_attenuator(const Material* material,
                              double front_thickness, double side_thickness,
                              double z_start, double z_end) {
    attenuators_.push_back({material, front_thickness, side_thickness, z_start, z_end, false});
}

void Geometry::add_collimator(const Material* material, double side_thickness,
                              double z_start, double z_end) {
    attenuators_.push_back({material, 0.0, side_thickness, z_start, z_end, true});
}

double Geometry::detector_radius() const {
    assert(shape_ == DetectorShape::Cylinder);
    return radius_;
}

double Geometry::detector_half_x() const {
    assert(shape_ == DetectorShape::Box);
    return half_x_;
}

double Geometry::detector_half_y() const {
    assert(shape_ == DetectorShape::Box);
    return half_y_;
}

double Geometry::detector_length() const {
    return length_;
}

double Geometry::outer_bounding_radius() const {
    double r = 0.0;

    if (shape_ == DetectorShape::Cylinder) {
        r = radius_;
    } else {
        r = std::sqrt(half_x_ * half_x_ + half_y_ * half_y_);
    }

    // Add dead layer side thickness
    if (dead_layer_) {
        r += dead_layer_->side;
    }

    // Add all attenuator side thicknesses
    for (const auto& att : attenuators_) {
        r += att.side_thickness;
    }

    return r;
}

std::pair<double, double> Geometry::outer_z_extent() const {
    double z_min = 0.0;  // Detector front face
    double z_max = length_;

    // Dead layer extends the z range
    if (dead_layer_) {
        z_min -= dead_layer_->front;
        z_max += dead_layer_->back;
    }

    // Attenuators may extend further
    for (const auto& att : attenuators_) {
        z_min = std::min(z_min, att.z_start - att.front_thickness);
        z_max = std::max(z_max, att.z_end);
    }

    return {z_min, z_max};
}

bool Geometry::ray_hits_outer_boundary(const Eigen::Vector3d& origin,
                                       const Eigen::Vector3d& direction) const {
    double r = outer_bounding_radius();
    auto [z_min, z_max] = outer_z_extent();

    if (shape_ == DetectorShape::Cylinder) {
        auto hit = intersect_cylinder(origin, direction, r, z_min, z_max);
        return hit.has_value() && hit->valid();
    } else {
        double hx = half_x_, hy = half_y_;
        if (dead_layer_) {
            hx += dead_layer_->side;
            hy += dead_layer_->side;
        }
        for (const auto& att : attenuators_) {
            hx += att.side_thickness;
            hy += att.side_thickness;
        }
        auto hit = intersect_box(origin, direction, hx, hy, z_min, z_max);
        return hit.has_value() && hit->valid();
    }
}


// ---- Full ray trace through nested geometry ----

std::vector<PathSegment> Geometry::trace_ray(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction) const
{
    std::vector<PathSegment> segments;
    trace_ray(origin, direction, segments);
    return segments;
}

void Geometry::trace_ray(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    std::vector<PathSegment>& segments) const
{
    segments.clear();
    const size_t target_capacity = 8 + attenuators_.size() * 4;
    if (segments.capacity() < target_capacity) {
        segments.reserve(target_capacity);
    }

    if (shape_ == DetectorShape::Cylinder) {
        trace_cylinder_geometry(origin, direction, segments);
    } else {
        trace_box_geometry(origin, direction, segments);
    }
}


void Geometry::trace_cylinder_geometry(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    std::vector<PathSegment>& segments) const
{
    // Direct shell subtraction for active crystal, dead layer, and attenuators.
    // Keep this path allocation-light; it runs in the innermost transport loops.
    double crystal_r = radius_;
    double crystal_z_min = 0.0;
    double crystal_z_max = length_;

    // Dead layer shell (if present)
    double dl_r = crystal_r;
    double dl_z_min = crystal_z_min;
    double dl_z_max = crystal_z_max;

    if (dead_layer_) {
        // The dead layer sits around the active crystal.
        // Active crystal shrinks inward by dead layer thicknesses.
        double active_r = crystal_r - dead_layer_->side;
        double active_z_min = crystal_z_min + dead_layer_->front;
        double active_z_max = crystal_z_max - dead_layer_->back;

        // Active crystal uses the reduced dimensions
        dl_r = active_r;
        dl_z_min = active_z_min;
        dl_z_max = active_z_max;
    }

    // Rounded front outer edge, if configured.  Sharp crystals -- the common
    // case -- pay only this predicate and take the original code paths below,
    // bit for bit.
    const bool bulletized = (bullet_radius_ > 0.0);
    FrontFillet outer_fillet{0.0, 0.0, 0.0, 0.0};
    FrontFillet active_fillet{0.0, 0.0, 0.0, 0.0};

    if (bulletized) {
        outer_fillet = make_front_fillet(crystal_r, crystal_z_min, bullet_radius_);

        // The active surface is the crystal surface offset inward by the dead
        // layer, so its fillet is the arc tangent to both offset faces.  For
        // equal front/side thicknesses that is exactly the fillet of radius
        // r_b - t about the same centre; for unequal ones we take the smaller
        // (more conservative) radius.  Dead layers are ~0.1 mm against an
        // ~8 mm fillet, so the convention is far below any measurable effect.
        double r_b_active = bullet_radius_;
        if (dead_layer_) {
            r_b_active -= std::max(dead_layer_->front, dead_layer_->side);
            if (r_b_active < 0.0) r_b_active = 0.0;
        }
        active_fillet = make_front_fillet(dl_r, dl_z_min, r_b_active);
        assert(dl_r > 0.0 && r_b_active < dl_r &&
               dl_z_min + r_b_active <= dl_z_max &&
               "Dead layer leaves no room for the bulletized active volume");
    }

    // Active crystal (innermost)
    if (bore_hole_) {
        double bore_z_start = crystal_z_max - bore_hole_->depth;
        RayHit bore_segments[2];
        int n_seg = (bulletized || bore_hole_->rounded_tip)
            ? intersect_shaped_bored_cylinder(
                  origin, direction,
                  dl_r, bore_hole_->radius,
                  dl_z_min, dl_z_max,
                  bore_z_start, crystal_z_max,
                  active_fillet, bore_hole_->rounded_tip,
                  bore_segments)
            : intersect_bored_cylinder(
                  origin, direction,
                  dl_r, bore_hole_->radius,
                  dl_z_min, dl_z_max,
                  bore_z_start, crystal_z_max,
                  bore_segments);

        for (int i = 0; i < n_seg; ++i) {
            if (bore_segments[i].valid()) {
                double t0 = std::max(bore_segments[i].t_enter, 0.0);
                double t1 = bore_segments[i].t_exit;
                segments.push_back({t0, t1, detector_material_, true});
            }
        }
    } else {
        auto hit = bulletized
            ? intersect_bulletized_cylinder(origin, direction,
                                            dl_r, dl_z_min, dl_z_max, active_fillet)
            : intersect_cylinder(origin, direction, dl_r, dl_z_min, dl_z_max);
        if (hit && hit->valid()) {
            double t0 = std::max(hit->t_enter, 0.0);
            segments.push_back({t0, hit->t_exit, detector_material_, true});
        }
    }

    // Dead layer
    if (dead_layer_) {
        auto outer_hit = bulletized
            ? intersect_bulletized_cylinder(origin, direction, crystal_r,
                                            crystal_z_min, crystal_z_max, outer_fillet)
            : intersect_cylinder(origin, direction,
                                 crystal_r, crystal_z_min, crystal_z_max);
        auto inner_hit = bulletized
            ? intersect_bulletized_cylinder(origin, direction, dl_r,
                                            dl_z_min, dl_z_max, active_fillet)
            : intersect_cylinder(origin, direction, dl_r, dl_z_min, dl_z_max);

        if (outer_hit && outer_hit->valid()) {
            double t0_out = std::max(outer_hit->t_enter, 0.0);
            double t1_out = outer_hit->t_exit;

            if (inner_hit && inner_hit->valid()) {
                double t0_in = std::max(inner_hit->t_enter, 0.0);
                double t1_in = inner_hit->t_exit;

                // Dead layer segments: [t0_out, t0_in] and [t1_in, t1_out]
                if (t0_out < t0_in - 1e-12) {
                    segments.push_back({t0_out, t0_in, detector_material_, false});
                }
                if (t1_in < t1_out - 1e-12) {
                    segments.push_back({t1_in, t1_out, detector_material_, false});
                }
            } else {
                // Inner crystal missed — entire outer hit is dead layer
                segments.push_back({t0_out, t1_out, detector_material_, false});
            }
        }
    }

    // Attenuator shells (from innermost to outermost)
    double inner_r = crystal_r;
    double inner_z_min = crystal_z_min;
    double inner_z_max = crystal_z_max;

    for (const auto& att : attenuators_) {
        double shell_r = inner_r + att.side_thickness;
        double front_z_min = inner_z_min - att.front_thickness;
        double shell_z_min = std::min(front_z_min, att.z_start);
        double shell_z_max = std::max(inner_z_max, att.z_end);

        auto outer_hit = intersect_cylinder(origin, direction,
                                            shell_r, shell_z_min, shell_z_max);

        // For side_only (collimator/tube) attenuators, the inner bore extends
        // the full z range of the outer shell, creating an open tube.
        // For regular (cup) attenuators, the inner subtraction matches the
        // previous shell extent.
        double sub_z_min = att.side_only ? shell_z_min : inner_z_min;
        double sub_z_max = att.side_only ? shell_z_max : inner_z_max;
        auto inner_hit = intersect_cylinder(origin, direction,
                                            inner_r, sub_z_min, sub_z_max);

        if (outer_hit && outer_hit->valid()) {
            double t0_out = std::max(outer_hit->t_enter, 0.0);
            double t1_out = outer_hit->t_exit;

            if (inner_hit && inner_hit->valid()) {
                double t0_in = std::max(inner_hit->t_enter, 0.0);
                double t1_in = inner_hit->t_exit;

                if (t0_out < t0_in - 1e-12) {
                    segments.push_back({t0_out, t0_in, att.material, false});
                }
                if (t1_in < t1_out - 1e-12) {
                    segments.push_back({t1_in, t1_out, att.material, false});
                }
            } else {
                // Inner shell missed — full outer hit is this attenuator
                segments.push_back({t0_out, t1_out, att.material, false});
            }
        }

        inner_r = shell_r;
        inner_z_min = shell_z_min;
        inner_z_max = shell_z_max;
    }

    // Sort segments by t_start
    std::sort(segments.begin(), segments.end(),
              [](const PathSegment& a, const PathSegment& b) {
                  return a.t_start < b.t_start;
              });
}


void Geometry::trace_box_geometry(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& direction,
    std::vector<PathSegment>& segments) const
{
    // Similar approach to cylinder but using box intersection.
    // No bore hole support for box detectors.

    double active_hx = half_x_;
    double active_hy = half_y_;
    double active_z_min = 0.0;
    double active_z_max = length_;

    double crystal_hx = half_x_;
    double crystal_hy = half_y_;
    double crystal_z_min = 0.0;
    double crystal_z_max = length_;

    if (dead_layer_) {
        active_hx = half_x_ - dead_layer_->side;
        active_hy = half_y_ - dead_layer_->side;
        active_z_min = dead_layer_->front;
        active_z_max = length_ - dead_layer_->back;
    }

    // Active crystal
    if (active_hx > 0.0 && active_hy > 0.0 && active_z_min < active_z_max) {
        auto hit = intersect_box(origin, direction,
                                 active_hx, active_hy, active_z_min, active_z_max);
        if (hit && hit->valid()) {
            double t0 = std::max(hit->t_enter, 0.0);
            segments.push_back({t0, hit->t_exit, detector_material_, true});
        }
    }

    // Dead layer
    if (dead_layer_) {
        auto outer_hit = intersect_box(origin, direction,
                                       crystal_hx, crystal_hy,
                                       crystal_z_min, crystal_z_max);
        auto inner_hit = intersect_box(origin, direction,
                                       active_hx, active_hy,
                                       active_z_min, active_z_max);

        if (outer_hit && outer_hit->valid()) {
            double t0_out = std::max(outer_hit->t_enter, 0.0);
            double t1_out = outer_hit->t_exit;

            if (inner_hit && inner_hit->valid()) {
                double t0_in = std::max(inner_hit->t_enter, 0.0);
                double t1_in = inner_hit->t_exit;

                if (t0_out < t0_in - 1e-12) {
                    segments.push_back({t0_out, t0_in, detector_material_, false});
                }
                if (t1_in < t1_out - 1e-12) {
                    segments.push_back({t1_in, t1_out, detector_material_, false});
                }
            } else {
                segments.push_back({t0_out, t1_out, detector_material_, false});
            }
        }
    }

    // Attenuator shells
    double inner_hx = crystal_hx;
    double inner_hy = crystal_hy;
    double inner_z_min = crystal_z_min;
    double inner_z_max = crystal_z_max;

    for (const auto& att : attenuators_) {
        double shell_hx = inner_hx + att.side_thickness;
        double shell_hy = inner_hy + att.side_thickness;
        double front_z_min = inner_z_min - att.front_thickness;
        double shell_z_min = std::min(front_z_min, att.z_start);
        double shell_z_max = std::max(inner_z_max, att.z_end);

        auto outer_hit = intersect_box(origin, direction,
                                       shell_hx, shell_hy, shell_z_min, shell_z_max);
        auto inner_hit = intersect_box(origin, direction,
                                       inner_hx, inner_hy, inner_z_min, inner_z_max);

        if (outer_hit && outer_hit->valid()) {
            double t0_out = std::max(outer_hit->t_enter, 0.0);
            double t1_out = outer_hit->t_exit;

            if (inner_hit && inner_hit->valid()) {
                double t0_in = std::max(inner_hit->t_enter, 0.0);
                double t1_in = inner_hit->t_exit;

                if (t0_out < t0_in - 1e-12) {
                    segments.push_back({t0_out, t0_in, att.material, false});
                }
                if (t1_in < t1_out - 1e-12) {
                    segments.push_back({t1_in, t1_out, att.material, false});
                }
            } else {
                segments.push_back({t0_out, t1_out, att.material, false});
            }
        }

        inner_hx = shell_hx;
        inner_hy = shell_hy;
        inner_z_min = shell_z_min;
        inner_z_max = shell_z_max;
    }

    // Sort segments by t_start
    std::sort(segments.begin(), segments.end(),
              [](const PathSegment& a, const PathSegment& b) {
                  return a.t_start < b.t_start;
              });
}


// ---- Free functions ----

double compute_transmission(const std::vector<PathSegment>& segments,
                            double energy_MeV) {
    double total_mu_L = 0.0;

    for (const auto& seg : segments) {
        if (seg.is_scoring) continue;  // Only attenuating segments
        if (!seg.material) continue;   // Skip vacuum
        double mu = seg.material->mu_total(energy_MeV);
        total_mu_L += mu * seg.length();
    }

    return std::exp(-total_mu_L);
}

double compute_active_path_length(const std::vector<PathSegment>& segments) {
    double total = 0.0;
    for (const auto& seg : segments) {
        if (seg.is_scoring) {
            total += seg.length();
        }
    }
    return total;
}

double compute_transmission_no_rayleigh(const std::vector<PathSegment>& segments,
                                         double energy_MeV) {
    double total_mu_L = 0.0;

    for (const auto& seg : segments) {
        if (seg.is_scoring) continue;
        if (!seg.material) continue;
        MacroscopicXS xs = seg.material->macroscopic_xs(energy_MeV);
        double mu_no_rs = xs.mu_pe + xs.mu_cs + xs.mu_pp;
        total_mu_L += mu_no_rs * seg.length();
    }

    return std::exp(-total_mu_L);
}

RayleighTransmissionResult compute_transmission_stochastic_rayleigh(
    const Eigen::Vector3d& start_position,
    const Eigen::Vector3d& start_direction,
    double energy_keV,
    const Geometry& geometry,
    std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::poisson_distribution<int> poisson; // parameterized per-segment

    double energy_MeV = energy_keV * 1e-3;
    Eigen::Vector3d position = start_position;
    Eigen::Vector3d direction = start_direction;
    double weight = 1.0;

    constexpr int kMaxRetraces = 10;
    std::vector<PathSegment> segments;

    for (int retrace = 0; retrace <= kMaxRetraces; ++retrace) {
        geometry.trace_ray(position, direction, segments);

        if (segments.empty()) {
            return {weight, position, direction, false};
        }

        bool did_scatter = false;

        for (const auto& seg : segments) {
            if (seg.is_scoring && seg.material) {
                // Reached scoring volume — return entry point
                Eigen::Vector3d entry = position + direction * std::max(seg.t_start, 0.0);
                return {weight, entry, direction, true};
            }

            if (!seg.material) continue;

            double L = seg.length();
            if (L <= 0.0) continue;

            MacroscopicXS xs = seg.material->macroscopic_xs(energy_MeV);
            double mu_no_rs = xs.mu_pe + xs.mu_cs + xs.mu_pp;
            double mu_rs = xs.mu_rs;

            // Analytical weight for energy-changing processes
            weight *= std::exp(-mu_no_rs * L);

            if (weight < 1e-30) {
                return {0.0, position, direction, false};
            }

            // Stochastic Rayleigh scattering
            double lambda_rs = mu_rs * L;
            int n_scatter = 0;

            if (lambda_rs > 0.01) {
                poisson = std::poisson_distribution<int>(lambda_rs);
                n_scatter = poisson(rng);
            } else if (lambda_rs > 1e-10) {
                // Bernoulli shortcut for thin layers
                n_scatter = (uniform(rng) < lambda_rs) ? 1 : 0;
            }

            if (n_scatter > 0) {
                // Sample scatter position: uniform within the segment
                double frac = uniform(rng);
                double t_scatter = std::max(seg.t_start, 0.0) + frac * L;
                position = position + direction * t_scatter;

                // Select element for Rayleigh scattering (use first element by
                // XS weight — for compounds, select proportionally)
                int Z_rs = seg.material->select_element(energy_MeV, 2 /*RS*/, uniform(rng));

                // Apply all n_scatter Rayleigh deflections
                for (int i = 0; i < n_scatter; ++i) {
                    double cos_theta = sample_rayleigh_cos_theta(Z_rs, energy_keV, rng);
                    double phi = kTwoPi * uniform(rng);
                    direction = rotate_direction(direction, cos_theta, phi);
                }

                // Must re-trace from new position/direction
                did_scatter = true;
                break;
            }
        }

        if (!did_scatter) {
            // Traversed all segments without Rayleigh scatter, no scoring hit
            return {weight, position, direction, false};
        }
    }

    // Exceeded max retraces — treat as miss
    return {weight, position, direction, false};
}

} // namespace ceelo
