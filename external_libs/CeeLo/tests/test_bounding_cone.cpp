/* CeeLo: aperture sampling-cone coverage tests.
 *
 * make_aperture_quadrature() cones around the ACTIVE CRYSTAL envelope rather than the outer
 * shell's circumscribed sphere. The cone is a variance-reduction device, so the ONE property
 * that must never break is coverage: every direction that reaches the active crystal has to lie
 * inside it. A cone slightly too wide only wastes rays; a cone slightly too narrow silently
 * deletes real paths and biases the efficiency low, with no error anywhere.
 *
 * These tests therefore verify coverage by brute force -- densely sampling the FULL sphere,
 * ray-tracing each direction, and asserting that every direction with a non-zero active chord is
 * inside the cone -- across the source placements and detector shapes where a cone construction
 * could plausibly go wrong.
 */

#define BOOST_TEST_MODULE BoundingConeTests
#include <boost/test/unit_test.hpp>

#include "geometry/Geometry.h"
#include "io/ResponseKernel.h"
#include "materials/Material.h"

#include <Eigen/Core>
#include <cmath>
#include <string>
#include <vector>

using namespace ceelo;

namespace {

constexpr double kPi = 3.14159265358979323846;

const Material& mat_NaI() { static Material m = make_NaI();      return m; }
const Material& mat_Ge()  { static Material m = make_HPGe();       return m; }
const Material& mat_Al()  { static Material m = make_Aluminum(); return m; }

/// Bare 3"x3" NaI cylinder.
Geometry cyl_bare() {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
    return g;
}

/// Coaxial HPGe with every feature that deforms the active volume: bulletized front edge, bore
/// hole, dead layer, and an Al can that makes the OUTER shell bigger than the crystal.
Geometry cyl_full_featured() {
    Geometry g;
    g.set_detector(DetectorShape::Cylinder, &mat_Ge(), {2.915, 6.89});
    g.set_bullet_radius(0.8);
    g.set_bore_hole(0.5, 4.0);
    g.set_dead_layer(0.0005, 0.05, 0.05);
    g.add_attenuator(&mat_Al(), 0.1, 0.1, -0.5, 7.5);
    return g;
}

/// Square box (cubic-ish) detector.
Geometry box_cubic() {
    Geometry g;
    g.set_detector(DetectorShape::Box, &mat_NaI(), {2.5, 2.5, 5.0});
    return g;
}

/// Strongly rectangular box - the case where an axisymmetric assumption would be wrong.
Geometry box_rect() {
    Geometry g;
    g.set_detector(DetectorShape::Box, &mat_NaI(), {4.0, 1.0, 3.0});
    g.add_attenuator(&mat_Al(), 0.1, 0.1, -0.3, 3.4);
    return g;
}

/// Brute-force coverage check at one source position.
///
/// Returns {n_active_found, n_outside_cone, worst_overshoot_rad}. `n_probe` directions are laid
/// over the FULL sphere by a Fibonacci spiral that is deliberately NOT the one the quadrature
/// uses (offset and a different count), so the test cannot pass by sampling the same directions
/// the sampler happened to pick.
struct Coverage { size_t active = 0, outside = 0; double worst_rad = 0.0; };

Coverage check_coverage(const Geometry& geom, const Eigen::Vector3d& src, int n_probe = 40000) {
    const ApertureQuadrature q = make_aperture_quadrature(geom, src, 256);
    Coverage cov;

    const bool full_sphere = (q.cone_omega_frac >= 1.0 - 1e-12);
    const double cos_alpha = 1.0 - 2.0 * q.cone_omega_frac;
    const double golden = kPi * (3.0 - std::sqrt(5.0));

    std::vector<PathSegment> segs;
    for (int i = 0; i < n_probe; ++i) {
        const double f = (i + 0.37) / n_probe;          // offset vs the sampler's (i+0.5)/n
        const double ct = 1.0 - 2.0 * f;
        const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
        const double ph = golden * i * 1.0001;
        const Eigen::Vector3d dir(st * std::cos(ph), st * std::sin(ph), ct);

        segs.clear();
        geom.trace_ray(src, dir, segs);
        double active = 0.0;
        for (const auto& s : segs)
            if (s.is_scoring) active += s.length();
        if (active <= 0.0) continue;

        ++cov.active;
        if (full_sphere) continue;

        const double c = q.cone_axis.dot(dir);
        if (c < cos_alpha) {
            ++cov.outside;
            cov.worst_rad = std::max(cov.worst_rad, std::acos(std::max(-1.0, std::min(1.0, c)))
                                                     - std::acos(std::max(-1.0, std::min(1.0, cos_alpha))));
        }
    }
    return cov;
}

struct Placement { const char* name; Eigen::Vector3d pos; };

/// The placements the cone has to survive. z=0 is the crystal front face, +z into the detector.
std::vector<Placement> placements(double R, double L, double hy) {
    return {
        {"near-field on-axis",        {0.0,    0.0,   -0.2}},
        {"contact on-axis",           {0.0,    0.0,   -0.02}},
        {"near-field off-axis",       {0.6*R,  0.0,   -0.5}},
        {"near-field at rim",         {R,      0.0,   -0.3}},
        {"extreme off-axis near",     {3.0*R,  0.0,   -0.2}},
        {"beside the crystal (side)", {2.0*R,  0.0,    0.5*L}},
        {"level with mid-crystal",    {1.05*R, 0.0,    0.5*L}},
        {"behind the detector",       {0.0,    0.0,    L + 1.0}},
        {"behind and off-axis",       {1.5*R,  0.0,    L + 0.5}},
        {"diagonal off-axis",         {0.8*R,  0.8*R, -0.4}},
        {"far-field on-axis",         {0.0,    0.0,   -60.0}},
        {"far-field off-axis",        {20.0,   0.0,   -60.0}},
        {"far-field extreme angle",   {60.0,   0.0,   -5.0}},
        {"far to the side",           {60.0,   0.0,    0.5*L}},

        // REGRESSION, all found by review after the first version of this test passed:
        //  - A box is NOT axisymmetric, so a source offset in BOTH x and y sits at an azimuth that
        //    is not a mirror plane. Every placement above has y == 0 or x == y, i.e. exactly the
        //    symmetric cases a box search can get right by accident. With the source beside a face
        //    at general azimuth an unguarded search settled on a cone missing ~43% of the active
        //    solid angle.
        {"beside +y face, general azi", {0.50*R, 1.05*hy, 0.50*L}},
        {"beside +y face, mid",         {0.25*R, 1.05*hy, 0.50*L}},
        {"beside +y face, near front",  {0.50*R, 1.05*hy, 0.07*L}},
        {"diagonal beside, general",    {0.60*R, 1.05*hy, 0.75*L}},
        {"behind, general azimuth",     {0.50*R, 0.70*hy, L + 0.05}},
        //  - Exactly ON a bounding plane: the true optimum is a 90-degree cone, and round-off
        //    (sin(pi) = 1.2e-16, not 0) can push the computed value a hair negative, which used to
        //    let the axis refinement run away by 30 degrees.
        {"exactly on back-face plane",  {0.0,    0.0,    L}},
        {"on back-face plane, off-axis",{0.5*R,  0.0,    L}},
        {"exactly on front-face plane", {0.0,    0.0,    0.0}},
        {"exactly on the side surface", {R,      0.0,    0.5*L}},
    };
}

/// Direct test of the bound's DEFINING claim: the cone must contain the direction to every point
/// of the crystal envelope.
///
/// This is far sharper than the ray-traced probe. Uniformly sampling 4*pi cannot resolve a leak
/// band narrower than its ~1-degree spacing, and for a distant or edge-on crystal the probe finds
/// only a handful of active directions at all -- it would never see a thin sliver escaping. Here
/// the sample points are placed ON the envelope surface, so the silhouette (where any leak must
/// occur) is sampled densely no matter how small the subtended angle is.
size_t envelope_outside_cone(const Geometry& geom, const Eigen::Vector3d& src,
                             bool cylinder, double R, double hx, double hy, double L) {
    const ApertureQuadrature q = make_aperture_quadrature(geom, src, 256);
    if (q.cone_omega_frac >= 1.0 - 1e-12) return 0;         // full sphere contains everything
    const double cos_alpha = 1.0 - 2.0 * q.cone_omega_frac;

    size_t outside = 0;
    const auto probe = [&](const Eigen::Vector3d& p) {
        const Eigen::Vector3d d = p - src;
        if (d.norm() <= 1e-12) return;
        if (q.cone_axis.dot(d.normalized()) < cos_alpha - 1e-12) ++outside;
    };

    constexpr int kN = 400;
    for (int i = 0; i <= kN; ++i) {
        const double u = double(i) / kN;
        const double z = u * L;
        for (int j = 0; j < kN; ++j) {
            const double ph = 2.0 * kPi * j / kN;
            if (cylinder) {
                probe({R * std::cos(ph), R * std::sin(ph), z});          // side wall
                probe({u * R * std::cos(ph), u * R * std::sin(ph), 0.0});// front face
                probe({u * R * std::cos(ph), u * R * std::sin(ph), L});  // back face
            } else {
                const double t = -1.0 + 2.0 * double(j) / kN;
                probe({ hx, t * hy, z}); probe({-hx, t * hy, z});        // +/-x faces
                probe({ t * hx,  hy, z}); probe({t * hx, -hy, z});       // +/-y faces
                probe({ t * hx, u * hy * (t >= 0 ? 1.0 : -1.0), 0.0});   // front face
                probe({ t * hx, u * hy * (t >= 0 ? 1.0 : -1.0), L});     // back face
            }
        }
    }
    return outside;
}

void run_suite(const Geometry& geom, const std::string& what, double R, double L) {
    const bool cylinder = (geom.shape() == DetectorShape::Cylinder);
    const double hx = cylinder ? 0.0 : geom.detector_half_x();
    const double hy = cylinder ? 0.0 : geom.detector_half_y();

    const double hy_eff = cylinder ? geom.detector_radius() : hy;
    for (const Placement& p : placements(R, L, hy_eff)) {
        const size_t env_out = envelope_outside_cone(geom, p.pos, cylinder,
                                                     cylinder ? geom.detector_radius() : 0.0,
                                                     hx, hy, geom.detector_length());
        BOOST_CHECK_MESSAGE(env_out == 0,
            what << " / " << p.name << ": " << env_out << " points of the crystal ENVELOPE lie"
            " outside the sampling cone - the cone is not a bound.");

        // Whenever a cone IS returned it must be under 90 degrees. A source outside a convex body
        // always admits one (separating hyperplane), and both extremum evaluators are valid only
        // for cos >= 0, so a non-positive cos_alpha means the construction vouched for something
        // it could not measure.
        const ApertureQuadrature qq = make_aperture_quadrature(geom, p.pos, 64);
        if (qq.cone_omega_frac < 1.0 - 1e-12) {
            BOOST_CHECK_MESSAGE(qq.cone_omega_frac < 0.5,
                what << " / " << p.name << ": returned a cone with Omega/4pi="
                << qq.cone_omega_frac << " (half-angle >= 90 deg), outside the domain where the"
                " extremum routines are exact.");
        }

        const Coverage cov = check_coverage(geom, p.pos);
        BOOST_TEST_MESSAGE("  " << what << " / " << p.name << ": active dirs=" << cov.active
                           << ", outside cone=" << cov.outside);

        BOOST_CHECK_MESSAGE(cov.outside == 0,
            what << " / " << p.name << ": " << cov.outside << " of " << cov.active
            << " directions that reach the active crystal fall OUTSIDE the sampling cone (worst by "
            << cov.worst_rad << " rad). The cone is not a bound - efficiency will be biased low.");
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE(cone_covers_bare_cylinder) {
    run_suite(cyl_bare(), "bare cylinder", 3.81, 7.62);
}

BOOST_AUTO_TEST_CASE(cone_covers_bulletized_bored_cylinder) {
    run_suite(cyl_full_featured(), "bulletized+bored HPGe", 2.915, 6.89);
}

/// Extreme aspect ratios and the geometry features the first version of this test never built.
/// These are where the axis search is hardest, so they are what a future change breaks first.
BOOST_AUTO_TEST_CASE(cone_covers_extreme_shapes) {
    Geometry needle;  needle.set_detector(DetectorShape::Cylinder, &mat_NaI(), {0.5, 20.0});
    Geometry pancake; pancake.set_detector(DetectorShape::Cylinder, &mat_NaI(), {10.0, 0.5});
    Geometry wafer;   wafer.set_detector(DetectorShape::Cylinder, &mat_NaI(), {8.0, 0.05});

    Geometry round_bore;   // rounded bore tip - never exercised anywhere before
    round_bore.set_detector(DetectorShape::Cylinder, &mat_Ge(), {2.915, 6.89});
    round_bore.set_bore_hole(0.6, 5.5, true);
    round_bore.set_bullet_radius(0.8);

    Geometry blade; blade.set_detector(DetectorShape::Box, &mat_NaI(), {6.0, 0.2, 4.0});
    Geometry bar;   bar.set_detector(DetectorShape::Box, &mat_NaI(), {0.4, 0.4, 12.0});

    run_suite(needle,     "needle cylinder",   0.5,   20.0);
    run_suite(pancake,    "pancake cylinder",  10.0,   0.5);
    run_suite(wafer,      "wafer cylinder",     8.0,   0.05);
    run_suite(round_bore, "rounded-tip bore",   2.915,  6.89);
    run_suite(blade,      "blade box",          6.0,    4.0);
    run_suite(bar,        "bar box",            0.4,   12.0);
}

/// REGRESSION: `set_dead_layer` accepts NEGATIVE thicknesses (nothing rejects them and they
/// survive an XML round-trip), and the tracer then scores `R - side` x [front, L - back], which
/// for a negative component lies OUTSIDE the crystal. The cone's entire premise is that the active
/// volume is inside the envelope it bounds, so the envelope has to be derived from those same
/// expressions rather than assuming the signs. Before that was handled, a -0.5 cm front dead layer
/// silently deleted 1.66% of the active solid angle at an ordinary far-field query.
///
/// Uses the ray-traced probe deliberately: it asks what actually SCORES, so it does not inherit
/// the envelope assumption under test.
BOOST_AUTO_TEST_CASE(negative_dead_layer_stays_inside_the_cone) {
    struct Cfg { const char* name; double front, side, back; };
    for (const Cfg& c : {Cfg{"front -0.5", -0.5, 0.0, 0.0}, Cfg{"side -0.3", 0.0, -0.3, 0.0},
                         Cfg{"back -0.4", 0.0, 0.0, -0.4}, Cfg{"all negative", -0.2, -0.2, -0.2}}) {
        Geometry g;
        g.set_detector(DetectorShape::Cylinder, &mat_NaI(), {3.81, 7.62});
        g.set_dead_layer(c.front, c.side, c.back);

        for (const Eigen::Vector3d& src : {Eigen::Vector3d{0.0, 0.0, -60.0},
                                           Eigen::Vector3d{0.0, 0.0, -1.0},
                                           Eigen::Vector3d{4.5, 0.0, 3.0}}) {
            const Coverage cov = check_coverage(g, src, 200000);
            BOOST_TEST_MESSAGE("  dead layer " << c.name << " src(" << src.x() << "," << src.y()
                               << "," << src.z() << "): active=" << cov.active
                               << " outside=" << cov.outside);
            BOOST_CHECK_MESSAGE(cov.outside == 0,
                "dead layer " << c.name << ": " << cov.outside << " of " << cov.active
                << " scoring directions fall outside the cone - the active volume escaped the"
                " envelope the cone was built from.");
        }
    }
}

BOOST_AUTO_TEST_CASE(cone_covers_cubic_box) {
    run_suite(box_cubic(), "cubic box", 2.5, 5.0);
}

BOOST_AUTO_TEST_CASE(cone_covers_rectangular_box) {
    run_suite(box_rect(), "rectangular box", 4.0, 3.0);
}

/// A source INSIDE the crystal envelope has no bounding cone; the sampler must say so rather than
/// return a cone that quietly clips the geometry.
BOOST_AUTO_TEST_CASE(inside_the_crystal_falls_back_to_full_sphere) {
    const Geometry g = cyl_bare();
    for (const Eigen::Vector3d& p : {Eigen::Vector3d(0.0, 0.0, 3.0),
                                     Eigen::Vector3d(1.0, 1.0, 1.0)}) {
        const ApertureQuadrature q = make_aperture_quadrature(g, p, 256);
        BOOST_CHECK_MESSAGE(q.cone_omega_frac >= 1.0 - 1e-12,
            "source inside the crystal should sample the full sphere, got cone_omega_frac="
            << q.cone_omega_frac);
    }
}

/// The cone must never be looser than the old bounding-sphere construction it replaced, or the
/// change is a regression dressed as an optimization.
BOOST_AUTO_TEST_CASE(never_looser_than_the_bounding_sphere) {
    const Geometry g = cyl_full_featured();
    const double r_out = g.outer_bounding_radius();
    const std::pair<double, double> z = g.outer_z_extent();
    const double half_z = 0.5 * (z.second - z.first);
    const double r_sphere = std::sqrt(r_out * r_out + half_z * half_z);
    const Eigen::Vector3d centre(0.0, 0.0, 0.5 * (z.first + z.second));

    for (const Placement& p : placements(2.915, 6.89, 2.915)) {
        const ApertureQuadrature q = make_aperture_quadrature(g, p.pos, 256);

        double sphere_frac = 1.0;
        const double dist = (centre - p.pos).norm();
        if (dist > r_sphere * 1.0001) {
            const double sin_a = r_sphere / dist;
            sphere_frac = 0.5 * (1.0 - std::sqrt(std::max(0.0, 1.0 - sin_a * sin_a)));
        }

        BOOST_CHECK_MESSAGE(q.cone_omega_frac <= sphere_frac * 1.02 + 1e-12,
            p.name << ": crystal cone Omega/4pi=" << q.cone_omega_frac
            << " is LOOSER than the old bounding-sphere cone " << sphere_frac);
    }
}

/// Coverage is the safety property; yield is the point. Record the improvement at contact.
BOOST_AUTO_TEST_CASE(contact_yield_improves) {
    const Geometry g = cyl_full_featured();
    for (const double d : {0.25, 0.5, 1.0, 2.0}) {
        const ApertureQuadrature q = make_aperture_quadrature(g, {0.0, 0.0, -d}, 4096);
        size_t active = 0;
        for (const KernelRay& r : q.rays) active += (r.active_len > 0.0f) ? 1u : 0u;
        BOOST_TEST_MESSAGE("  d=" << d << " cm: cone_omega_frac=" << q.cone_omega_frac
                           << "  active=" << active << "/" << q.n_rays_total
                           << "  (" << 100.0 * active / q.n_rays_total << "%)");
        BOOST_CHECK_MESSAGE(q.cone_omega_frac < 1.0,
            "d=" << d << " cm still falls back to the full sphere");
    }
}
