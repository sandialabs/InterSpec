/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
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

/** Parity suite for the migration of the double-valued ray-tracing primitives in
 GammaInteractionCalc onto the shared templated implementations in GammaInteractionCalc_imp.hpp.

 Three oracles are in play, and the distinction matters when reading a failure:

   1. LegacyGeom::*        - a FROZEN copy of the double bodies that used to live in
                             src/GammaInteractionCalc.cpp (see LegacyGeometryRef.h).  "What shipped."
   2. GeomRef::*           - independent reference implementations using *different formulations*,
                             plus formulation-free invariant checks (see GeometryReference.h).
                             "What is true."
   3. GammaInteractionCalc - the production functions, which are now thin `_imp<double>` wrappers.

 The suite proves two things:
   (a) The migration changed nothing outside three documented divergences.  Cases 1-3 sweep legacy
       against production and require agreement.
   (b) None of those three divergences is a regression.  Two of them (Cases 4 and 5, both in
       rectangle_exit_location) are outright bug FIXES: those cases check production against GeomRef
       and *also* assert that the legacy copy still exhibits the bug, so if someone ever "fixes" the
       frozen oracle they fail loudly and its "known wrong" contract stays honest.  The third (the
       sphere, Cases 6 and 7) turned out to be neither better nor worse - just a different
       formulation, agreeing with truth to ~1e-13 of the sphere radius either way.  See the long
       note on Case 7; it corrects the premise this migration was planned under.

 Case 8 covers the point-source chord unification (the hand-rolled detector position and per-layer
 `switch` in ShieldingSourceChi2Fcn::energy_chi_contributions being replaced by
 detector_geom_from_config + center_ray_exit_distance).

 Case 9 instantiates the templates at `ceres::Jet` and requires the scalar lane to reproduce the
 double answer - the whole premise of collapsing the two copies is that both types run the same
 code, and nothing else here exercises that.  Case 10 covers the two caller contracts no sweep
 reaches: passing the same array as source and exit point, and exit_point_of_sphere_z's
 input-validation exits.

 There is deliberately no golden-value file here: every number is either an agreement between two
 independent implementations or an assertion about a geometric invariant, so nothing needs
 re-recording when a compiler changes.

 This suite retires together with LegacyGeometryRef.h once the migration has shipped a release cycle.
 */

#include "InterSpec_config.h"

#ifdef _WIN32
#include "winsock2.h"
#include "Windows.h"
#endif

#include <array>
#include <cmath>
#include <cfloat>
#include <random>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "ceres/jet.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

#include "LegacyGeometryRef.h"
#include "GeometryReference.h"

#define BOOST_TEST_MODULE GeometryLegacyParity_suite
#include <boost/test/included/unit_test.hpp>

using namespace std;
using namespace GammaInteractionCalc;


namespace
{
/** Cap on how many individual mismatches a case will report, so a systematic break produces a
 readable failure instead of 100k lines. */
const int sm_max_reported = 3;

/** The legacy and templated cylinder / rectangle-intersection bodies are line-for-line identical, so
 in practice these comparisons are bit-exact.  We allow ~4.5 ulp anyway, purely to absorb a possible
 difference in FMA contraction between the header TU and GammaInteractionCalc.cpp (clang defaults to
 -ffp-contract=on, and the two are inlined into different contexts).  A real logic divergence is many
 orders of magnitude larger than this, so the tolerance costs no diagnostic power.  If a failure ever
 lands just above it, suspect contraction before suspecting logic. */
const double sm_bitwise_tol = 1.0E-15;

struct Mismatch
{
  int count = 0;

  /** Compare two scalars, reporting at most sm_max_reported failures with the caller's context. */
  void check( const double lhs, const double rhs, const double rel_tol,
              const std::string &context )
  {
    const double scale = std::max( 1.0, std::max( std::fabs(lhs), std::fabs(rhs) ) );
    if( std::fabs(lhs - rhs) <= rel_tol*scale )
      return;

    ++count;
    if( count <= sm_max_reported )
    {
      BOOST_ERROR( context << ": " << std::setprecision(17) << lhs << " vs " << rhs
                   << " (diff " << (lhs - rhs) << ", tol " << rel_tol*scale << ")" );
    }
  }//check(...)

  /** Compare a distance plus all three exit-point components. */
  void check_ray( const double dist_a, const double pt_a[3],
                  const double dist_b, const double pt_b[3],
                  const double rel_tol, const std::string &context )
  {
    check( dist_a, dist_b, rel_tol, context + " dist" );
    check( pt_a[0], pt_b[0], rel_tol, context + " x" );
    check( pt_a[1], pt_b[1], rel_tol, context + " y" );
    check( pt_a[2], pt_b[2], rel_tol, context + " z" );
  }
};//struct Mismatch


std::string vec_str( const double v[3] )
{
  std::ostringstream strm;
  strm << std::setprecision(17) << "{" << v[0] << "," << v[1] << "," << v[2] << "}";
  return strm.str();
}

}//anonymous namespace


// ------------------------------------------------------------------------------------------------
//  Case 1 - the cylinder intersector.  Legacy and templated bodies are line-for-line identical
//  (after commits c1f9983e and 4c85390a), so any disagreement is a transcription slip.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( CylinderLegacyParityBitwise )
{
  std::mt19937 rng( 20260726 );
  std::uniform_real_distribution<double> uniform( 0.0, 1.0 );

  // Same dimension families the existing CylinderIntersectionReferenceSweep uses: a small lab
  //  source, a realistic detector-sized can, an extreme aspect ratio, and a unit cylinder.
  const std::vector<std::pair<double,double>> dims{
    { 2.0, 0.15 }, { 6.35, 106.68 }, { 225000.0, 1000.0 }, { 0.5, 0.5 }
  };

  Mismatch mm;
  int num_checked = 0, num_ref_checked = 0;

  for( const std::pair<double,double> &rl : dims )
  {
    const double radius = rl.first, half_len = rl.second;
    const double scale = std::max( radius, half_len );

    for( int iter = 0; iter < 500; ++iter )
    {
      // Half the sources inside the volume, half outside - the intersector must handle both.
      const bool src_inside = ((iter % 2) == 0);
      const double r_frac = src_inside ? (0.999*std::sqrt(uniform(rng)))
                                       : (1.5 + 3.5*uniform(rng));
      const double theta = 2.0*M_PI*uniform(rng);
      const double z_frac = src_inside ? (1.998*uniform(rng) - 0.999)
                                       : (6.0*uniform(rng) - 3.0);

      const double source[3] = { r_frac*radius*std::cos(theta),
                                 r_frac*radius*std::sin(theta),
                                 z_frac*half_len };

      // Detector: end-on, side-on, and arbitrary directions
      double detector[3];
      switch( iter % 4 )
      {
        case 0:  //end-on, on axis
          detector[0] = 0.0; detector[1] = 0.0; detector[2] = 20.0*scale;
          break;
        case 1:  //side-on, on axis
          detector[0] = 20.0*scale; detector[1] = 0.0; detector[2] = 0.0;
          break;
        case 2:  //end-on, laterally offset
          detector[0] = 0.3*scale; detector[1] = -0.4*scale; detector[2] = 20.0*scale;
          break;
        default:
        {
          const double dtheta = 2.0*M_PI*uniform(rng);
          const double dphi = M_PI*uniform(rng);
          const double ddist = (5.0 + 25.0*uniform(rng))*scale;
          detector[0] = ddist*std::sin(dphi)*std::cos(dtheta);
          detector[1] = ddist*std::sin(dphi)*std::sin(dtheta);
          detector[2] = ddist*std::cos(dphi);
          break;
        }
      }//switch( iter % 4 )

      if( (source[0] == detector[0]) && (source[1] == detector[1]) && (source[2] == detector[2]) )
        continue;

      for( const CylExitDir dir : { CylExitDir::TowardDetector, CylExitDir::AwayFromDetector } )
      {
        double legacy_pt[3], prod_pt[3];
        const double legacy_d = LegacyGeom::cylinder_line_intersection( radius, half_len, source,
                                                              detector, dir, legacy_pt );
        const double prod_d = cylinder_line_intersection( radius, half_len, source, detector,
                                                          dir, prod_pt );

        std::ostringstream ctx;
        ctx << "cyl R=" << radius << " L=" << half_len << " src=" << vec_str(source)
            << " det=" << vec_str(detector)
            << " dir=" << ((dir == CylExitDir::TowardDetector) ? "toward" : "away");
        mm.check_ray( legacy_d, legacy_pt, prod_d, prod_pt, sm_bitwise_tol, ctx.str() );
        ++num_checked;

        // Independently: the production answer must also match the slab-interval reference, for
        //  the interior-source case that reference is valid for.
        if( src_inside && (prod_d > 0.0) )
        {
          double ref_pt[3], ref_d = 0.0;
          const bool hit = GeomRef::reference_cyl_exit( radius, half_len, source, detector,
                                        (dir == CylExitDir::TowardDetector), ref_pt, ref_d );
          if( hit )
          {
            mm.check( prod_d, ref_d, 1.0E-9, ctx.str() + " vs reference dist" );
            ++num_ref_checked;
          }
        }//if( interior source that actually intersected )
      }//for( each direction )
    }//for( iterations )
  }//for( each dimension family )

  BOOST_TEST_MESSAGE( "CylinderLegacyParityBitwise: " << num_checked << " legacy/production pairs, "
                      << num_ref_checked << " reference cross-checks" );
  BOOST_CHECK_EQUAL( mm.count, 0 );
}//BOOST_AUTO_TEST_CASE( CylinderLegacyParityBitwise )


// ------------------------------------------------------------------------------------------------
//  Case 2 - rectangle_intersections.  Templated version is an algorithmically identical
//  transcription (it did NOT receive the zero-extent fix), so this is also bit-exact.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( RectIntersectionsLegacyParityBitwise )
{
  std::mt19937 rng( 20260727 );
  std::uniform_real_distribution<double> uniform( 0.0, 1.0 );

  const std::vector<std::array<double,3>> boxes{
    { 1.0, 1.0, 1.0 }, { 0.5, 20.0, 3.0 }, { 100.0, 0.1, 100.0 }
  };

  Mismatch mm;
  int num_checked = 0, num_hits = 0;

  for( const std::array<double,3> &half : boxes )
  {
    const double scale = std::max( half[0], std::max( half[1], half[2] ) );

    // A point on a sphere of radius `mult*scale`, guaranteed outside the box.
    auto outside_point = [&]( const double mult, double pt[3] ){
      const double th = 2.0*M_PI*uniform(rng);
      const double ph = std::acos( 2.0*uniform(rng) - 1.0 );
      const double d = mult*scale;
      pt[0] = d*std::sin(ph)*std::cos(th);
      pt[1] = d*std::sin(ph)*std::sin(th);
      pt[2] = d*std::cos(ph);
    };

    for( int iter = 0; iter < 1400; ++iter )
    {
      double source[3], detector[3];

      // Near-grazing rays (just off a face plane) are included on purpose - they are the hardest
      //  branch-selection cases.  Rays lying *exactly* in a face plane are not: both implementations
      //  share the `0*DBL_MAX == 0` formulation so they agree trivially, and the result is wrong in
      //  a way the ported contract assert (rightly) aborts a Debug build on.  That case is
      //  documented on its own in RectIntersectionsGrazingIsKnownWrong below.
      if( (iter % 7) == 0 )
      {
        // Near-face ray: aim down an axis, offset a hair off a face plane.
        const int axis = iter % 3;
        const int other = (axis + 1) % 3;
        source[0] = source[1] = source[2] = 0.0;
        detector[0] = detector[1] = detector[2] = 0.0;
        source[axis] = -5.0*scale;
        detector[axis] = 5.0*scale;
        source[other] = half[other]*(1.0 - 1.0E-9);
        detector[other] = half[other]*(1.0 - 1.0E-9);
      }else if( (iter % 11) == 0 )
      {
        // Near-corner ray
        source[0] = -5.0*scale; source[1] = half[1]*(1.0 - 1.0E-9); source[2] = half[2]*(1.0 - 1.0E-9);
        detector[0] = 5.0*scale; detector[1] = source[1]; detector[2] = source[2];
      }else
      {
        outside_point( 2.0 + 8.0*uniform(rng), source );
        outside_point( 2.0 + 8.0*uniform(rng), detector );
      }

      if( (source[0] == detector[0]) && (source[1] == detector[1]) && (source[2] == detector[2]) )
        continue;

      // Seed both output pairs identically: on the (documented-unreachable) trailing else branch
      //  neither implementation writes `exit_point`, and comparing two uninitialized arrays would
      //  be a spurious failure rather than a finding.
      double legacy_enter[3] = { -1.0E300, -1.0E300, -1.0E300 };
      double legacy_exit[3] = { -1.0E300, -1.0E300, -1.0E300 };
      double prod_enter[3] = { -1.0E300, -1.0E300, -1.0E300 };
      double prod_exit[3] = { -1.0E300, -1.0E300, -1.0E300 };
      const bool legacy_hit = LegacyGeom::rectangle_intersections( half[0], half[1], half[2],
                                        source, detector, legacy_enter, legacy_exit );
      const bool prod_hit = rectangle_intersections( half[0], half[1], half[2], source, detector,
                                        prod_enter, prod_exit );

      std::ostringstream ctx;
      ctx << "rect-inter half=" << vec_str(half.data()) << " src=" << vec_str(source)
          << " det=" << vec_str(detector);

      if( legacy_hit != prod_hit )
      {
        ++mm.count;
        if( mm.count <= sm_max_reported )
          BOOST_ERROR( ctx.str() << ": hit flag " << legacy_hit << " vs " << prod_hit );
        continue;
      }
      ++num_checked;

      if( !prod_hit )
        continue;
      ++num_hits;

      mm.check_ray( 0.0, legacy_enter, 0.0, prod_enter, sm_bitwise_tol, ctx.str() + " enter" );
      mm.check_ray( 0.0, legacy_exit, 0.0, prod_exit, sm_bitwise_tol, ctx.str() + " exit" );

      // Invariants, and the independent slab-traversal reference.
      BOOST_CHECK_MESSAGE( GeomRef::in_box( prod_enter, half.data(), 1.0E-9*scale ),
                           ctx.str() << ": enter point outside box " << vec_str(prod_enter) );
      BOOST_CHECK_MESSAGE( GeomRef::in_box( prod_exit, half.data(), 1.0E-9*scale ),
                           ctx.str() << ": exit point outside box " << vec_str(prod_exit) );

      double ref_enter[3], ref_exit[3];
      if( GeomRef::reference_rect_intersections( half.data(), source, detector, ref_enter, ref_exit ) )
      {
        mm.check_ray( 0.0, prod_enter, 0.0, ref_enter, 1.0E-12, ctx.str() + " enter vs reference" );
        mm.check_ray( 0.0, prod_exit, 0.0, ref_exit, 1.0E-12, ctx.str() + " exit vs reference" );
      }
    }//for( iterations )
  }//for( each box )

  BOOST_TEST_MESSAGE( "RectIntersectionsLegacyParityBitwise: " << num_checked << " pairs, "
                      << num_hits << " actual intersections" );
  BOOST_CHECK_EQUAL( mm.count, 0 );
}//BOOST_AUTO_TEST_CASE( RectIntersectionsLegacyParityBitwise )


// ------------------------------------------------------------------------------------------------
//  Case 2b - a LATENT bug both implementations share, documented so it is not rediscovered.
//
//  rectangle_intersections_imp did NOT receive the zero-extent fix that rectangle_exit_location_imp
//  did: it still forms `t = (intersect - source)*(1/norm)`.  When the ray lies exactly in a face
//  plane (`norm[i] == 0` and `|source[i]| == half[i]`) that evaluates (0)*DBL_MAX == 0, so t == 0 is
//  the winning "crossing" and the reported entry/exit points are just the source.
//
//  This is NOT reachable from production today - the sole caller (DistributedSrcCalc::eval_rect)
//  passes a real box with both endpoints well outside it - so it is out of scope for this migration.
//  The case exists so that whoever fixes rectangle_intersections_imp gets told to update it, and so
//  the sweep above is honest about why it excludes exactly-grazing rays.
//
//  Release-only: the trailing `|exit_point| <= half + 1e-9` contract assert inside
//  rectangle_intersections_imp catches this bad result and aborts a Debug build.  That assert is
//  doing its job, and is a second, independent record that the bug is real - so the check below
//  simply cannot run with asserts live.
// ------------------------------------------------------------------------------------------------
#ifdef NDEBUG
BOOST_AUTO_TEST_CASE( RectIntersectionsGrazingIsKnownWrong )
{
  const double half[3] = { 1.0, 1.0, 1.0 };
  const double source[3] = { -5.0, 1.0, 0.0 };     //y is exactly on the +y face plane
  const double detector[3] = { 5.0, 1.0, 0.0 };    //ray is purely +x, parallel to the y planes

  double enter_pt[3], exit_pt[3];
  const bool hit = rectangle_intersections( half[0], half[1], half[2], source, detector,
                                            enter_pt, exit_pt );

  double ref_enter[3], ref_exit[3];
  const bool ref_hit = GeomRef::reference_rect_intersections( half, source, detector,
                                                              ref_enter, ref_exit );

  // The truth: the ray grazes the +y face from x=-1 to x=+1.
  BOOST_REQUIRE( ref_hit );
  BOOST_CHECK_CLOSE( ref_enter[0], -1.0, 1.0E-10 );
  BOOST_CHECK_CLOSE( ref_exit[0], 1.0, 1.0E-10 );

  // What production actually does: reports a hit whose "exit" is the source point.
  BOOST_CHECK( hit );
  BOOST_CHECK_MESSAGE( exit_pt[0] == source[0],
      "rectangle_intersections no longer exhibits the grazing-ray `0*DBL_MAX == 0` bug (exit x = "
      << exit_pt[0] << ", expected the buggy " << source[0] << ").  If rectangle_intersections_imp"
      " was given the same short-circuit fix as rectangle_exit_location_imp, that is good - update"
      " this case to assert the correct answer, and drop the `grazing` skip in"
      " RectIntersectionsLegacyParityBitwise and the TODO in GammaInteractionCalc_imp.hpp." );
}//BOOST_AUTO_TEST_CASE( RectIntersectionsGrazingIsKnownWrong )
#endif //NDEBUG


// ------------------------------------------------------------------------------------------------
//  Case 3 - rectangle_exit_location over the NON-degenerate domain (strictly positive half-extents,
//  source strictly inside).  Legacy computes t as (delta)*(1/norm); the templated version as
//  (delta)/norm - one rounding apart, so not bit-exact but very close.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( RectExitLegacyParityNonDegenerate )
{
  std::mt19937 rng( 20260728 );
  std::uniform_real_distribution<double> uniform( 0.0, 1.0 );

  const std::vector<std::array<double,3>> boxes{
    { 1.0, 1.0, 1.0 }, { 0.5, 20.0, 3.0 }, { 100.0, 0.1, 100.0 }, { 2.54, 2.54, 15.24 }
  };

  Mismatch mm;
  int num_checked = 0;

  for( const std::array<double,3> &half : boxes )
  {
    const double scale = std::max( half[0], std::max( half[1], half[2] ) );

    for( int iter = 0; iter < 1000; ++iter )
    {
      const double source[3] = { half[0]*(1.98*uniform(rng) - 0.99),
                                 half[1]*(1.98*uniform(rng) - 0.99),
                                 half[2]*(1.98*uniform(rng) - 0.99) };

      double detector[3];
      if( (iter % 3) == 0 )
      {
        detector[0] = 0.0; detector[1] = 0.0; detector[2] = 25.0*scale;   //on-axis, the usual case
      }else
      {
        const double th = 2.0*M_PI*uniform(rng);
        const double ph = std::acos( 2.0*uniform(rng) - 1.0 );
        const double d = (5.0 + 25.0*uniform(rng))*scale;
        detector[0] = d*std::sin(ph)*std::cos(th);
        detector[1] = d*std::sin(ph)*std::sin(th);
        detector[2] = d*std::cos(ph);
      }

      double legacy_pt[3], prod_pt[3];
      const double legacy_d = LegacyGeom::rectangle_exit_location( half[0], half[1], half[2],
                                                    source, detector, legacy_pt );
      const double prod_d = rectangle_exit_location( half[0], half[1], half[2], source, detector,
                                                     prod_pt );

      std::ostringstream ctx;
      ctx << "rect-exit half=" << vec_str(half.data()) << " src=" << vec_str(source)
          << " det=" << vec_str(detector);

      mm.check_ray( legacy_d, legacy_pt, prod_d, prod_pt, 1.0E-14, ctx.str() );
      ++num_checked;

      // Both must agree with the independent slab traversal, and the exit point must genuinely be
      //  on the box surface, on the ray, at the returned distance.
      double ref_pt[3];
      const double ref_d = GeomRef::reference_rect_exit( half.data(), source, detector, ref_pt );
      mm.check( prod_d, ref_d, 1.0E-12, ctx.str() + " production vs reference" );
      mm.check( legacy_d, ref_d, 1.0E-12, ctx.str() + " legacy vs reference" );

      BOOST_CHECK_MESSAGE( GeomRef::on_box( prod_pt, half.data(), 1.0E-9*scale ),
                           ctx.str() << ": exit not on box surface " << vec_str(prod_pt) );
      BOOST_CHECK_MESSAGE( GeomRef::collinear( source, prod_pt, detector, 1.0E-9 ),
                           ctx.str() << ": exit not on the source->detector ray" );
      BOOST_CHECK_MESSAGE( GeomRef::dist_matches( prod_d, source, prod_pt, 1.0E-9 ),
                           ctx.str() << ": returned distance != |exit - source|" );
    }//for( iterations )
  }//for( each box )

  BOOST_TEST_MESSAGE( "RectExitLegacyParityNonDegenerate: " << num_checked << " pairs" );
  BOOST_CHECK_EQUAL( mm.count, 0 );
}//BOOST_AUTO_TEST_CASE( RectExitLegacyParityNonDegenerate )


// ------------------------------------------------------------------------------------------------
//  Case 4 - KNOWN DIVERGENCE #1: zero half-extent.  Legacy computes (0 - 0)*DBL_MAX == 0 and picks
//  the degenerate plane as the nearest exit, returning a zero chord.  Do NOT compare the two here;
//  instead check production against the reference, and pin the legacy bug so the frozen oracle
//  cannot be silently "fixed".
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( RectExitZeroExtentNewIsCorrect )
{
  const double nominal[3] = { 3.0, 4.0, 5.0 };
  const double origin[3] = { 0.0, 0.0, 0.0 };

  Mismatch mm;
  int num_checked = 0;

  // All 7 non-empty subsets of {w,h,d} zeroed
  for( int mask = 1; mask < 8; ++mask )
  {
    double half[3];
    for( int i = 0; i < 3; ++i )
      half[i] = (mask & (1 << i)) ? 0.0 : nominal[i];

    const std::vector<std::array<double,3>> detectors{
      { 0.0, 0.0, 100.0 },      //on the +z axis (the usual point-source ray)
      { 100.0, 0.0, 0.0 },      //on the +x axis (side-on cylinders map here)
      { 0.0, 100.0, 0.0 },
      { 30.0, -40.0, 100.0 }    //off-axis
    };

    for( const std::array<double,3> &det : detectors )
    {
      double prod_pt[3], ref_pt[3];
      const double prod_d = rectangle_exit_location( half[0], half[1], half[2], origin,
                                                     det.data(), prod_pt );
      const double ref_d = GeomRef::reference_rect_exit( half, origin, det.data(), ref_pt );

      std::ostringstream ctx;
      ctx << "rect-exit-degenerate half=" << vec_str(half) << " det=" << vec_str(det.data());

      mm.check( prod_d, ref_d, 1.0E-14, ctx.str() + " production vs reference" );
      ++num_checked;

      BOOST_CHECK_MESSAGE( GeomRef::in_box( prod_pt, half, 1.0E-9*nominal[2] ),
                           ctx.str() << ": exit outside the (degenerate) box " << vec_str(prod_pt) );
      BOOST_CHECK_MESSAGE( GeomRef::collinear( origin, prod_pt, det.data(), 1.0E-9 ),
                           ctx.str() << ": exit not on the source->detector ray" );
    }//for( each detector )

    // On-axis with a zero transverse extent: the chord is exactly the along-ray half-extent, and
    //  this is the case center_ray_exit_distance depends on being exact (a zero-thickness shield
    //  layer must contribute a chord equal to the remaining dimension, with derivative 1).
    if( (half[2] > 0.0) && (half[0] == 0.0) && (half[1] == 0.0) )
    {
      const double det_z[3] = { 0.0, 0.0, 100.0 };
      double pt[3];
      const double d = rectangle_exit_location( half[0], half[1], half[2], origin, det_z, pt );
      BOOST_CHECK_CLOSE( d, half[2], 1.0E-10 );
    }
  }//for( each zeroed-dimension subset )

  // --- Pin the legacy bug.  If this ever fails, someone edited the frozen oracle. ---
  {
    const double det_z[3] = { 0.0, 0.0, 100.0 };
    double legacy_pt[3], prod_pt[3];
    const double legacy_d = LegacyGeom::rectangle_exit_location( 0.0, 4.0, 5.0, origin, det_z,
                                                                 legacy_pt );
    const double prod_d = rectangle_exit_location( 0.0, 4.0, 5.0, origin, det_z, prod_pt );

    BOOST_CHECK_MESSAGE( legacy_d == 0.0,
        "LegacyGeometryRef.h has been modified: the frozen rectangle_exit_location is supposed to"
        " still exhibit the zero-transverse-dimension bug (return 0 for a zero half-width), but"
        " returned " << legacy_d );
    BOOST_CHECK_CLOSE( prod_d, 5.0, 1.0E-10 );
  }

  BOOST_TEST_MESSAGE( "RectExitZeroExtentNewIsCorrect: " << num_checked << " degenerate configs" );
  BOOST_CHECK_EQUAL( mm.count, 0 );
}//BOOST_AUTO_TEST_CASE( RectExitZeroExtentNewIsCorrect )


// ------------------------------------------------------------------------------------------------
//  Case 5 - KNOWN DIVERGENCE #2: strictly positive dimensions, but the ray is parallel to one axis'
//  planes AND the source sits exactly on that face.  Legacy again forms 0*DBL_MAX == 0 and returns a
//  zero chord.  This is reachable in production from eval_rect's outer-shell loop, where the source
//  IS the previous shell's exit point - i.e. exactly on a face.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( RectExitSourceOnParallelFace )
{
  const double half[3] = { 1.0, 1.0, 5.0 };
  const double source[3] = { 1.0, 0.0, 0.0 };      //exactly on the +x face
  const double detector[3] = { 1.0, 0.0, 20.0 };   //ray is purely +z, parallel to the x planes

  double legacy_pt[3], prod_pt[3], ref_pt[3];
  const double legacy_d = LegacyGeom::rectangle_exit_location( half[0], half[1], half[2],
                                                source, detector, legacy_pt );
  const double prod_d = rectangle_exit_location( half[0], half[1], half[2], source, detector,
                                                 prod_pt );
  const double ref_d = GeomRef::reference_rect_exit( half, source, detector, ref_pt );

  // The truth: the ray runs the full remaining half-depth.
  BOOST_CHECK_CLOSE( ref_d, 5.0, 1.0E-10 );
  BOOST_CHECK_CLOSE( prod_d, ref_d, 1.0E-10 );
  BOOST_CHECK_SMALL( std::fabs(prod_pt[2] - 5.0), 1.0E-12 );

  BOOST_CHECK_MESSAGE( legacy_d == 0.0,
      "LegacyGeometryRef.h has been modified: the frozen rectangle_exit_location is supposed to"
      " still return a zero chord for a source exactly on a face with the ray parallel to that"
      " face, but returned " << legacy_d );
}//BOOST_AUTO_TEST_CASE( RectExitSourceOnParallelFace )


// ------------------------------------------------------------------------------------------------
//  Case 6 - DIVERGENCE #3 (healthy regime): the sphere intersector.  This is a reformulation rather
//  than a fix (see Case 7), so both should be fine well away from the surface; this establishes
//  that, and sets the accuracy bar.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( SphereExitBothVsReference )
{
  std::mt19937 rng( 20260729 );
  std::uniform_real_distribution<double> uniform( 0.0, 1.0 );

  const std::vector<double> radii{ 0.5, 2.54, 100.0, 225000.0 };
  const std::vector<double> dist_mults{ 2.0, 20.0, 1000.0 };

  Mismatch mm;
  int num_checked = 0;

  for( const double S : radii )
  {
    for( const double mult : dist_mults )
    {
      const double obs_dist = mult*S;

      for( int iter = 0; iter < 400; ++iter )
      {
        // Uniform-ish in the ball, capped at 0.99 of the radius (the "healthy" regime)
        const double r = 0.99*S*std::cbrt( uniform(rng) );
        const double theta = 2.0*M_PI*uniform(rng);
        const double phi = std::acos( 2.0*uniform(rng) - 1.0 );
        const double source[3] = { r*std::sin(phi)*std::cos(theta),
                                   r*std::sin(phi)*std::sin(theta),
                                   r*std::cos(phi) };

        for( const bool positive : { true, false } )
        {
          double legacy_pt[3], prod_pt[3], ref_pt[3];
          const double legacy_d = LegacyGeom::exit_point_of_sphere_z( source, legacy_pt, S,
                                                                      obs_dist, positive );
          const double prod_d = exit_point_of_sphere_z( source, prod_pt, S, obs_dist, positive );
          const double ref_d = GeomRef::reference_sphere_exit( source, S, obs_dist, positive,
                                                               ref_pt );

          std::ostringstream ctx;
          ctx << "sphere S=" << S << " R=" << obs_dist << " |P|/S=" << (r/S)
              << " src=" << vec_str(source) << (positive ? " forward" : " backward");

          mm.check( prod_d, ref_d, 1.0E-12, ctx.str() + " production vs reference" );
          mm.check( legacy_d, ref_d, 1.0E-9, ctx.str() + " legacy vs reference" );
          ++num_checked;

          BOOST_CHECK_MESSAGE( GeomRef::on_sphere( prod_pt, S, 1.0E-12 ),
                               ctx.str() << ": exit not on sphere " << vec_str(prod_pt) );
          BOOST_CHECK_MESSAGE( GeomRef::dist_matches( prod_d, source, prod_pt, 1.0E-9 ),
                               ctx.str() << ": returned distance != |exit - source|" );
        }//for( each root )
      }//for( iterations )
    }//for( each observation distance )
  }//for( each radius )

  BOOST_TEST_MESSAGE( "SphereExitBothVsReference: " << num_checked << " configs" );
  BOOST_CHECK_EQUAL( mm.count, 0 );
}//BOOST_AUTO_TEST_CASE( SphereExitBothVsReference )


// ------------------------------------------------------------------------------------------------
//  Case 7 - the sphere intersector with the source approaching the shell surface, and with a
//  vanishing radius.
//
//  A CORRECTION TO THE MIGRATION'S ORIGINAL PREMISE, measured here rather than assumed: the
//  radius-scaled solver is NOT more accurate than the legacy closed form in this regime - but it is
//  not meaningfully worse either.  The two trade wins across the grid below (legacy is the more
//  accurate one perhaps twice as often), and their *worst* errors are within a few percent of each
//  other: both ~1e-13 of the sphere radius.  The worst case for either is a transverse near-surface
//  source, where the loser is a small multiple - not orders of magnitude - worse.  The
//  BOOST_TEST_MESSAGE at the end of this case prints the measured numbers; read them rather than
//  trusting this paragraph.
//
//  Neither formulation avoids the cancellation that sets that floor.  #exit_point_of_sphere_z_imp
//  forms `q = source/sphere_rad` first and the scaled solver then evaluates
//  `G = (q.u)^2 + 1 - |q|^2`; as |q| -> 1 the `1 - |q|^2` term cancels just as badly as the
//  `S^2 - r^2` it replaced, and some of the information was already lost in the division.  Factoring
//  it as `(S-r)*(S+r)` - what GeomRef::reference_sphere_exit does - is what actually avoids the
//  cancellation, which is why that reference is ~100x more accurate than either production form even
//  though `long double` is only binary64 here.
//
//  The migration is still right, because what the scaled form was written for is the ceres::Jet
//  DERIVATIVE staying finite as sphere_rad -> 0 (see its comment in GammaInteractionCalc_imp.hpp,
//  and test_ShieldingDimLimit), and because the absolute error either way is ~1e-12 of the sphere
//  radius - utterly negligible against attenuation physics.  So this case asserts the claim that
//  actually matters (accuracy at the scale of the sphere) and merely *reports* the comparison,
//  instead of asserting a "never worse" property that is not true.
//
//  If someone later wants the accuracy back, the fix is to have exit_point_of_sphere_z_imp hand the
//  scaled solver a pre-factored `(S-r)*(S+r)/S^2` instead of letting it recompute `1 - |q|^2`.
//  That would move the TShieldingGeomGolden values, so it is a deliberate separate change.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( SphereExitNearSurfaceAccuracy )
{
  const std::vector<double> ratios{ 1.0 - 1.0E-4, 1.0 - 1.0E-6, 1.0 - 1.0E-8,
                                    1.0 - 1.0E-10, 1.0 - 1.0E-12 };
  const std::vector<double> radii{ 0.5, 2.54, 100.0, 225000.0 };
  const std::vector<double> dist_mults{ 2.0, 20.0, 1000.0 };

  // The accuracy bar: the error in a chord through a sphere is only meaningful against the size of
  //  that sphere, not against the chord itself.  A chord of 1e-9*S is a numerically meaningless
  //  sliver, and demanding relative accuracy on it would just be measuring the reference's own noise
  //  (`long double` is binary64 on arm64, so the reference is only double-accurate anyway).
  //  1e-11*S is ~100x looser than the worst production error observed here, and many orders tighter
  //  than any genuine formulation defect would be.
  const double radius_rel_tol = 1.0E-11;

  int num_checked = 0, num_prod_better = 0, num_legacy_better = 0;
  double worst_prod_over_S = 0.0, worst_legacy_over_S = 0.0;

  for( const double S : radii )
  {
    for( const double mult : dist_mults )
    {
      const double obs_dist = mult*S;

      for( const double ratio : ratios )
      {
        // A few directions, including straight along +z (toward the detector), straight away, and
        //  transverse - the transverse one is where the chord is longest and conditioning worst.
        const std::vector<std::array<double,3>> dirs{
          { 0.0, 0.0, 1.0 }, { 0.0, 0.0, -1.0 }, { 1.0, 0.0, 0.0 },
          { 0.57735026918962576, 0.57735026918962576, 0.57735026918962576 }
        };

        for( const std::array<double,3> &dir : dirs )
        {
          const double source[3] = { ratio*S*dir[0], ratio*S*dir[1], ratio*S*dir[2] };

          for( const bool positive : { true, false } )
          {
            double legacy_pt[3], prod_pt[3], ref_pt[3];
            const double legacy_d = LegacyGeom::exit_point_of_sphere_z( source, legacy_pt, S,
                                                                        obs_dist, positive );
            const double prod_d = exit_point_of_sphere_z( source, prod_pt, S, obs_dist, positive );
            const double ref_d = GeomRef::reference_sphere_exit( source, S, obs_dist, positive,
                                                                 ref_pt );

            const double err_prod = std::fabs( prod_d - ref_d );
            const double err_legacy = std::fabs( legacy_d - ref_d );

            std::ostringstream ctx;
            ctx << "sphere-near-surface S=" << S << " R=" << obs_dist << " |P|/S=" << ratio
                << " dir=" << vec_str(dir.data()) << (positive ? " forward" : " backward");

            // The claim that matters: the chord is right to ~1e-11 of the sphere radius.
            BOOST_CHECK_MESSAGE( err_prod <= radius_rel_tol*S,
                ctx.str() << ": production absolute error " << err_prod << " exceeds "
                << (radius_rel_tol*S) << " (ref=" << ref_d << ", legacy err=" << err_legacy << ")" );

            // A formulation-free invariant, and one that holds to full double precision regardless
            //  of what the reference thinks the distance is.
            BOOST_CHECK_MESSAGE( GeomRef::on_sphere( prod_pt, S, 1.0E-11 ),
                                 ctx.str() << ": exit not on sphere " << vec_str(prod_pt) );

            worst_prod_over_S = std::max( worst_prod_over_S, err_prod/S );
            worst_legacy_over_S = std::max( worst_legacy_over_S, err_legacy/S );
            if( err_legacy > 2.0*err_prod )
              ++num_prod_better;
            else if( err_prod > 2.0*err_legacy )
              ++num_legacy_better;
            ++num_checked;
          }//for( each root )
        }//for( each direction )
      }//for( each surface ratio )
    }//for( each observation distance )
  }//for( each radius )

  // The small-radius family: q = source/S is fixed at 0.5, so only the radius shrinks.
  for( const double S : { 1.0E-3, 1.0E-6, 1.0E-9 } )
  {
    const double source[3] = { 0.0, 0.0, 0.5*S };
    double prod_pt[3], ref_pt[3];
    const double prod_d = exit_point_of_sphere_z( source, prod_pt, S, 100.0*S, true );
    const double ref_d = GeomRef::reference_sphere_exit( source, S, 100.0*S, true, ref_pt );

    BOOST_CHECK_MESSAGE( std::fabs(prod_d - ref_d) <= 1.0E-12*std::max(S, ref_d),
        "sphere small-radius S=" << S << ": " << prod_d << " vs reference " << ref_d );
    BOOST_CHECK_MESSAGE( std::isfinite(prod_d), "sphere small-radius S=" << S << " gave " << prod_d );
    ++num_checked;
  }

  BOOST_TEST_MESSAGE( "SphereExitNearSurfaceAccuracy: " << num_checked << " configs."
                      "  Worst error/radius: production " << worst_prod_over_S
                      << ", legacy " << worst_legacy_over_S
                      << ".  Production >2x more accurate in " << num_prod_better
                      << ", legacy >2x more accurate in " << num_legacy_better
                      << " (see the note above - the legacy wins are expected, and negligible)." );
}//BOOST_AUTO_TEST_CASE( SphereExitNearSurfaceAccuracy )


// ------------------------------------------------------------------------------------------------
//  Case 8 - the point-source chord unification.  ShieldingSourceChi2Fcn::energy_chi_contributions
//  used to hand-roll the detector position and switch over the geometry inline; it now calls
//  detector_geom_from_config + center_ray_exit_distance, the same helpers the Ceres path uses.
//  This case pins the sign conventions of that mapping - the easiest thing to get silently wrong.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( CenterRayChordVsHandBuiltDetector )
{
  const std::vector<GeometryType> geometries{
    GeometryType::Spherical, GeometryType::CylinderEndOn,
    GeometryType::CylinderSideOn, GeometryType::Rectangular
  };

  const std::vector<std::pair<double,double>> offsets{
    { 0.0, 0.0 }, { 3.0, 0.0 }, { 0.0, -4.0 }, { 3.0, -4.0 }
  };

  const std::vector<double> distances{ 25.0, 100.0 };

  const double det_radius = 2.54, det_setback = 0.0;

  for( const GeometryType geom : geometries )
  {
    for( const std::pair<double,double> &off : offsets )
    {
      for( const double dist : distances )
      {
        // --- The literal pre-refactor expression, copied from GammaInteractionCalc.cpp:5716-5720 ---
        const bool is_side_on_pt = (geom == GeometryType::CylinderSideOn);
        const double det_pos_pt[3] = {
          (is_side_on_pt ? dist : -off.first),
          (is_side_on_pt ? -off.first : -off.second),
          (is_side_on_pt ? -off.second : dist)
        };

        const DetectorGeomT<double> det = detector_geom_from_config<double>( geom, dist, det_radius,
                                                       det_setback, off.first, off.second );

        std::ostringstream ctx;
        ctx << "det-geom geom=" << static_cast<int>(geom) << " dist=" << dist
            << " offsets={" << off.first << "," << off.second << "}";

        BOOST_CHECK_MESSAGE( det.position[0] == det_pos_pt[0],
                             ctx.str() << ": x " << det.position[0] << " vs " << det_pos_pt[0] );
        BOOST_CHECK_MESSAGE( det.position[1] == det_pos_pt[1],
                             ctx.str() << ": y " << det.position[1] << " vs " << det_pos_pt[1] );
        BOOST_CHECK_MESSAGE( det.position[2] == det_pos_pt[2],
                             ctx.str() << ": z " << det.position[2] << " vs " << det_pos_pt[2] );

        // --- The chord itself, against the pre-refactor switch (nudge included) ---
        const std::vector<std::array<double,3>> dim_sets{
          { 5.0, 7.0, 9.0 },        //all non-degenerate
          { 0.0, 7.0, 9.0 },        //zero radius / width
          { 5.0, 0.0, 9.0 },        //zero half-length / height
          { 5.0, 7.0, 0.0 }         //zero depth
        };

        for( const std::array<double,3> &dims : dim_sets )
        {
          const bool degenerate = ((dims[0] == 0.0) || (dims[1] == 0.0) || (dims[2] == 0.0));

          const std::array<double,3> cum = dims;
          const double chord = center_ray_exit_distance( geom, cum, det );

          std::ostringstream cctx;
          cctx << ctx.str() << " dims=" << vec_str(dims.data());

          BOOST_CHECK_MESSAGE( std::isfinite(chord) && (chord >= 0.0),
                               cctx.str() << ": non-finite/negative chord " << chord );

          if( degenerate )
            continue;   //the pre-refactor switch is not a valid oracle here - that is the bug fixed

          // Replicate the old inline switch exactly (see the plan's Phase 3 notes).
          const double origin_pt[3] = { 0.0, 0.0, 0.0 };
          double old_exit_dist = 0.0, old_exit_point[3];
          switch( geom )
          {
            case GeometryType::Spherical:
              old_exit_dist = dims[0];
              break;

            case GeometryType::CylinderEndOn:
            case GeometryType::CylinderSideOn:
              old_exit_dist = LegacyGeom::cylinder_line_intersection( dims[0], dims[1], origin_pt,
                                          det_pos_pt, CylExitDir::TowardDetector, old_exit_point );
              break;

            case GeometryType::Rectangular:
            {
              const double tiny = 1.0E-6 * PhysicalUnits::mm;
              const double w = (dims[0] > 0.0) ? dims[0] : tiny;
              const double h = (dims[1] > 0.0) ? dims[1] : tiny;
              const double d = (dims[2] > 0.0) ? dims[2] : tiny;
              old_exit_dist = LegacyGeom::rectangle_exit_location( w, h, d, origin_pt, det_pos_pt,
                                                                   old_exit_point );
              break;
            }

            case GeometryType::NumGeometryType:
              break;
          }//switch( geom )

          BOOST_CHECK_MESSAGE( std::fabs(chord - old_exit_dist)
                               <= 1.0E-12*std::max(1.0, std::fabs(old_exit_dist)),
                               cctx.str() << ": chord " << chord << " vs pre-refactor switch "
                               << old_exit_dist );
        }//for( each dimension set )
      }//for( each distance )
    }//for( each offset )
  }//for( each geometry )

  // The behavioral change worth pinning: on-axis, a cylinder with one dimension pinned at zero used
  //  to get a ZERO chord from the intersector (its `radius<=0 || half_length<=0` guard), which
  //  silently dropped that shielding layer.  center_ray_exit_distance returns the along-ray
  //  dimension instead.
  {
    const DetectorGeomT<double> det_end = detector_geom_from_config<double>(
                                    GeometryType::CylinderEndOn, 100.0, det_radius, det_setback );
    const std::array<double,3> zero_rad{ 0.0, 7.0, 0.0 };
    BOOST_CHECK_CLOSE( center_ray_exit_distance( GeometryType::CylinderEndOn, zero_rad, det_end ),
                       7.0, 1.0E-10 );

    const DetectorGeomT<double> det_side = detector_geom_from_config<double>(
                                    GeometryType::CylinderSideOn, 100.0, det_radius, det_setback );
    const std::array<double,3> zero_len{ 5.0, 0.0, 0.0 };
    BOOST_CHECK_CLOSE( center_ray_exit_distance( GeometryType::CylinderSideOn, zero_len, det_side ),
                       5.0, 1.0E-10 );

    // ...and the legacy intersector really did return zero for both.
    double pt[3];
    const double origin_pt[3] = { 0.0, 0.0, 0.0 };
    const double det_z[3] = { 0.0, 0.0, 100.0 };
    BOOST_CHECK_MESSAGE( LegacyGeom::cylinder_line_intersection( 0.0, 7.0, origin_pt, det_z,
                                            CylExitDir::TowardDetector, pt ) == 0.0,
        "LegacyGeometryRef.h has been modified: the frozen cylinder_line_intersection is supposed"
        " to still return 0 for a zero-radius cylinder." );
  }
}//BOOST_AUTO_TEST_CASE( CenterRayChordVsHandBuiltDetector )


// ------------------------------------------------------------------------------------------------
//  Case 9 - the `T = ceres::Jet` instantiation.  The point of collapsing the double bodies onto the
//  templates is that both types run the *same* code, so the contract asserts and the
//  DEBUG_RAYTRACE_CALCS tracing now cover the Ceres path too.  Nothing above actually instantiates
//  a template at Jet, so this case does: the scalar lane must reproduce the double answer, and the
//  derivative lane must stay finite.  (On a Debug build this is also what runs the newly-ported
//  asserts against Jet arguments.)
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( TemplatesAgreeAtJetAndDouble )
{
  typedef ceres::Jet<double,1> Jet1;

  auto seeded = []( const double v ) -> Jet1 { Jet1 j(v); j.v[0] = 1.0; return j; };
  auto plain  = []( const double v ) -> Jet1 { return Jet1(v); };

  Mismatch mm;
  int num_checked = 0;

  // --- cylinder ---
  {
    const double radius = 6.35, half_len = 106.68;
    const double src[3] = { 2.0, -1.5, 20.0 };
    const double det[3] = { 30.0, -40.0, 500.0 };

    const Jet1 j_src[3] = { plain(src[0]), plain(src[1]), plain(src[2]) };
    const Jet1 j_det[3] = { plain(det[0]), plain(det[1]), plain(det[2]) };

    for( const CylExitDir dir : { CylExitDir::TowardDetector, CylExitDir::AwayFromDetector } )
    {
      double d_pt[3];
      const double d_dist = cylinder_line_intersection( radius, half_len, src, det, dir, d_pt );

      Jet1 j_pt[3];
      const Jet1 j_dist = cylinder_line_intersection_imp<Jet1>( seeded(radius), plain(half_len),
                                                                j_src, j_det, dir, j_pt );

      const std::string ctx = std::string("jet cyl ")
                              + ((dir == CylExitDir::TowardDetector) ? "toward" : "away");
      const double j_pt_a[3] = { j_pt[0].a, j_pt[1].a, j_pt[2].a };
      mm.check_ray( d_dist, d_pt, j_dist.a, j_pt_a, sm_bitwise_tol, ctx );
      BOOST_CHECK_MESSAGE( std::isfinite( j_dist.v[0] ),
                           ctx << ": d(dist)/d(radius) is not finite (" << j_dist.v[0] << ")" );
      ++num_checked;
    }//for( each direction )
  }

  // --- rectangle exit ---
  {
    const double half[3] = { 2.54, 2.54, 15.24 };
    const double src[3] = { 0.5, -1.0, 3.0 };
    const double det[3] = { 0.0, 0.0, 100.0 };

    const Jet1 j_src[3] = { plain(src[0]), plain(src[1]), plain(src[2]) };
    const Jet1 j_det[3] = { plain(det[0]), plain(det[1]), plain(det[2]) };

    double d_pt[3];
    const double d_dist = rectangle_exit_location( half[0], half[1], half[2], src, det, d_pt );

    Jet1 j_pt[3];
    const Jet1 j_dist = rectangle_exit_location_imp<Jet1>( plain(half[0]), plain(half[1]),
                                                seeded(half[2]), j_src, j_det, j_pt );

    const double j_pt_a[3] = { j_pt[0].a, j_pt[1].a, j_pt[2].a };
    mm.check_ray( d_dist, d_pt, j_dist.a, j_pt_a, sm_bitwise_tol, "jet rect-exit" );

    // The ray leaves through the +z face, so the chord is t = (half_depth - source_z)/norm_z and the
    //  derivative w.r.t. the half-depth is exactly 1/norm_z.  (It is 1 only for a ray exactly along
    //  z - i.e. the on-axis center ray center_ray_exit_distance uses - and this source is off-axis,
    //  which is the point: it checks the derivative actually flows through the division.)
    const double ray[3] = { det[0]-src[0], det[1]-src[1], det[2]-src[2] };
    const double ray_len = std::sqrt( ray[0]*ray[0] + ray[1]*ray[1] + ray[2]*ray[2] );
    BOOST_CHECK_CLOSE( j_dist.v[0], ray_len/ray[2], 1.0E-9 );
    ++num_checked;
  }

  // --- sphere ---
  {
    const double S = 100.0, R = 2000.0;
    const double src[3] = { 10.0, -20.0, 30.0 };
    const Jet1 j_src[3] = { plain(src[0]), plain(src[1]), plain(src[2]) };

    for( const bool positive : { true, false } )
    {
      double d_pt[3];
      const double d_dist = exit_point_of_sphere_z( src, d_pt, S, R, positive );

      Jet1 j_pt[3];
      const Jet1 j_dist = exit_point_of_sphere_z_imp<Jet1>( j_src, j_pt, seeded(S), plain(R),
                                                            positive );

      const std::string ctx = std::string("jet sphere ") + (positive ? "forward" : "backward");
      const double j_pt_a[3] = { j_pt[0].a, j_pt[1].a, j_pt[2].a };
      mm.check_ray( d_dist, d_pt, j_dist.a, j_pt_a, sm_bitwise_tol, ctx );
      BOOST_CHECK_MESSAGE( std::isfinite( j_dist.v[0] ),
                           ctx << ": d(dist)/d(radius) is not finite (" << j_dist.v[0] << ")" );
      ++num_checked;
    }//for( each root )
  }

  BOOST_TEST_MESSAGE( "TemplatesAgreeAtJetAndDouble: " << num_checked << " double/Jet pairs" );
  BOOST_CHECK_EQUAL( mm.count, 0 );
}//BOOST_AUTO_TEST_CASE( TemplatesAgreeAtJetAndDouble )


// ------------------------------------------------------------------------------------------------
//  Case 10 - the two contracts the production callers rely on but no sweep above exercises: passing
//  the same array as source and exit point (DistributedSrcCalc::eval_spherical and eval_rect both
//  do this), and exit_point_of_sphere_z's three input-validation exits.
// ------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE( AliasingAndValidationContracts )
{
  // --- exit_point == source.  eval_spherical does exit_point_of_sphere_z(src,src,...) and
  //     eval_rect does rectangle_exit_location(...,exit_point,det,exit_point). ---
  {
    const double S = 100.0, R = 2000.0;
    const double src[3] = { 10.0, -20.0, 30.0 };

    double sep_pt[3];
    const double sep_d = exit_point_of_sphere_z( src, sep_pt, S, R );

    double alias[3] = { src[0], src[1], src[2] };
    const double alias_d = exit_point_of_sphere_z( alias, alias, S, R );

    BOOST_CHECK_EQUAL( alias_d, sep_d );
    for( int i = 0; i < 3; ++i )
      BOOST_CHECK_EQUAL( alias[i], sep_pt[i] );
  }

  {
    const double half[3] = { 2.54, 2.54, 15.24 };
    const double src[3] = { 0.5, -1.0, 3.0 };
    const double det[3] = { 0.0, 0.0, 100.0 };

    double sep_pt[3];
    const double sep_d = rectangle_exit_location( half[0], half[1], half[2], src, det, sep_pt );

    double alias[3] = { src[0], src[1], src[2] };
    const double alias_d = rectangle_exit_location( half[0], half[1], half[2], alias, det, alias );

    BOOST_CHECK_EQUAL( alias_d, sep_d );
    for( int i = 0; i < 3; ++i )
      BOOST_CHECK_EQUAL( alias[i], sep_pt[i] );
  }

  // --- exit_point_of_sphere_z input validation.  The refactor moved these into the template and
  //     folded the offending values into the exception message instead of printing them to cerr. ---
  {
    const double src[3] = { 0.0, 0.0, 1.0 };
    double pt[3];

    // Detector inside the sphere.
    BOOST_CHECK_THROW( exit_point_of_sphere_z( src, pt, 10.0, 5.0 ), std::runtime_error );

    // The values reach whoever catches it (this is what replaced the cerr diagnostics).
    try
    {
      exit_point_of_sphere_z( src, pt, 10.0, 5.0 );
      BOOST_ERROR( "exit_point_of_sphere_z(obs_dist < sphere_rad) did not throw" );
    }catch( const std::runtime_error &e )
    {
      const std::string msg = e.what();
      BOOST_CHECK_MESSAGE( msg.find("obs_dist=") != std::string::npos, "message lacks obs_dist: " << msg );
      BOOST_CHECK_MESSAGE( msg.find("sphere_rad=") != std::string::npos, "message lacks sphere_rad: " << msg );
    }

    // Source outside the sphere by more than the rounding tolerance.
    const double way_out[3] = { 0.0, 0.0, 20.0 };
    BOOST_CHECK_THROW( exit_point_of_sphere_z( way_out, pt, 10.0, 100.0 ), std::runtime_error );

    // ...but a source *just* outside is treated as a rounding error: zero chord, exit == source.
    const double S = 10.0;
    const double just_out[3] = { 0.0, 0.0, S*(1.0 + 1.0E-8) };
    double near_pt[3] = { -1.0, -1.0, -1.0 };
    const double near_d = exit_point_of_sphere_z( just_out, near_pt, S, 100.0 );
    BOOST_CHECK_EQUAL( near_d, 0.0 );
    for( int i = 0; i < 3; ++i )
      BOOST_CHECK_EQUAL( near_pt[i], just_out[i] );
  }

  // --- the distance() forwarder, which nothing else covers ---
  {
    const double a[3] = { 1.0, 2.0, 3.0 };
    const double b[3] = { 4.0, 6.0, 3.0 };
    BOOST_CHECK_EQUAL( distance( a, b ), 5.0 );
    BOOST_CHECK_EQUAL( distance( a, a ), 0.0 );
  }
}//BOOST_AUTO_TEST_CASE( AliasingAndValidationContracts )
