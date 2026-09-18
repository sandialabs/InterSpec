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

#define BOOST_TEST_MODULE test_DetectorGeometryCrossSection_suite
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "io/DetectorResponse.h"
#include "materials/Material.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/DetectorGeometryCrossSection.h"

using namespace std;

using Model = DetectorGeometryCrossSection::Model;
using Plane = DetectorGeometryCrossSection::Plane;
using Region = DetectorGeometryCrossSection::Region;

namespace
{
  const double kPi = 3.14159265358979323846;

  ceelo::LayerSpec layer( const int material, const double front, const double side, const double L )
  {
    ceelo::LayerSpec spec;
    spec.material_index = material;
    spec.front_thickness_cm = front;
    spec.side_thickness_cm = side;
    spec.z_start_cm = 0.0;
    spec.z_end_cm = L;
    return spec;
  }

  /** A 3x3 NaI in a 1 mm Al can. */
  ceelo::GeometryDescriptor nai_cylinder()
  {
    ceelo::GeometryDescriptor gd;
    gd.set_dimensions( ceelo::CylinderDims{ 3.81, 7.62 } );
    gd.materials = { ceelo::MaterialSpec::from( ceelo::make_NaI() ),
                     ceelo::MaterialSpec::from( ceelo::make_Aluminum() ) };
    gd.crystal_material_index = 0;
    gd.layers.push_back( layer( 1, 0.1, 0.1, 7.62 ) );
    return gd;
  }

  /** A coaxial HPGe: fillet, rounded bore, dead layer, two Al layers. */
  ceelo::GeometryDescriptor hpge_coax()
  {
    ceelo::GeometryDescriptor gd;
    gd.set_dimensions( ceelo::CylinderDims{ 3.0, 6.0 } );
    gd.materials = { ceelo::MaterialSpec::from( ceelo::make_HPGe() ),
                     ceelo::MaterialSpec::from( ceelo::make_Aluminum() ) };
    gd.crystal_material_index = 0;
    gd.bullet_radius_cm = 0.8;
    gd.bore = ceelo::BoreHoleConfig{ 0.5, 5.0, true };
    gd.dead_layer = ceelo::DeadLayerConfig{ 0.07, 0.07, 0.0 };
    gd.layers.push_back( layer( 1, 0.15, 0.15, 6.0 ) );
    gd.layers.push_back( layer( 1, 0.1, 0.1, 6.0 ) );
    return gd;
  }

  ceelo::GeometryDescriptor czt_box()
  {
    ceelo::GeometryDescriptor gd;
    gd.set_dimensions( ceelo::BoxDims{ 0.5, 0.75, 1.0 } );
    gd.symmetry = ceelo::ResponseSymmetry::Quadrant;
    gd.materials = { ceelo::MaterialSpec::from( ceelo::make_CZT() ),
                     ceelo::MaterialSpec::from( ceelo::make_Aluminum() ) };
    gd.crystal_material_index = 0;
    gd.layers.push_back( layer( 1, 0.05, 0.05, 1.0 ) );
    return gd;
  }

  ceelo::GeometryDescriptor collimated_nai()
  {
    ceelo::GeometryDescriptor gd = nai_cylinder();
    ceelo::CollimatorSpec coll;
    coll.material_index = 1;
    coll.side_thickness_cm = 0.5;
    coll.z_start_cm = -1.0;
    coll.z_end_cm = 7.62;
    gd.collimator = coll;
    return gd;
  }

  const Region *find_region( const Model &model, const string &id )
  {
    for( const Region &r : model.regions )
      if( r.id == id )
        return &r;
    return nullptr;
  }

  size_t region_index( const Model &model, const string &id )
  {
    for( size_t i = 0; i < model.regions.size(); ++i )
      if( model.regions[i].id == id )
        return i;
    return model.regions.size();
  }

  /** Every profile: non-decreasing z, 0 <= rmin <= rmax, at least two planes. */
  void check_profile( const Region &r )
  {
    BOOST_REQUIRE_MESSAGE( r.profile.size() >= 2, "Region " + r.id + " has too few planes" );
    for( size_t i = 0; i < r.profile.size(); ++i )
    {
      const Plane &p = r.profile[i];
      BOOST_CHECK_MESSAGE( p.rmin >= 0.0, r.id + ": rmin < 0" );
      BOOST_CHECK_MESSAGE( p.rmax >= p.rmin - 1e-12, r.id + ": rmax < rmin" );
      if( i > 0 )
        BOOST_CHECK_MESSAGE( p.z >= r.profile[i-1].z - 1e-12, r.id + ": z not ascending" );
    }
    BOOST_CHECK( !r.tip_title.empty() );
    BOOST_CHECK( !r.tip_lines.empty() );
  }

  /** The volume of the solid `r.profile` describes, by revolving it: each consecutive plane pair
   is a conical frustum annulus, V = pi/3 * dz * [(a0^2+a0*a1+a1^2) - (b0^2+b0*b1+b1^2)].

   This is the check that the picture and the number in its own tooltip agree - the arithmetic in
   `buildModel` is closed-form and shares nothing with the profile construction, so a disagreement
   means one of them is wrong.  Arcs are drawn as 24 chords, so a shape with a fillet or a rounded
   bore tip revolves to slightly less than the true solid.
   */
  double revolved_volume( const Region &r )
  {
    double volume = 0.0;
    for( size_t i = 1; i < r.profile.size(); ++i )
    {
      const Plane &p0 = r.profile[i-1];
      const Plane &p1 = r.profile[i];
      const double dz = p1.z - p0.z;
      if( dz <= 0.0 )
        continue;  //a step: two planes at the same z

      const double outer = p0.rmax*p0.rmax + p0.rmax*p1.rmax + p1.rmax*p1.rmax;
      const double inner = p0.rmin*p0.rmin + p0.rmin*p1.rmin + p1.rmin*p1.rmin;
      volume += (kPi/3.0) * dz * (outer - inner);
    }

    return volume;
  }//revolved_volume(...)


  /** Every region's reported volume must be the volume of the shape it draws.

   Cylinders only: a box's profile is its x-z section, and the y half-extent it is extruded through
   is not in the drawing at all, so revolving a box section says nothing about its volume.
   */
  void check_volume_matches_profile( const Model &model, const double tol_percent )
  {
    if( model.box )
      return;

    for( const Region &r : model.regions )
    {
      const double drawn = revolved_volume( r );
      BOOST_CHECK_MESSAGE( std::fabs(drawn - r.volume_cm3)
                             <= 0.01*tol_percent*std::max(drawn, r.volume_cm3),
                           "Region '" + r.id + "' reports " + std::to_string(r.volume_cm3)
                           + " cm3 but draws " + std::to_string(drawn) + " cm3" );
    }
  }//check_volume_matches_profile(...)


  bool has_knot( const vector<double> &knots, const double v )
  {
    for( const double k : knots )
      if( std::fabs(k - v) < 1e-9 )
        return true;
    return false;
  }

  void check_knots_sorted_unique( const vector<double> &knots )
  {
    for( size_t i = 1; i < knots.size(); ++i )
      BOOST_CHECK( knots[i] > knots[i-1] + 1e-12 );
  }
}//namespace


BOOST_AUTO_TEST_CASE( plain_cylinder_regions_and_volumes )
{
  const ceelo::GeometryDescriptor gd = nai_cylinder();
  const Model model = DetectorGeometryCrossSection::buildModel( gd );

  BOOST_CHECK( !model.box );
  BOOST_CHECK( !find_region( model, "dead" ) );
  BOOST_CHECK( !find_region( model, "bore" ) );
  BOOST_CHECK( !find_region( model, "collimator" ) );

  const Region *crystal = find_region( model, "crystal" );
  const Region *can = find_region( model, "layer0" );
  BOOST_REQUIRE( crystal && can );
  BOOST_CHECK( region_index( model, "layer0" ) < region_index( model, "crystal" ) );  //outermost first

  for( const Region &r : model.regions )
    check_profile( r );

  check_volume_matches_profile( model, 1e-6 );  //no arcs here: the revolve is exact

  BOOST_CHECK_CLOSE( crystal->volume_cm3, kPi*3.81*3.81*7.62, 1e-9 );
  BOOST_CHECK_CLOSE( crystal->mass_g, crystal->volume_cm3 * gd.materials[0].density_g_per_cm3, 1e-9 );
  BOOST_CHECK_EQUAL( crystal->kind, "crystal" );

  // The can is a cup: the outer cylinder less the cavity it wraps.
  const double expected_can = kPi*3.91*3.91*(7.62 + 0.1) - kPi*3.81*3.81*7.62;
  BOOST_CHECK_CLOSE( can->volume_cm3, expected_can, 1e-9 );
  BOOST_CHECK_CLOSE( can->mass_g, expected_can * gd.materials[1].density_g_per_cm3, 1e-9 );
  BOOST_CHECK_EQUAL( can->kind, "layer" );
  BOOST_CHECK( can->tip_title.find( gd.materials[1].name ) != string::npos );

  // The crystal spans the full length; the can starts 1 mm ahead of the face.
  BOOST_CHECK_CLOSE( crystal->profile.front().z, 0.0, 1e-9 );
  BOOST_CHECK_CLOSE( crystal->profile.back().z, 7.62, 1e-9 );
  BOOST_CHECK_CLOSE( can->profile.front().z, -0.1, 1e-9 );

  BOOST_CHECK_CLOSE( model.z_min, -0.1, 1e-9 );
  BOOST_CHECK_CLOSE( model.z_max, 7.62, 1e-9 );
  BOOST_CHECK_CLOSE( model.r_max, 3.91, 1e-9 );

  check_knots_sorted_unique( model.z_knots );
  check_knots_sorted_unique( model.r_knots );
  BOOST_CHECK( has_knot( model.z_knots, -0.1 ) );
  BOOST_CHECK( has_knot( model.z_knots, 0.0 ) );
  BOOST_CHECK( has_knot( model.z_knots, 7.62 ) );
  BOOST_CHECK( has_knot( model.r_knots, 0.0 ) );
  BOOST_CHECK( has_knot( model.r_knots, 3.81 ) );
  BOOST_CHECK( has_knot( model.r_knots, 3.91 ) );
}//plain_cylinder_regions_and_volumes


BOOST_AUTO_TEST_CASE( coax_volumes_add_up )
{
  const ceelo::GeometryDescriptor gd = hpge_coax();
  const Model model = DetectorGeometryCrossSection::buildModel( gd );

  const Region *crystal = find_region( model, "crystal" );
  const Region *dead = find_region( model, "dead" );
  const Region *bore = find_region( model, "bore" );
  BOOST_REQUIRE( crystal && dead && bore );
  BOOST_REQUIRE( find_region( model, "layer0" ) && find_region( model, "layer1" ) );

  for( const Region &r : model.regions )
    check_profile( r );

  // Draw order: outer layer, inner layer, dead, crystal, bore.
  BOOST_CHECK( region_index( model, "layer1" ) < region_index( model, "layer0" ) );
  BOOST_CHECK( region_index( model, "layer0" ) < region_index( model, "dead" ) );
  BOOST_CHECK( region_index( model, "dead" ) < region_index( model, "crystal" ) );
  BOOST_CHECK( region_index( model, "crystal" ) < region_index( model, "bore" ) );

  const double R = 3.0, L = 6.0, r_b = 0.8, R_bore = 0.5, D = 5.0;
  const double fillet = 2.0*kPi*(R - r_b)*r_b*r_b*(1.0 - kPi/4.0) + kPi*r_b*r_b*r_b/3.0;
  const double solid = kPi*R*R*L - fillet;
  const double bore_v = kPi*R_bore*R_bore*(D - R_bore) + (2.0/3.0)*kPi*R_bore*R_bore*R_bore;

  // The fillet and the hemispherical bore tip are drawn as 24 chords each, so the revolve sits
  //  just inside the true solid; anything worse than a fraction of a percent is a real disagreement.
  check_volume_matches_profile( model, 0.5 );

  BOOST_CHECK_CLOSE( bore->volume_cm3, bore_v, 1e-9 );
  BOOST_CHECK_EQUAL( bore->mass_g, 0.0 );
  BOOST_CHECK_EQUAL( bore->kind, "void" );
  BOOST_CHECK( dead->volume_cm3 > 0.0 );
  BOOST_CHECK( crystal->volume_cm3 < solid );
  BOOST_CHECK_CLOSE( crystal->volume_cm3 + dead->volume_cm3 + bore->volume_cm3, solid, 1e-9 );
  BOOST_CHECK_CLOSE( dead->mass_g, dead->volume_cm3 * gd.materials[0].density_g_per_cm3, 1e-9 );

  // The active volume starts behind the front dead layer; the bore runs from its apex to the back.
  BOOST_CHECK_CLOSE( crystal->profile.front().z, 0.07, 1e-9 );
  BOOST_CHECK_CLOSE( crystal->profile.back().z, L, 1e-9 );
  BOOST_CHECK_CLOSE( bore->profile.front().z, L - D, 1e-9 );
  BOOST_CHECK_CLOSE( bore->profile.back().z, L, 1e-9 );
  BOOST_CHECK_CLOSE( bore->profile.back().rmax, R_bore, 1e-9 );
  BOOST_CHECK_CLOSE( dead->profile.front().z, 0.0, 1e-9 );
  BOOST_CHECK_CLOSE( dead->profile.back().z, L, 1e-9 );

  // The fillet narrows the crystal at the very front, to R - r_b.
  BOOST_CHECK_CLOSE( dead->profile.front().rmax, R - r_b, 1e-9 );
  BOOST_CHECK( crystal->profile.front().rmax < (R - 0.07) );

  // Layers stack outward from the crystal radius.
  const Region *inner = find_region( model, "layer0" );
  const Region *outer = find_region( model, "layer1" );
  BOOST_CHECK_CLOSE( inner->profile.back().rmin, R, 1e-9 );
  BOOST_CHECK_CLOSE( inner->profile.back().rmax, R + 0.15, 1e-9 );
  BOOST_CHECK_CLOSE( outer->profile.back().rmin, R + 0.15, 1e-9 );
  BOOST_CHECK_CLOSE( outer->profile.back().rmax, R + 0.25, 1e-9 );
  BOOST_CHECK_CLOSE( outer->profile.front().z, -0.25, 1e-9 );
  BOOST_CHECK_CLOSE( model.r_max, R + 0.25, 1e-9 );

  check_knots_sorted_unique( model.z_knots );
  check_knots_sorted_unique( model.r_knots );
  BOOST_CHECK( has_knot( model.z_knots, 0.07 ) );
  BOOST_CHECK( has_knot( model.z_knots, r_b ) );
  BOOST_CHECK( has_knot( model.z_knots, L - D ) );
  BOOST_CHECK( has_knot( model.r_knots, R_bore ) );
  BOOST_CHECK( has_knot( model.r_knots, R - 0.07 ) );
}//coax_volumes_add_up


BOOST_AUTO_TEST_CASE( box_crystal )
{
  const ceelo::GeometryDescriptor gd = czt_box();
  const Model model = DetectorGeometryCrossSection::buildModel( gd );

  BOOST_CHECK( model.box );
  const Region *crystal = find_region( model, "crystal" );
  const Region *can = find_region( model, "layer0" );
  BOOST_REQUIRE( crystal && can );
  BOOST_CHECK( !find_region( model, "dead" ) );
  BOOST_CHECK( !find_region( model, "bore" ) );

  for( const Region &r : model.regions )
    check_profile( r );

  // (no revolve check: a box is extruded through y, not revolved - see check_volume_matches_profile)
  BOOST_CHECK_CLOSE( crystal->volume_cm3, 4.0*0.5*0.75*1.0, 1e-9 );

  // What the drawing *can* be held to for a box: the section is the x-z one, at the half-extents.
  BOOST_CHECK_CLOSE( crystal->profile.front().z, 0.0, 1e-9 );
  BOOST_CHECK_CLOSE( crystal->profile.back().z, 1.0, 1e-9 );
  BOOST_CHECK_CLOSE( can->profile.back().rmax, 0.55, 1e-9 );
  const double expected_can = 4.0*0.55*0.80*1.05 - 4.0*0.5*0.75*1.0;
  BOOST_CHECK_CLOSE( can->volume_cm3, expected_can, 1e-9 );
  BOOST_CHECK_CLOSE( crystal->profile.back().rmax, 0.5, 1e-9 );  //the x half-extent stands in for r
}//box_crystal


/** A flat-bottomed bore and a back dead layer - the branches the coax fixture does not reach, and
 where the share of the bore that falls inside the active crystal is not simply its depth.
 */
BOOST_AUTO_TEST_CASE( flat_bore_and_back_dead_layer )
{
  ceelo::GeometryDescriptor gd = hpge_coax();
  gd.bore = ceelo::BoreHoleConfig{ 0.5, 5.0, false };          //flat bottom
  gd.dead_layer = ceelo::DeadLayerConfig{ 0.07, 0.07, 0.3 };   //and a back dead layer

  const Model model = DetectorGeometryCrossSection::buildModel( gd );

  const Region *crystal = find_region( model, "crystal" );
  const Region *dead = find_region( model, "dead" );
  const Region *bore = find_region( model, "bore" );
  BOOST_REQUIRE( crystal && dead && bore );

  for( const Region &r : model.regions )
    check_profile( r );
  check_volume_matches_profile( model, 0.5 );   //the fillet is still an arc

  const double R = 3.0, L = 6.0, r_b = 0.8, R_bore = 0.5, D = 5.0;
  const double fillet = 2.0*kPi*(R - r_b)*r_b*r_b*(1.0 - kPi/4.0) + kPi*r_b*r_b*r_b/3.0;
  const double solid = kPi*R*R*L - fillet;

  BOOST_CHECK_CLOSE( bore->volume_cm3, kPi*R_bore*R_bore*D, 1e-9 );   //flat: a plain cylinder
  BOOST_CHECK_CLOSE( crystal->volume_cm3 + dead->volume_cm3 + bore->volume_cm3, solid, 1e-9 );

  // The active crystal stops short of the back dead layer, and the bore runs on through it.
  BOOST_CHECK_CLOSE( crystal->profile.back().z, L - 0.3, 1e-9 );
  BOOST_CHECK_CLOSE( bore->profile.back().z, L, 1e-9 );
  BOOST_CHECK( has_knot( model.z_knots, L - 0.3 ) );

  // Only the part of the bore inside the active slab comes out of the active volume.
  const double active_len = (L - 0.3) - 0.07;
  const double r_ba = std::max( 0.0, r_b - 0.07 );
  const double fillet_a = 2.0*kPi*((R - 0.07) - r_ba)*r_ba*r_ba*(1.0 - kPi/4.0) + kPi*r_ba*r_ba*r_ba/3.0;
  const double bore_in_active = kPi*R_bore*R_bore*((L - 0.3) - (L - D));
  BOOST_CHECK_CLOSE( crystal->volume_cm3,
                     kPi*(R - 0.07)*(R - 0.07)*active_len - fillet_a - bore_in_active, 1e-9 );
}//flat_bore_and_back_dead_layer


/** A material name with characters that would break a JavaScript string literal must survive into
 the JSON as data - the tooltip text is user-typed and reaches the browser through doJavaScript.
 */
/** A layer whose material is a vacuum is empty space, so it draws unfilled (kind "void", which the
    CSS gives `fill: none`) and its tooltip quotes no material or mass.  It is still a real gap, so
    its edges must stay in the knot lists - that is what keeps a thin one at the min-pixel floor.
 */
BOOST_AUTO_TEST_CASE( vacuum_layers_are_voids )
{
  auto has_knot = []( const std::vector<double> &knots, const double v ) -> bool {
    for( const double k : knots )
      if( fabs(k - v) < 1.0E-7 )
        return true;
    return false;
  };

  auto joined_lines = []( const Region &r ) -> string {
    string all = r.tip_title;
    for( const string &line : r.tip_lines )
      all += "\n" + line;
    return all;
  };

  // The user's case: a 5 mm vacuum gap between two Al layers of a coax HPGe.  Also check the two
  //  ways a gap can arrive - named as a vacuum, or carrying a density no material has.
  ceelo::GeometryDescriptor gd = hpge_coax();
  gd.materials.push_back( CeeLoUtils::vacuumMaterialSpec() );          //index 2: named "vacuum"
  ceelo::MaterialSpec unnamed_gap = CeeLoUtils::vacuumMaterialSpec();
  unnamed_gap.name.clear();                                            //index 3: blank on the form
  gd.materials.push_back( unnamed_gap );
  gd.layers.clear();
  gd.layers.push_back( layer( 1, 0.1, 0.1, 6.0 ) );    //Al
  gd.layers.push_back( layer( 2, 0.5, 0.5, 6.0 ) );    //vacuum
  gd.layers.push_back( layer( 3, 0.2, 0.2, 6.0 ) );    //blank -> vacuum
  gd.layers.push_back( layer( 1, 0.1, 0.1, 6.0 ) );    //Al

  const Model model = DetectorGeometryCrossSection::buildModel( gd );
  check_volume_matches_profile( model, 0.5 );   //a coax: the fillet and bore tip are arcs
  for( const Region &r : model.regions )
    check_profile( r );

  const Region *al_inner = find_region( model, "layer0" );
  const Region *vac = find_region( model, "layer1" );
  const Region *blank = find_region( model, "layer2" );
  const Region *al_outer = find_region( model, "layer3" );
  BOOST_REQUIRE( al_inner && vac && blank && al_outer );

  BOOST_CHECK_EQUAL( al_inner->kind, "layer" );
  BOOST_CHECK_EQUAL( al_outer->kind, "layer" );
  BOOST_CHECK_EQUAL( vac->kind, "void" );
  BOOST_CHECK_EQUAL( blank->kind, "void" );

  // No material line ("g/cm3") and no mass - it says it is empty instead.
  for( const Region *r : { vac, blank } )
  {
    const string txt = joined_lines( *r );
    BOOST_CHECK( txt.find( "Mass" ) == string::npos );
    BOOST_CHECK( txt.find( "g/cm" ) == string::npos );
    BOOST_CHECK( txt.find( "(empty)" ) != string::npos );
    BOOST_CHECK( txt.find( "Volume" ) != string::npos );
  }
  // A blank material still gets a name in the tooltip rather than "Layer 3: ".
  BOOST_CHECK( blank->tip_title.find( "vacuum" ) != string::npos );
  // A real layer keeps both.
  BOOST_CHECK( joined_lines( *al_inner ).find( "Mass" ) != string::npos );
  BOOST_CHECK( joined_lines( *al_inner ).find( "g/cm" ) != string::npos );

  // The gap's edges are knots, so the min-pixel map still gives it its own width.  Radii: crystal
  //  3.0 -> Al 3.1 -> vacuum 3.6 -> blank 3.8 -> Al 3.9; fronts step by the same thicknesses.
  BOOST_CHECK( has_knot( model.r_knots, 3.1 ) );
  BOOST_CHECK( has_knot( model.r_knots, 3.6 ) );
  BOOST_CHECK( has_knot( model.r_knots, 3.8 ) );
  BOOST_CHECK( has_knot( model.r_knots, 3.9 ) );
  BOOST_CHECK( has_knot( model.z_knots, -0.1 ) );
  BOOST_CHECK( has_knot( model.z_knots, -0.6 ) );
  BOOST_CHECK( has_knot( model.z_knots, -0.8 ) );
  BOOST_CHECK( has_knot( model.z_knots, -0.9 ) );

  // The bore stays a void, and the JSON carries the kinds the CSS keys on.
  const Region *bore = find_region( model, "bore" );
  BOOST_REQUIRE( bore );
  BOOST_CHECK_EQUAL( bore->kind, "void" );

  const nlohmann::json j = nlohmann::json::parse( DetectorGeometryCrossSection::toJson( model ) );
  int n_void = 0;
  for( const nlohmann::json &r : j["regions"] )
    n_void += (r["kind"].get<string>() == "void");
  BOOST_CHECK_EQUAL( n_void, 3 );   //two gaps and the bore
}//vacuum_layers_are_voids


BOOST_AUTO_TEST_CASE( json_escapes_material_names )
{
  ceelo::GeometryDescriptor gd = nai_cylinder();
  const string nasty = "Al\"; alert(1); //\n</script>\\";
  gd.materials[1].name = nasty;

  const Model model = DetectorGeometryCrossSection::buildModel( gd );
  const string json_txt = DetectorGeometryCrossSection::toJson( model );

  // It parses back (an unescaped quote or newline would have ended the string early), and the name
  //  comes through as data rather than as syntax.
  const nlohmann::json j = nlohmann::json::parse( json_txt );
  BOOST_CHECK( json_txt.find( "\n" ) == string::npos );  //a literal newline would break the JS string

  bool found = false;
  for( const nlohmann::json &r : j["regions"] )
  {
    if( r["tip"]["title"].get<string>().find( nasty ) != string::npos )
      found = true;
  }
  BOOST_CHECK( found );
}//json_escapes_material_names


BOOST_AUTO_TEST_CASE( collimator_and_json )
{
  const ceelo::GeometryDescriptor gd = collimated_nai();
  const Model model = DetectorGeometryCrossSection::buildModel( gd );

  const Region *coll = find_region( model, "collimator" );
  BOOST_REQUIRE( coll );
  check_volume_matches_profile( model, 1e-6 );
  BOOST_CHECK_EQUAL( region_index( model, "collimator" ), size_t(0) );  //outermost
  for( const Region &r : model.regions )
    check_profile( r );

  // Outside the can (r 3.91 -> 4.41), from its own extension past the face to the back.
  const double r_in = 3.91, r_out = 4.41, z_cs = -1.0;
  BOOST_CHECK_CLOSE( coll->profile.front().z, z_cs, 1e-9 );
  BOOST_CHECK_CLOSE( coll->profile.front().rmin, r_in, 1e-9 );
  BOOST_CHECK_CLOSE( coll->profile.front().rmax, r_out, 1e-9 );
  BOOST_CHECK_CLOSE( coll->volume_cm3, kPi*(r_out*r_out - r_in*r_in)*(7.62 - z_cs), 1e-9 );
  BOOST_CHECK_CLOSE( model.z_min, z_cs, 1e-9 );
  BOOST_CHECK_CLOSE( model.r_max, r_out, 1e-9 );

  const string json_txt = DetectorGeometryCrossSection::toJson( model );
  const nlohmann::json j = nlohmann::json::parse( json_txt );
  BOOST_CHECK( j.contains("box") && j.contains("zMin") && j.contains("zMax") && j.contains("rMax") );
  BOOST_REQUIRE( j.contains("knots") && j["knots"].contains("z") && j["knots"].contains("r") );
  BOOST_CHECK_EQUAL( j["knots"]["z"].size(), model.z_knots.size() );
  BOOST_REQUIRE( j.contains("labels") && j["labels"].contains("front") );
  BOOST_CHECK( !j["labels"]["front"].get<string>().empty() );
  BOOST_REQUIRE( j.contains("regions") && j["regions"].is_array() );
  BOOST_CHECK_EQUAL( j["regions"].size(), model.regions.size() );
  for( const nlohmann::json &r : j["regions"] )
  {
    BOOST_CHECK( r.contains("id") && r.contains("kind") && r.contains("profile") && r.contains("tip") );
    BOOST_CHECK( r["profile"].is_array() && (r["profile"].size() >= 2) );
    BOOST_CHECK_EQUAL( r["profile"][0].size(), size_t(3) );
    BOOST_CHECK( r["tip"].contains("title") && r["tip"].contains("lines") );
    BOOST_CHECK( !r["tip"]["title"].get<string>().empty() );
  }
}//collimator_and_json
