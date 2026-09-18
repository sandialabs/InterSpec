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

#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <functional>

#include <Wt/WString.h>
#include <Wt/WApplication.h>

#include <nlohmann/json.hpp>

#include "SpecUtils/StringAlgo.h"

#include "io/DetectorResponse.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorGeometryDiagram.h"

using namespace Wt;
using namespace std;

namespace
{
  const double kPi = 3.14159265358979323846;

  /** The English text of every `dgd-*` key, so #buildModel can be exercised (unit tests) with no
   Wt session; with a session the message bundle wins. */
  const char *english_text( const char *key )
  {
    static const std::pair<const char *,const char *> table[] = {
      { "dgd-front",        "front (source side)" },
      { "dgd-crystal",      "Active crystal volume" },
      { "dgd-dead",         "Dead layer (inactive crystal)" },
      { "dgd-bore",         "Bore hole (void)" },
      { "dgd-layer",        "Layer {1}: {2}" },
      { "dgd-collimator",   "Collimator: {1}" },
      { "dgd-material",     "{1} ({2} g/cm\xc2\xb3)" },
      { "dgd-dims-cyl",     "Diameter {1}, length {2}" },
      { "dgd-dims-box",     "{1} \xc3\x97 {2} \xc3\x97 {3} (width \xc3\x97 height \xc3\x97 length)" },
      { "dgd-fillet",       "Front edge fillet radius {1}" },
      { "dgd-bore-dims",    "Diameter {1}, depth {2} from the back{3}" },
      { "dgd-bore-rounded", " (rounded tip)" },
      { "dgd-dead-dims",    "Front {1}, side {2}" },
      { "dgd-dead-dims-back", "Front {1}, side {2}, back {3}" },
      { "dgd-layer-dims",   "Front {1}, side {2}; outer diameter {3}" },
      { "dgd-coll-dims",    "Thickness {1}, extends {2} past the face; inner diameter {3}" },
      { "dgd-volume",       "Volume: {1} cm\xc2\xb3" },
      { "dgd-mass",         "Mass: {1}" },
      { "dgd-void",         "(empty)" },
      { "dgd-vacuum",       "vacuum" }
    };

    for( const auto &entry : table )
    {
      if( strcmp( entry.first, key ) == 0 )
        return entry.second;
    }

    return key;
  }//english_text(...)


  /** `WString::tr(key)` in a session, else the English text. */
  WString text( const char *key )
  {
    if( WApplication::instance() )
      return WString::tr( key );
    return WString::fromUTF8( english_text(key) );
  }//text(...)


  WString length_str( const double cm )
  {
    return WString::fromUTF8( PhysicalUnits::printToBestLengthUnitsCompact( cm * PhysicalUnits::cm, 4 ) );
  }

  WString mass_str( const double grams )
  {
    return WString::fromUTF8( PhysicalUnits::printToBestMassUnits( grams * PhysicalUnits::gram ) );
  }

  WString compact_str( const double value )
  {
    return WString::fromUTF8( SpecUtils::printCompact( value, 3 ) );
  }


  /** Volume a quarter-torus fillet of radius `r_b` removes from the front edge of a cylinder of
   radius `R` (Pappus: the (1 - pi/4) r_b^2 corner cross-section swept around its centroid, plus
   the correction for the centroid sitting inside the sweep radius). */
  double fillet_volume( const double R, const double r_b )
  {
    if( (r_b <= 0.0) || (R <= r_b) )
      return 0.0;
    return 2.0*kPi*(R - r_b)*r_b*r_b*(1.0 - kPi/4.0) + kPi*r_b*r_b*r_b/3.0;
  }

  /** Volume of the bore solid between the axial positions `z0` and `z1`, given a bore of radius
   `R_b` whose apex (tip) is at `z_apex` and which runs back to `z_back`.

   A z-range rather than a depth, because the active crystal and the dead layer each take only
   part of the bore: with a hemispherical tip the bore's cross-section varies over the first
   `R_b` of its length, so "the length of the overlap" is not enough to say how much of the bore
   lies in a given slab.
   */
  double bore_volume_between( const double R_b, const double z_apex, const double z_back,
                              const bool rounded, const double z0, const double z1 )
  {
    const double lo = std::max( z0, z_apex );
    const double hi = std::min( z1, z_back );
    if( (R_b <= 0.0) || (hi <= lo) )
      return 0.0;

    // Hemispherical tip: r(z)^2 = R_b^2 - (z_apex + R_b - z)^2 over [z_apex, z_apex + R_b].
    double volume = 0.0;
    double straight_from = lo;
    if( rounded )
    {
      const double cap_end = std::min( hi, z_apex + R_b );
      if( cap_end > lo )
      {
        // Integral of pi*(R_b^2 - u^2) du, u = (z_apex + R_b) - z, taken between the two ends.
        const double u_lo = (z_apex + R_b) - cap_end;   //u decreases as z grows
        const double u_hi = (z_apex + R_b) - lo;
        auto anti = [R_b]( const double u ){ return R_b*R_b*u - u*u*u/3.0; };
        volume += kPi * ( anti(u_hi) - anti(u_lo) );
        straight_from = cap_end;
      }
    }//if( rounded )

    if( hi > straight_from )
      volume += kPi*R_b*R_b*(hi - straight_from);

    return volume;
  }


  /** Appends the sorted, unique members of `values` that lie within [lo, hi]. */
  void add_knots( std::vector<double> &knots, const std::vector<double> &values,
                  const double lo, const double hi )
  {
    for( const double v : values )
    {
      if( (v >= lo - 1e-12) && (v <= hi + 1e-12) )
        knots.push_back( v );
    }
  }

  void finish_knots( std::vector<double> &knots )
  {
    std::sort( begin(knots), end(knots) );
    knots.erase( std::unique( begin(knots), end(knots),
                              []( double a, double b ){ return std::fabs(a - b) < 1e-9; } ),
                 end(knots) );
  }
}//namespace


DetectorGeometryDiagram::DetectorGeometryDiagram()
  : WContainerWidget(),
    m_geometryJson(),
    m_stale( false ),
    m_jsDefined( false )
{
  addStyleClass( "DetectorGeometryDiagram" );

  InterSpec *viewer = InterSpec::instance();
  if( viewer )
    viewer->useMessageResourceBundle( "DetectorGeometryDiagram" );

  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/DetectorGeometryDiagram.js" );
  wApp->useStyleSheet( "InterSpec_resources/DetectorGeometryDiagram.css" );
}//DetectorGeometryDiagram constructor


DetectorGeometryDiagram::~DetectorGeometryDiagram()
{
}


void DetectorGeometryDiagram::doXsJs( const std::string &method_call )
{
  // See the note on this function in the header: the element and the object may both be absent.
  doJavaScript( "{const c=" + jsRef() + ";if(c&&c.xs){c.xs." + method_call + ";}}" );
}//doXsJs(...)


void DetectorGeometryDiagram::defineJavaScript()
{
  m_jsDefined = true;

  // minPx: the thinnest a layer may draw, so a 10-um window is still a visible line.
  // minRoomPx / minRoomHeightPx: below this much room in the row it shares, the drawing hides
  //   itself rather than being squeezed, pushed underneath, or drawn as a sliver.
  //
  // Guarded: `jsRef()` is a getElementById, which is null while this widget is only a stub on a tab
  //  that has not been shown.  Wt re-applies JavaScript members when the real element replaces the
  //  stub, so the object gets built then instead.
  setJavaScriptMember( "xs", "(function(){const e=" + jsRef()
                             + ";return e ? new DetectorGeometryDiagram(e,"
                               " {minPx: 2, minRoomPx: 200, minRoomHeightPx: 120}) : null;})()" );

  setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + id() + "') ){"
          "const c=" + jsRef() + ";"
          "if(c && c.xs)"
            "c.xs.handleResize();"
        "}"
      "}"
    "});"
  );

  doJavaScript( "{const e=" + jsRef() + ";if(e&&e.resizeObserver)e.resizeObserver.observe(e);}" );

  // Re-send what the drawing is currently showing, rather than replaying a one-shot queue: this runs
  //  again whenever the client-side object is (re)built, and the object that is built is empty.
  if( !m_geometryJson.empty() )
    doXsJs( "setData(" + m_geometryJson + ")" );
}//defineJavaScript()


void DetectorGeometryDiagram::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = flags.test( Wt::RenderFlag::Full );

  WContainerWidget::render( flags );

  // `!m_jsDefined` as well as the Full flag: a first render does not always carry it, and without
  //  the client-side object every call hits `undefined`.
  if( renderFull || !m_jsDefined )
    defineJavaScript();
}//render(...)


void DetectorGeometryDiagram::setGeometry( const ceelo::GeometryDescriptor &gd )
{
  m_geometryJson = toJson( buildModel(gd) );   //kept so defineJavaScript can re-send it
  doXsJs( "setData(" + m_geometryJson + ")" );
  setStale( false );
}//setGeometry(...)


void DetectorGeometryDiagram::setStale( const bool stale )
{
  if( stale == m_stale )
    return;
  m_stale = stale;
  toggleStyleClass( "DgdStale", stale );
}//setStale(...)


DetectorGeometryDiagram::Model DetectorGeometryDiagram::buildModel( const ceelo::GeometryDescriptor &gd )
{
  Model model;

  const bool box = (gd.shape == ceelo::DetectorShape::Box);
  const vector<double> &dims = gd.dimensions_cm;
  if( dims.size() < (box ? 3u : 2u) )
    return model;   //nothing sensible to draw

  model.box = box;

  // Transverse extents are stored as halves, the length in full (CeeLo's convention); the x
  //  half-extent plays the role of the radius for a box.
  const double R = dims[0];
  const double hy = box ? dims[1] : 0.0;
  const double L = box ? dims[2] : dims[1];
  if( (R <= 0.0) || (L <= 0.0) )
    return model;

  auto material_name = [&gd]( const int index ) -> string {
    if( (index >= 0) && (index < static_cast<int>(gd.materials.size())) )
      return gd.materials[static_cast<size_t>(index)].name;
    return string();
  };
  auto material_density = [&gd]( const int index ) -> double {
    if( (index >= 0) && (index < static_cast<int>(gd.materials.size())) )
      return gd.materials[static_cast<size_t>(index)].density_g_per_cm3;
    return 0.0;
  };

  // An empty gap rather than a material: named as one (the geometry form lets a layer material be
  //  left blank), or carrying a density no real material has (an imported geometry's filler).
  auto void_material = [&material_name, &material_density]( const int index ) -> bool {
    return CeeLoUtils::isVacuumMaterialName( material_name(index) )
           || (material_density(index) < 1.0E-6);
  };

  // --- the crystal: fillet, dead layer and bore (cylinders only) --------------------------------
  const double r_b = box ? 0.0 : std::max( 0.0, gd.bullet_radius_cm );
  double t_f = 0.0, t_s = 0.0, t_b = 0.0;
  if( !box && gd.dead_layer )
  {
    t_f = std::max( 0.0, gd.dead_layer->front );
    t_s = std::max( 0.0, gd.dead_layer->side );
    t_b = std::max( 0.0, gd.dead_layer->back );
  }
  double R_bore = 0.0, D_bore = 0.0;
  bool rounded = false;
  if( !box && gd.bore && (gd.bore->radius > 0.0) && (gd.bore->depth > 0.0) )
  {
    R_bore = gd.bore->radius;
    D_bore = std::min( gd.bore->depth, L );
    rounded = gd.bore->rounded_tip;
  }

  const double rho_c = R - r_b;                           //fillet arc centre radius
  const double R_a = std::max( 0.0, R - t_s );            //active radius
  const double z_a0 = std::min( t_f, L );                 //active volume front
  const double z_a1 = std::max( z_a0, L - t_b );          //active volume back
  const double r_ba = std::max( 0.0, r_b - std::max(t_f, t_s) );  //fillet of the active volume
  const double rho_ca = R_a - r_ba;
  const double z_apex = L - D_bore;                       //bore apex (tip), from the front

  // Outer radius of the crystal at z (the front-edge fillet is a quarter circle).
  auto rmax_at = [=]( const double z ) -> double {
    if( (r_b <= 0.0) || (z >= r_b) )
      return R;
    if( z <= 0.0 )
      return rho_c;
    const double dz = r_b - z;
    return rho_c + std::sqrt( std::max( 0.0, r_b*r_b - dz*dz ) );
  };
  // Outer radius of the active volume (the fillet, offset inward by the dead layer).
  auto rmax_active_at = [=]( const double z ) -> double {
    if( (r_ba <= 0.0) || (z >= z_a0 + r_ba) )
      return R_a;
    if( z <= z_a0 )
      return rho_ca;
    const double dz = (z_a0 + r_ba) - z;
    return rho_ca + std::sqrt( std::max( 0.0, r_ba*r_ba - dz*dz ) );
  };
  // Inner radius carved by the bore at z (zero ahead of the apex; a hemisphere for a rounded tip).
  auto rmin_bore_at = [=]( const double z ) -> double {
    if( (R_bore <= 0.0) || (z < z_apex) )
      return 0.0;
    if( !rounded || (z >= z_apex + R_bore) )
      return R_bore;
    const double dz = (z_apex + R_bore) - z;
    return std::sqrt( std::max( 0.0, R_bore*R_bore - dz*dz ) );
  };

  // The z values a profile over [z0,z1] is sampled at: its ends, any interior edges, and the arcs
  //  (24 segments each) - a step at an interior z comes out as two planes (just before / just after).
  auto sample_zs = [=]( const double z0, const double z1, const std::vector<double> &edges ) -> std::vector<double> {
    std::vector<double> zs{ z0, z1 };
    add_knots( zs, edges, z0, z1 );
    auto add_arc = [&]( const double a, const double b ){
      if( b <= a )
        return;
      for( int i = 0; i <= 24; ++i )
      {
        const double z = a + (b - a)*i/24.0;
        if( (z >= z0) && (z <= z1) )
          zs.push_back( z );
      }
    };
    if( r_b > 0.0 )
      add_arc( 0.0, r_b );
    if( r_ba > 0.0 )
      add_arc( z_a0, z_a0 + r_ba );
    if( rounded && (R_bore > 0.0) )
      add_arc( z_apex, z_apex + R_bore );
    finish_knots( zs );
    return zs;
  };

  // Builds a profile from radius functions, emitting a plane pair wherever either radius steps.
  auto make_profile = [&]( const std::vector<double> &zs,
                           const std::function<double(double)> &rmin_fn,
                           const std::function<double(double)> &rmax_fn ) -> std::vector<Plane> {
    std::vector<Plane> profile;
    const double eps = 1e-9;
    for( size_t i = 0; i < zs.size(); ++i )
    {
      const double z = zs[i];
      const bool first = (i == 0), last = (i + 1 == zs.size());
      const bool interior = !first && !last;
      const double zb = interior ? (z - eps) : z;   //just before (interior planes only)
      const double za = interior ? (z + eps) : z;   //just after
      const Plane before{ z, std::max(0.0, rmin_fn(zb)), std::max(0.0, rmax_fn(zb)) };
      const Plane after{ z, std::max(0.0, rmin_fn(za)), std::max(0.0, rmax_fn(za)) };
      const bool steps = (std::fabs(before.rmin - after.rmin) > 1e-7)
                         || (std::fabs(before.rmax - after.rmax) > 1e-7);
      if( first )
        profile.push_back( after );
      else if( last )
        profile.push_back( before );
      else if( steps )
      {
        profile.push_back( before );
        profile.push_back( after );
      }else
      {
        profile.push_back( after );
      }
    }//for( size_t i = 0; i < zs.size(); ++i )

    for( Plane &p : profile )
      p.rmax = std::max( p.rmax, p.rmin );
    return profile;
  };

  const string crystal_name = material_name( gd.crystal_material_index );
  const double crystal_density = material_density( gd.crystal_material_index );

  // Crystal volumes (cm3): the whole solid, the bore, the active part, and the dead remainder.
  double v_solid = 0.0, v_active = 0.0, v_bore = 0.0, v_dead = 0.0;
  if( box )
  {
    v_solid = 4.0*R*hy*L;
    v_active = v_solid;
  }else
  {
    v_solid = kPi*R*R*L - fillet_volume( R, r_b );
    v_bore = bore_volume_between( R_bore, z_apex, L, rounded, 0.0, L );
    const double active_len = std::max( 0.0, z_a1 - z_a0 );
    // Only the part of the bore that falls inside the active slab comes out of the active volume;
    //  whatever is left of it is inside the dead layer, which is the remainder below.
    v_active = kPi*R_a*R_a*active_len - fillet_volume( R_a, r_ba )
               - bore_volume_between( R_bore, z_apex, L, rounded, z_a0, z_a1 );
    v_active = std::max( 0.0, v_active );
    v_dead = std::max( 0.0, v_solid - v_bore - v_active );
  }

  // --- regions, outermost first so the crystal and bore draw on top ---------------------------
  std::vector<Region> layer_regions;
  Region collimator_region;
  bool have_collimator = false;

  int drawn_layers = 0;   //what the tooltip numbers: a skipped zero-thickness row must not leave a gap
  double r_in = R;        //inner radius of the next layer (layers start at the crystal radius)
  double hy_in = hy;      //and the box's y half-extent
  double z_if = 0.0;      //inner front z of the next layer (the crystal face, then each layer's front)
  for( size_t k = 0; k < gd.layers.size(); ++k )
  {
    const ceelo::LayerSpec &spec = gd.layers[k];
    const double f = std::max( 0.0, spec.front_thickness_cm );
    const double s = std::max( 0.0, spec.side_thickness_cm );
    if( (f <= 0.0) && (s <= 0.0) )
      continue;

    const double r_out = r_in + s;
    const double hy_out = hy_in + s;
    const double z_of = z_if - f;

    Region region;
    region.id = "layer" + std::to_string( k );
    region.kind = void_material( spec.material_index ) ? "void" : "layer";
    if( (f > 0.0) && (s > 0.0) )
      region.profile = { {z_of, 0.0, r_out}, {z_if, 0.0, r_out}, {z_if, r_in, r_out}, {L, r_in, r_out} };
    else if( f > 0.0 )
      region.profile = { {z_of, 0.0, r_in}, {z_if, 0.0, r_in} };          //a front disk only
    else
      region.profile = { {z_if, r_in, r_out}, {L, r_in, r_out} };         //a side tube only

    // A cup: the outer cylinder (or box) less the cavity it wraps.
    if( box )
      region.volume_cm3 = 4.0*r_out*hy_out*(L - z_of) - 4.0*r_in*hy_in*(L - z_if);
    else
      region.volume_cm3 = kPi*r_out*r_out*(L - z_of) - kPi*r_in*r_in*(L - z_if);
    region.volume_cm3 = std::max( 0.0, region.volume_cm3 );
    region.mass_g = region.volume_cm3 * material_density( spec.material_index );

    const bool is_void = (region.kind == "void");
    const string mat = material_name( spec.material_index );
    const WString mat_txt = (is_void && mat.empty()) ? text("dgd-vacuum") : WString::fromUTF8(mat);
    region.tip_title = text("dgd-layer").arg( ++drawn_layers ).arg( mat_txt ).toUTF8();
    if( !is_void )
      region.tip_lines.push_back( text("dgd-material").arg( mat_txt )
                                    .arg( compact_str( material_density(spec.material_index) ) ).toUTF8() );
    region.tip_lines.push_back( text("dgd-layer-dims").arg( length_str(f) ).arg( length_str(s) )
                                  .arg( length_str(2.0*r_out) ).toUTF8() );
    region.tip_lines.push_back( text("dgd-volume").arg( compact_str(region.volume_cm3) ).toUTF8() );
    // A void has no mass worth quoting (`vacuumMaterialSpec` gives it 1e-25 g/cm3 so the transport
    //  code has something to work with), so say it is empty instead - as the bore does.
    region.tip_lines.push_back( is_void ? text("dgd-void").toUTF8()
                                        : text("dgd-mass").arg( mass_str(region.mass_g) ).toUTF8() );
    layer_regions.push_back( region );

    r_in = r_out;
    hy_in = hy_out;
    z_if = z_of;
  }//for( each layer )

  if( gd.collimator && (gd.collimator->side_thickness_cm > 0.0) )
  {
    const ceelo::CollimatorSpec &spec = *gd.collimator;
    const double t_c = spec.side_thickness_cm;
    const double r_out = r_in + t_c;
    const double hy_out = hy_in + t_c;
    // The tube always spans the whole endcap: it starts at its own extension past the face, or at
    //  the outermost layer's front, whichever is further forward (as the ray trace has it).
    const double z_cs = std::min( spec.z_start_cm, z_if );

    collimator_region.id = "collimator";
    collimator_region.kind = void_material( spec.material_index ) ? "void" : "collimator";
    collimator_region.profile = { {z_cs, r_in, r_out}, {L, r_in, r_out} };
    if( box )
      collimator_region.volume_cm3 = 4.0*(r_out*hy_out - r_in*hy_in)*(L - z_cs);
    else
      collimator_region.volume_cm3 = kPi*(r_out*r_out - r_in*r_in)*(L - z_cs);
    collimator_region.mass_g = collimator_region.volume_cm3 * material_density( spec.material_index );

    const bool coll_void = (collimator_region.kind == "void");
    const string mat = material_name( spec.material_index );
    const WString mat_txt = (coll_void && mat.empty()) ? text("dgd-vacuum") : WString::fromUTF8(mat);
    collimator_region.tip_title = text("dgd-collimator").arg( mat_txt ).toUTF8();
    if( !coll_void )
      collimator_region.tip_lines.push_back( text("dgd-material").arg( mat_txt )
                                               .arg( compact_str( material_density(spec.material_index) ) ).toUTF8() );
    // The extension the user entered, measured from the detector face - not from z_cs, which is
    //  pulled forward to the endcap front when the endcap reaches further than the collimator.
    collimator_region.tip_lines.push_back( text("dgd-coll-dims").arg( length_str(t_c) )
                                             .arg( length_str( std::max(0.0, -spec.z_start_cm) ) )
                                             .arg( length_str(2.0*r_in) ).toUTF8() );
    collimator_region.tip_lines.push_back( text("dgd-volume").arg( compact_str(collimator_region.volume_cm3) ).toUTF8() );
    collimator_region.tip_lines.push_back( coll_void ? text("dgd-void").toUTF8()
                            : text("dgd-mass").arg( mass_str(collimator_region.mass_g) ).toUTF8() );
    have_collimator = true;

    r_in = r_out;
    z_if = z_cs;
  }//if( collimator )

  if( have_collimator )
    model.regions.push_back( collimator_region );
  for( size_t k = layer_regions.size(); k > 0; --k )
    model.regions.push_back( layer_regions[k-1] );   //outermost first

  // Dead layer: an L (front slab + side band, and a back slab when there is one) as one profile -
  //  valid because [rmin, rmax] is unique at every z.
  if( !box && ((t_f > 0.0) || (t_s > 0.0) || (t_b > 0.0)) )
  {
    Region dead;
    dead.id = "dead";
    dead.kind = "dead";
    auto dead_rmin = [=]( const double z ) -> double {
      if( (z < z_a0) || (z > z_a1) )
        return rmin_bore_at( z );      //front / back slab: the bore, if it reaches this far
      return rmax_active_at( z );      //side band: outside the active volume
    };
    const std::vector<double> zs = sample_zs( 0.0, L, { z_a0, z_a1, z_apex, z_apex + R_bore, r_b } );
    dead.profile = make_profile( zs, dead_rmin, rmax_at );
    dead.volume_cm3 = v_dead;
    dead.mass_g = v_dead * crystal_density;
    dead.tip_title = text("dgd-dead").toUTF8();
    dead.tip_lines.push_back( text("dgd-material").arg( WString::fromUTF8(crystal_name) )
                                .arg( compact_str(crystal_density) ).toUTF8() );
    if( t_b > 0.0 )
      dead.tip_lines.push_back( text("dgd-dead-dims-back").arg( length_str(t_f) ).arg( length_str(t_s) )
                                  .arg( length_str(t_b) ).toUTF8() );
    else
      dead.tip_lines.push_back( text("dgd-dead-dims").arg( length_str(t_f) ).arg( length_str(t_s) ).toUTF8() );
    dead.tip_lines.push_back( text("dgd-volume").arg( compact_str(v_dead) ).toUTF8() );
    dead.tip_lines.push_back( text("dgd-mass").arg( mass_str(dead.mass_g) ).toUTF8() );
    model.regions.push_back( dead );
  }//if( dead layer )

  // The active crystal.
  {
    Region crystal;
    crystal.id = "crystal";
    crystal.kind = "crystal";
    if( box )
    {
      crystal.profile = { {0.0, 0.0, R}, {L, 0.0, R} };
    }else
    {
      const std::vector<double> zs = sample_zs( z_a0, z_a1, { z_apex, z_apex + R_bore } );
      crystal.profile = make_profile( zs, rmin_bore_at, rmax_active_at );
    }
    crystal.volume_cm3 = v_active;
    crystal.mass_g = v_active * crystal_density;
    crystal.tip_title = text("dgd-crystal").toUTF8();
    crystal.tip_lines.push_back( text("dgd-material").arg( WString::fromUTF8(crystal_name) )
                                   .arg( compact_str(crystal_density) ).toUTF8() );
    if( box )
      crystal.tip_lines.push_back( text("dgd-dims-box").arg( length_str(2.0*R) ).arg( length_str(2.0*hy) )
                                     .arg( length_str(L) ).toUTF8() );
    else
      crystal.tip_lines.push_back( text("dgd-dims-cyl").arg( length_str(2.0*R) ).arg( length_str(L) ).toUTF8() );
    if( r_b > 0.0 )
      crystal.tip_lines.push_back( text("dgd-fillet").arg( length_str(r_b) ).toUTF8() );
    crystal.tip_lines.push_back( text("dgd-volume").arg( compact_str(v_active) ).toUTF8() );
    crystal.tip_lines.push_back( text("dgd-mass").arg( mass_str(crystal.mass_g) ).toUTF8() );
    model.regions.push_back( crystal );
  }

  // The bore (empty).
  if( R_bore > 0.0 )
  {
    Region bore;
    bore.id = "bore";
    bore.kind = "void";
    const std::vector<double> zs = sample_zs( z_apex, L, { z_apex + R_bore } );
    bore.profile = make_profile( zs, []( double ){ return 0.0; }, rmin_bore_at );
    bore.volume_cm3 = v_bore;
    bore.mass_g = 0.0;
    bore.tip_title = text("dgd-bore").toUTF8();
    bore.tip_lines.push_back( text("dgd-bore-dims").arg( length_str(2.0*R_bore) ).arg( length_str(D_bore) )
                                .arg( rounded ? text("dgd-bore-rounded") : WString() ).toUTF8() );
    bore.tip_lines.push_back( text("dgd-volume").arg( compact_str(v_bore) ).toUTF8() );
    bore.tip_lines.push_back( text("dgd-void").toUTF8() );
    model.regions.push_back( bore );
  }//if( bore )

  // --- extents and knots -------------------------------------------------------------------------
  double z_min = 0.0, z_max = L, r_max = R;
  for( const Region &region : model.regions )
  {
    for( const Plane &p : region.profile )
    {
      z_min = std::min( z_min, p.z );
      z_max = std::max( z_max, p.z );
      r_max = std::max( r_max, p.rmax );
    }
  }//for( each region )

  // The arc samples are not knots (the JS scales them like any other point); only true edges are,
  //  and only of the features the geometry actually has.
  model.z_knots = { 0.0, L };
  model.r_knots = { 0.0, R };
  if( r_b > 0.0 )
  {
    model.z_knots.push_back( r_b );
    model.r_knots.push_back( rho_c );
  }
  if( t_f > 0.0 )
    model.z_knots.push_back( z_a0 );
  if( t_s > 0.0 )
    model.r_knots.push_back( R_a );
  if( r_ba > 0.0 )
  {
    model.z_knots.push_back( z_a0 + r_ba );
    model.r_knots.push_back( rho_ca );
  }
  if( t_b > 0.0 )
    model.z_knots.push_back( z_a1 );
  if( R_bore > 0.0 )
  {
    model.z_knots.push_back( z_apex );
    model.r_knots.push_back( R_bore );
    if( rounded )
      model.z_knots.push_back( z_apex + R_bore );
  }
  // Every layer/collimator edge is a true step in the geometry, so it seeds the min-pixel map (this
  //  is what keeps a 10 um window a visible 2 px).  The crystal/dead/bore profiles are sampled arcs
  //  whose points are not edges; their handful of real knots are the ones added by hand above.
  //  Iterate the source lists rather than filtering `model.regions` by kind: a vacuum layer's kind
  //  is "void", the same as the bore's.
  auto add_edge_knots = [&model]( const Region &region ){
    for( const Plane &p : region.profile )
    {
      model.z_knots.push_back( p.z );
      model.r_knots.push_back( p.rmin );
      model.r_knots.push_back( p.rmax );
    }
  };
  for( const Region &region : layer_regions )
    add_edge_knots( region );
  if( have_collimator )
    add_edge_knots( collimator_region );
  finish_knots( model.z_knots );
  finish_knots( model.r_knots );

  model.z_min = z_min;
  model.z_max = z_max;
  model.r_max = r_max;

  return model;
}//buildModel(...)


std::string DetectorGeometryDiagram::toJson( const Model &model )
{
  nlohmann::json j;
  j["box"] = model.box;
  j["zMin"] = model.z_min;
  j["zMax"] = model.z_max;
  j["rMax"] = model.r_max;
  j["knots"]["z"] = model.z_knots;
  j["knots"]["r"] = model.r_knots;
  j["labels"]["front"] = text("dgd-front").toUTF8();

  nlohmann::json regions = nlohmann::json::array();
  for( const Region &region : model.regions )
  {
    nlohmann::json r;
    r["id"] = region.id;
    r["kind"] = region.kind;
    nlohmann::json profile = nlohmann::json::array();
    for( const Plane &p : region.profile )
      profile.push_back( { p.z, p.rmin, p.rmax } );
    r["profile"] = profile;
    r["tip"]["title"] = region.tip_title;
    r["tip"]["lines"] = region.tip_lines;
    regions.push_back( r );
  }//for( each region )
  j["regions"] = regions;

  // `replace`: a material name that is not valid UTF-8 (an imported descriptor can carry one) would
  //  otherwise throw, and the drawing would silently freeze on its last content with nothing saying
  //  why.  A replacement character in a tooltip is the better failure.
  return j.dump( -1, ' ', false, nlohmann::json::error_handler_t::replace );
}//toJson(...)
