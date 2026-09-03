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

#include <array>
#include <cmath>
#include <string>
#include <vector>
#include <cstring>
#include <istream>
#include <algorithm>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/AngleOutxImport.h"

using namespace std;

/** Reader for the ANGLE 5 XML file family, as documented at
 https://www.angle.me/file-formats.html.

 Every ANGLE file type shares one grammar, so this file is organized by that
 grammar rather than by file type: one reader per element, and a single
 `AngleOutx::parse` that calls whichever ones the document happens to contain.
 */

namespace
{
using rapidxml::xml_node;
using rapidxml::xml_attribute;

/** A numeric attribute plus whether it was actually present, so a legitimate
 zero (ANGLE's "float0") is distinguishable from an absent attribute - which is
 how per-detector-type attribute sets are detected without branching. */
struct Attr
{
  double value = 0.0;
  bool present = false;
};//struct Attr


const xml_attribute<char> *find_attr( const xml_node<char> * const node, const char * const name )
{
  return node ? node->first_attribute( name, strlen(name), true ) : nullptr;
}//find_attr(...)


/** An integer attribute, bounded so `std::lround` cannot be handed something it
 has no defined answer for. */
int attr_int( const xml_node<char> * const node, const char * const name,
              const int max_value );

/** Numeric attribute; throws if present but not a number (a malformed file
 should not be silently read as zero). */
Attr attr_num( const xml_node<char> * const node, const char * const name )
{
  Attr answer;
  const xml_attribute<char> * const a = find_attr( node, name );
  if( !a || !a->value_size() )
    return answer;

  const string txt = SpecUtils::trim_copy( SpecUtils::xml_value_str(a) );

  // `SpecUtils::parse_double` stops at the first character it cannot use and
  //  still reports success, so "29,15" would silently become 29.  ANGLE's
  //  decimal separator is always a period, whatever the writer's locale, so
  //  anything else in the token is an error - not a value to salvage.
  for( const char c : txt )
  {
    const bool numeric = (c >= '0' && c <= '9') || (c == '.') || (c == '+')
                          || (c == '-') || (c == 'e') || (c == 'E');
    if( !numeric )
      throw runtime_error( "ANGLE: attribute '" + string(name) + "' is not a number ('"
                           + txt + "')." + ((c == ',') ? "  ANGLE files always use a"
                           " period as the decimal separator." : "") );
  }//for( const char c : txt )

  // parse_double, not parse_float: ANGLE writes its solid angles and summing
  //  corrections to ~15 significant figures, and float would quietly drop half
  //  of them.
  double val = 0.0;
  if( !SpecUtils::parse_double( txt.c_str(), txt.size(), val ) )
    throw runtime_error( "ANGLE: attribute '" + string(name) + "' is not a number ('"
                         + txt + "')." );

  if( IsNan(val) || IsInf(val) )
    throw runtime_error( "ANGLE: attribute '" + string(name) + "' is not a finite number ('"
                         + txt + "')." );

  answer.value = val;
  answer.present = true;
  return answer;
}//attr_num(...)


int attr_int( const xml_node<char> * const node, const char * const name,
              const int max_value )
{
  const Attr attr = attr_num( node, name );   //already rejects NaN / infinity
  if( !attr.present )
    return 0;

  if( (attr.value < 0.0) || (attr.value > max_value) )
    throw runtime_error( "ANGLE: attribute '" + string(name) + "' is out of range ("
                         + std::to_string(attr.value) + ")." );

  return static_cast<int>( std::lround( attr.value ) );
}//attr_int(...)


/** Length attribute, converted from the file's units into PhysicalUnits. */
Attr attr_len( const xml_node<char> * const node, const char * const name, const double unit )
{
  Attr answer = attr_num( node, name );
  answer.value *= unit;
  return answer;
}//attr_len(...)


/** ANGLE writes yes/no; accept the other usual truthy spellings so a
 hand-edited or future-version file isn't silently read as false. */
bool attr_bool( const xml_node<char> * const node, const char * const name )
{
  const string val = SpecUtils::xml_value_str( find_attr(node, name) );
  return SpecUtils::iequals_ascii( val, "yes" )
         || SpecUtils::iequals_ascii( val, "true" )
         || (val == "1");
}//attr_bool(...)


string attr_str( const xml_node<char> * const node, const char * const name )
{
  return SpecUtils::xml_value_str( find_attr(node, name) );
}//attr_str(...)


const xml_node<char> *child( const xml_node<char> * const node, const char * const name )
{
  return node ? node->first_node( name, strlen(name), true ) : nullptr;
}//child(...)


/** Text content of an element, as a number; returns false if absent/unparsable. */
bool node_num( const xml_node<char> * const node, double &value )
{
  if( !node || !node->value_size() )
    return false;
  return SpecUtils::parse_double( node->value(), node->value_size(), value );
}//node_num(...)


/** Scales `fracs` so they sum to one; a file whose mass fractions do not add up
 to 100 (rounding, or a hand edit) should still give a sane composition. */
template<class T>
bool normalize_fractions( vector<T> &fracs )
{
  double sum = 0.0;
  for( const T &f : fracs )
    sum += f.massFraction;

  if( sum <= 0.0 )
    return false;

  for( T &f : fracs )
    f.massFraction /= sum;

  // The spec requires the percentages to add to 100; rescaling keeps a slightly
  //  off file usable, but a badly off one is a hand edit worth reporting.
  return (std::fabs(sum - 1.0) <= 0.005);
}//normalize_fractions(...)


/** Reads a `<material>` element in any of its four forms: a bare predefined
 reference, a mixture of elements, a single compound, or a mixture of
 compounds.  This is the one place material syntax is understood; every
 context that can hold a material calls it. */
AngleMaterial read_material( const xml_node<char> * const mat_node )
{
  AngleMaterial answer;
  if( !mat_node )
    return answer;

  answer.name = attr_str( mat_node, "name" );
  answer.density_g_cm3 = attr_num( mat_node, "density" ).value;

  // <elements> - mass percentages summing to 100.
  const xml_node<char> * const elems_node = child( mat_node, "elements" );
  if( elems_node )
  {
    XML_FOREACH_CHILD( elem, elems_node, "element" )
    {
      AngleElementFraction frac;
      frac.symbol = attr_str( elem, "symbol" );
      frac.massFraction = attr_num( elem, "massFraction" ).value / 100.0;
      if( !frac.symbol.empty() && (frac.massFraction > 0.0) )
        answer.elements.push_back( frac );
    }//XML_FOREACH_CHILD( elem, elems_node, "element" )

    if( !normalize_fractions( answer.elements ) )
      answer.fractionsRescaled = true;
  }//if( elems_node )

  // A single <compound> (no massFraction attribute), or a <compounds> mixture.
  auto read_compound = []( const xml_node<char> * const cmp_node, const bool weighted ) -> AngleCompound
  {
    AngleCompound cmp;
    cmp.chemicalFormula = attr_str( cmp_node, "chemicalFormula" );
    const Attr frac = attr_num( cmp_node, "massFraction" );
    cmp.massFraction = (weighted && frac.present) ? (frac.value / 100.0) : 1.0;

    XML_FOREACH_CHILD( elem, cmp_node, "element" )
    {
      AngleCompoundElement ce;
      ce.symbol = attr_str( elem, "symbol" );
      ce.atoms = attr_int( elem, "atoms", 1000 );
      if( !ce.symbol.empty() && (ce.atoms > 0) )
        cmp.elements.push_back( ce );
    }//XML_FOREACH_CHILD( elem, cmp_node, "element" )

    return cmp;
  };//read_compound lambda

  const xml_node<char> * const cmps_node = child( mat_node, "compounds" );
  if( cmps_node )
  {
    XML_FOREACH_CHILD( cmp_node, cmps_node, "compound" )
    {
      const AngleCompound cmp = read_compound( cmp_node, true );
      if( !cmp.elements.empty() && (cmp.massFraction > 0.0) )
        answer.compounds.push_back( cmp );
    }//XML_FOREACH_CHILD( cmp_node, cmps_node, "compound" )

    if( !normalize_fractions( answer.compounds ) )
      answer.fractionsRescaled = true;
  }else
  {
    const xml_node<char> * const cmp_node = child( mat_node, "compound" );
    if( cmp_node )
    {
      const AngleCompound cmp = read_compound( cmp_node, false );
      if( !cmp.elements.empty() )
        answer.compounds.push_back( cmp );
    }
  }//if( cmps_node ) / else

  return answer;
}//read_material(...)


/** The common `<x><material/></x>` shape: a layer's material is a child
 element, not an attribute of the layer. */
AngleMaterial read_child_material( const xml_node<char> * const parent )
{
  return read_material( child(parent, "material") );
}//read_child_material(...)


/** Thickness of a layer on the crystal's OUTWARD-facing top surface.

 Most detector types spell it `topThickness`.  Well types (6, 8) have an annular
 crystal with two distinct axial faces and spell them `topUpperThickness` (the
 outward one) and `topLowerThickness` (the well floor, facing the sample down in
 the bore).  They are two different surfaces, so they must not be added
 together; only the outward one is on the path from a source in front of the
 detector. */
double outward_front( const xml_node<char> * const node, const double unit )
{
  const Attr top = attr_len( node, "topThickness", unit );
  return top.present ? top.value : attr_len( node, "topUpperThickness", unit ).value;
}//outward_front(...)


/** Counterpart of #outward_front for the side: `sideThickness`, or a well
 crystal's `sideOuterThickness` (the outer wall; `sideInnerThickness` lines the
 well itself). */
double outward_side( const xml_node<char> * const node, const double unit )
{
  const Attr side = attr_len( node, "sideThickness", unit );
  return side.present ? side.value : attr_len( node, "sideOuterThickness", unit ).value;
}//outward_side(...)


void add_layer( vector<AngleOutxContents::Layer> &layers, const AngleLayerRole role,
                const AngleMaterial &material, const double front, const double side,
                const bool behind = false )
{
  // The spec's "float0" means non-negative; a negative thickness would survive
  //  into the geometry as a layer that un-recesses the crystal.
  if( (front < 0.0) || (side < 0.0) )
    throw runtime_error( "ANGLE: negative thickness on the detector's "
                         + string(to_string(role)) + "." );

  if( (front <= 0.0) && (side <= 0.0) )
    return;

  AngleOutxContents::Layer layer;
  layer.role = role;
  layer.material = material;
  layer.front = front;
  layer.side = side;
  layer.behindCrystal = behind;
  layers.push_back( layer );
}//add_layer(...)


AngleContainer read_container( const xml_node<char> * const node, const double unit )
{
  AngleContainer answer;
  if( !node )
    return answer;

  answer.present = true;
  answer.type = attr_str( node, "type" );
  answer.name = attr_str( node, "name" );
  answer.description = attr_str( node, "description" );

  const xml_node<char> * const shape = child( node, "shape" );
  answer.innerRadius = attr_len( shape, "innerRadius", unit ).value;
  answer.sideThickness = attr_len( shape, "sideThickness", unit ).value;
  answer.footHeight = attr_len( shape, "footHeight", unit ).value;

  // The published spec spells this "botomThickness"; real files write
  //  "bottomThickness".  Accept either.
  const Attr bottom = attr_len( shape, "bottomThickness", unit );
  answer.bottomThickness = bottom.present ? bottom.value
                                          : attr_len( shape, "botomThickness", unit ).value;

  answer.cavityRadius = attr_len( shape, "cavityRadius", unit ).value;
  answer.cavityDepth = attr_len( shape, "cavityDepth", unit ).value;
  answer.sideInnerThickness = attr_len( shape, "sideInnerThickness", unit ).value;
  answer.bottomUpperThickness = attr_len( shape, "bottomUpperThickness", unit ).value;
  answer.bottomLowerThickness = attr_len( shape, "bottomLowerThickness", unit ).value;

  answer.material = read_child_material( node );

  const xml_node<char> * const coatings = child( node, "coatings" );
  XML_FOREACH_CHILD( coating, coatings, "coating" )
  {
    AngleContainer::Coating c;
    c.side = attr_len( coating, "sideThickness", unit ).value;
    c.bottom = attr_len( coating, "bottomThickness", unit ).value;
    c.bottomUpper = attr_len( coating, "bottomUpperThickness", unit ).value;
    c.bottomLower = attr_len( coating, "bottomLowerThickness", unit ).value;
    c.material = read_child_material( coating );
    answer.coatings.push_back( c );
  }//XML_FOREACH_CHILD( coating, coatings, "coating" )

  return answer;
}//read_container(...)


AngleHolderGeometry read_geometry( const xml_node<char> * const node, const double unit )
{
  AngleHolderGeometry answer;
  if( !node )
    return answer;

  answer.present = true;
  answer.type = attr_str( node, "type" );
  answer.name = attr_str( node, "name" );
  answer.description = attr_str( node, "description" );

  const xml_node<char> * const holder = child( node, "holder" );
  if( holder )
  {
    answer.holderOuterRadius = attr_len( holder, "outerRadius", unit ).value;
    answer.holderHeight = attr_len( holder, "height", unit ).value;

    const xml_node<char> * const cap = child( holder, "cap" );
    answer.capThickness = attr_len( cap, "thickness", unit ).value;
    answer.capMaterial = read_child_material( cap );

    const xml_node<char> * const wall = child( holder, "wall" );
    answer.wallThickness = attr_len( wall, "thickness", unit ).value;
    answer.wallMaterial = read_child_material( wall );
  }//if( holder )

  const xml_node<char> * const abs_layers = child( node, "absorbingLayers" );
  XML_FOREACH_CHILD( abs_layer, abs_layers, "absorbingLayer" )
  {
    AngleHolderGeometry::AbsorbingLayer layer;
    layer.top = attr_len( abs_layer, "topThickness", unit ).value;
    layer.bottom = attr_len( abs_layer, "bottomThickness", unit ).value;
    layer.side = attr_len( abs_layer, "sideThickness", unit ).value;
    layer.material = read_child_material( abs_layer );
    if( (layer.top > 0.0) || (layer.bottom > 0.0) || (layer.side > 0.0) )
      answer.absorbingLayers.push_back( layer );
  }//XML_FOREACH_CHILD( abs_layer, abs_layers, "absorbingLayer" )

  return answer;
}//read_geometry(...)


AngleSource read_source( const xml_node<char> * const node, const double unit )
{
  AngleSource answer;
  if( !node )
    return answer;

  answer.present = true;
  answer.radius = attr_len( node, "radius", unit ).value;
  answer.height = attr_len( node, "height", unit ).value;
  answer.material = read_child_material( node );

  return answer;
}//read_source(...)


/** Reads `<detector>` - the crystal, its dead layers and bore, and the stack of
 concentric layers around it.

 Every known attribute spelling is read unconditionally: ANGLE only writes the
 ones its detector type defines, and `Attr::present` tells them apart.  Which
 spellings mean what for the declared type is settled once, afterwards, in
 `applyTypeSemantics()`. */
void read_detector( const xml_node<char> * const det_node, const double unit,
                    AngleOutxContents &contents )
{
  contents.detName = attr_str( det_node, "name" );
  contents.detDescription = attr_str( det_node, "description" );
  contents.detType = attr_str( det_node, "type" );
  contents.detTypeEnum = angleDetectorTypeFromString( contents.detType );

  const xml_node<char> * const crystal = child( det_node, "crystal" );
  contents.crystalRadius = attr_len( crystal, "radius", unit ).value;
  contents.crystalLength = attr_len( crystal, "height", unit ).value;
  contents.bulletizingRadius = attr_len( crystal, "bulletizingRadius", unit ).value;

  const xml_node<char> * const core = child( det_node, "core" );
  if( core )
  {
    contents.hasCore = true;
    contents.coreRadius = attr_len( core, "radius", unit ).value;
    contents.coreDepth = attr_len( core, "height", unit ).value;
    contents.coreRounded = attr_bool( core, "rounded" );
  }//if( core )

  const xml_node<char> * const well = child( det_node, "well" );
  if( well )
  {
    contents.well.present = true;
    contents.well.depth = attr_len( well, "depth", unit ).value;
    contents.well.radius = attr_len( well, "radius", unit ).value;
  }//if( well )

  const xml_node<char> * const inactive = child( det_node, "inactiveGe" );
  contents.deadLayerFront = attr_len( inactive, "topThickness", unit ).value;
  contents.deadLayerSide = attr_len( inactive, "sideThickness", unit ).value;
  contents.deadLayerBack = attr_len( inactive, "bottomThickness", unit ).value;
  contents.deadLayerTopUpper = attr_len( inactive, "topUpperThickness", unit ).value;
  contents.deadLayerTopLower = attr_len( inactive, "topLowerThickness", unit ).value;
  contents.deadLayerSideInner = attr_len( inactive, "sideInnerThickness", unit ).value;
  contents.deadLayerSideOuter = attr_len( inactive, "sideOuterThickness", unit ).value;

  const xml_node<char> * const contact = child( det_node, "contact" );
  if( contact )
  {
    contents.contact.front = attr_len( contact, "topThickness", unit ).value;
    contents.contact.side = attr_len( contact, "sideThickness", unit ).value;
    contents.contact.material = read_child_material( contact );

    // An implanted or Li-diffused contact is dead crystal, not a separate
    //  attenuator: fold it into the dead layer.  A thick contact of some other
    //  material is a real layer.  10 um is comfortably above any implanted
    //  contact (they are ~0.3 um) and below any deposited one.
    const string &mat = contents.contact.material.name;
    const double thickest = std::max( contents.contact.front, contents.contact.side );
    contents.contact.readsAsDeadCrystal = mat.empty()
                    || SpecUtils::iequals_ascii( mat, "Germanium" )
                    || SpecUtils::iequals_ascii( mat, "Ge" )
                    || SpecUtils::iequals_ascii( mat, "Lithium" )
                    || SpecUtils::iequals_ascii( mat, "Li" )
                    || (thickest < (0.01 * PhysicalUnits::mm));
  }//if( contact )

  const xml_node<char> * const pin = child( det_node, "contactPin" );
  if( pin )
  {
    contents.contactPin.radius = attr_len( pin, "radius", unit ).value;
    contents.contactPin.material = read_child_material( pin );
  }//if( pin )

  // Vacuum: several attributes, whose split depends on whether an
  //  antimicrophonic shield sits inside the gap.  Keep the split, and the
  //  totals - only the totals matter when there is no shield.
  const xml_node<char> * const vacuum = child( det_node, "vacuum" );
  const xml_node<char> * const shield = child( det_node, "antimicrophonicShield" );
  {
    const double top = attr_len( vacuum, "topThickness", unit ).value;
    const double top_upper = attr_len( vacuum, "topUpperThickness", unit ).value;
    const double bottom = attr_len( vacuum, "bottomThickness", unit ).value;
    const double bottom_upper = attr_len( vacuum, "bottomUpperThickness", unit ).value;
    const double side = attr_len( vacuum, "sideThickness", unit ).value;
    const double side_inner = attr_len( vacuum, "sideInnerThickness", unit ).value;
    const double side_outer = attr_len( vacuum, "sideOuterThickness", unit ).value;

    // The spec's superscripts form a 2x2: "top"/"bottom" pair the way
    //  "sideOuter"/"sideInner" do, and the second of each pair is written ONLY
    //  when an antimicrophonic shield exists - which is what a shield does, it
    //  splits the gap.  So "top"/"sideOuter" is the gap outside the shield and
    //  "bottom"/"sideInner" the gap inside it.
    //
    //  Well types spell the same two positions "topUpper"/"bottomUpper"; their
    //  "topLower"/"bottomLower" are the WELL FLOOR's gap - a different surface,
    //  not part of the path from a source in front of the detector - so they are
    //  deliberately not added in (see outward_front()).
    contents.vacuumFrontOuter = top + top_upper;
    contents.vacuumFrontInner = bottom + bottom_upper;
    contents.vacuumSideOuter = side + side_outer;
    contents.vacuumSideInner = side_inner;

    if( !shield )
    {
      // No shield: the split is meaningless, everything is one gap.
      contents.vacuumFrontOuter += contents.vacuumFrontInner;
      contents.vacuumSideOuter += contents.vacuumSideInner;
      contents.vacuumFrontInner = contents.vacuumSideInner = 0.0;
    }

    contents.vacuumFront = contents.vacuumFrontInner + contents.vacuumFrontOuter;
    contents.vacuumSide = contents.vacuumSideInner + contents.vacuumSideOuter;
  }

  // --- the concentric stack, innermost -> outermost ---

  // A contact that is not dead crystal sits directly on the crystal.
  if( !contents.contact.readsAsDeadCrystal )
    add_layer( contents.layers, AngleLayerRole::Contact, contents.contact.material,
               contents.contact.front, contents.contact.side );

  // NaI reflecting layer, wrapped directly around the crystal.
  const xml_node<char> * const reflecting = child( det_node, "reflectingLayer" );
  if( reflecting )
  {
    // For a well crystal "topUpper"/"topLower" and "sideOuter"/"sideInner" are
    //  the OUTWARD-facing and WELL-facing surfaces, not a stack - see
    //  outward_front()/outward_side() below.
    const double front = outward_front( reflecting, unit );
    const double side = outward_side( reflecting, unit );
    add_layer( contents.layers, AngleLayerRole::ReflectingLayer,
               read_child_material(reflecting), front, side );
  }//if( reflecting )

  // Housing - the cup the crystal sits in, inside the cryostat vacuum.
  const xml_node<char> * const housing = child( det_node, "housing" );
  if( housing )
  {
    const xml_node<char> * const side_in = child( housing, "sideInner" );
    const xml_node<char> * const side_out = child( housing, "sideOuter" );
    const xml_node<char> * const top_lo = child( housing, "topLower" );
    const xml_node<char> * const top_up = child( housing, "topUpper" );

    add_layer( contents.layers, AngleLayerRole::HousingInner, read_child_material(top_lo),
               attr_len( top_lo, "thickness", unit ).value, 0.0 );
    add_layer( contents.layers, AngleLayerRole::HousingInner, read_child_material(side_in),
               0.0, attr_len( side_in, "thickness", unit ).value );
    add_layer( contents.layers, AngleLayerRole::HousingOuter, read_child_material(top_up),
               attr_len( top_up, "thickness", unit ).value, 0.0 );
    add_layer( contents.layers, AngleLayerRole::HousingOuter, read_child_material(side_out),
               0.0, attr_len( side_out, "thickness", unit ).value );
  }//if( housing )

  // The cryostat vacuum gap, with the antimicrophonic shield (if any) sandwiched
  //  inside it.  The gap is emitted as a layer with no material: it is a
  //  spacing, not an attenuator, but it has to hold its place in the stack so
  //  everything outside it sits at the right depth / radius.
  add_layer( contents.layers, AngleLayerRole::VacuumInner, AngleMaterial{},
             contents.vacuumFrontInner, contents.vacuumSideInner );

  if( shield )
  {
    const double front = outward_front( shield, unit );
    const double side = attr_len( shield, "sideThickness", unit ).value;
    add_layer( contents.layers, AngleLayerRole::AntimicrophonicShield,
               read_child_material(shield), front, side );
  }//if( shield )

  add_layer( contents.layers, AngleLayerRole::VacuumOuter, AngleMaterial{},
             contents.vacuumFrontOuter, contents.vacuumSideOuter );

  const xml_node<char> * const endcap = child( det_node, "endCap" );
  if( endcap )
  {
    // <window> (optional; not present in all files) - a thin full-face window
    //  let into the endcap front over a hole of radius `holeRadius`.  When
    //  holeRadius >= crystal radius the window replaces the endcap front
    //  on-axis (the endcap has been bored out there), so on-axis the true front
    //  is only the window - counting BOTH the endcap top and the window
    //  over-recesses the crystal.  In that case take the endcap SIDE only, and
    //  the window as the front; otherwise take the endcap front (an intact
    //  endcap, no on-axis hole).
    const xml_node<char> * const window = child( endcap, "window" );
    const double endcap_front = outward_front( endcap, unit );
    const double endcap_side = attr_len( endcap, "sideThickness", unit ).value;

    bool window_replaces_front = false;
    if( window )
    {
      contents.endCapWindowRadius = attr_len( window, "radius", unit ).value;
      contents.endCapWindowHoleRadius = attr_len( window, "holeRadius", unit ).value;
      contents.endCapWindowHolderThickness = attr_len( window, "holderThickness", unit ).value;

      // The end cap is bored out behind the window over `holeRadius`.  When the
      //  hole covers the whole crystal, the on-axis front is the window alone -
      //  counting the end cap too would recess the crystal twice.
      const double hole_radius = contents.endCapWindowHoleRadius;
      window_replaces_front = (hole_radius > 0.0) && (hole_radius >= contents.crystalRadius);

      // A hole narrower than the crystal makes the front thickness depend on
      //  radius, which a concentric layer cannot express: on-axis it is the
      //  window alone, off-axis the end cap as well.  Taking both is the
      //  conservative choice, but it over-attenuates exactly on-axis, which is
      //  where a response is anchored - so say so rather than hide it.
      if( (hole_radius > 0.0) && !window_replaces_front )
        contents.parseNotes.push_back( "This detector's end-cap window covers only the middle of"
                                       " the crystal, so how much end cap a gamma ray crosses"
                                       " depends on where it enters.  Both the end cap and the"
                                       " window were counted, which slightly over-estimates the"
                                       " on-axis attenuation." );

      if( contents.endCapWindowHolderThickness > 0.0 )
        contents.parseNotes.push_back( "The ring that clamps this detector's end-cap window is"
                                       " not modeled; it attenuates off-axis only." );
    }//if( window )

    add_layer( contents.layers, AngleLayerRole::EndCap, read_child_material(endcap),
               window_replaces_front ? 0.0 : endcap_front, endcap_side );

    if( window )
      add_layer( contents.layers, AngleLayerRole::EndCapWindow, read_child_material(window),
                 attr_len( window, "thickness", unit ).value, 0.0 );

    const xml_node<char> * const coatings = child( endcap, "coatings" );
    XML_FOREACH_CHILD( coating, coatings, "coating" )
    {
      const double front = outward_front( coating, unit );
      const double side = attr_len( coating, "sideThickness", unit ).value;
      add_layer( contents.layers, AngleLayerRole::EndCapCoating,
                 read_child_material(coating), front, side );
    }//XML_FOREACH_CHILD( coating, coatings, "coating" )
  }//if( endcap )

  // Pieces behind the crystal (NaI light path).  Carried for fidelity; they
  //  affect backscatter only, and CeeLo models no back layers.
  const xml_node<char> * const coupling = child( det_node, "opticalCoupling" );
  if( coupling )
    add_layer( contents.layers, AngleLayerRole::OpticalCoupling,
               read_child_material(coupling),
               attr_len( coupling, "thickness", unit ).value, 0.0, true );

  const xml_node<char> * const pmt = child( det_node, "photomultiplierTube" );
  if( pmt )
  {
    const xml_node<char> * const wall = child( pmt, "wall" );
    const xml_node<char> * const window = child( pmt, "window" );
    add_layer( contents.layers, AngleLayerRole::PmtWall, read_child_material(wall),
               0.0, attr_len( wall, "thickness", unit ).value, true );
    add_layer( contents.layers, AngleLayerRole::PmtWindow, read_child_material(window),
               attr_len( window, "thickness", unit ).value, 0.0, true );
  }//if( pmt )

  contents.hasGeometry = (contents.crystalRadius > 0.0) && (contents.crystalLength > 0.0);
}//read_detector(...)


/** Applies the per-detector-type meaning of the attributes read above, and
 decides whether the result can become a `ceelo::GeometryDescriptor`.  Kept as
 one pass so there is a single place to edit when the ANGLE spec moves. */
void applyTypeSemantics( AngleOutxContents &contents )
{
  const AngleDetectorType type = contents.detTypeEnum;

  if( contents.detType.empty() )
    contents.parseNotes.push_back( "The ANGLE file does not say what type of detector this is;"
                                   " it was read as a closed-end coaxial HPGe." );
  else if( type == AngleDetectorType::Unknown )
    contents.parseNotes.push_back( "Unrecognized ANGLE detector type '" + contents.detType
                                   + "'; it was read as a closed-end coaxial HPGe." );

  // Well detectors spell the dead layer with Upper/Lower/Inner/Outer suffixes;
  //  collapse them onto the front/side pair the rest of the code uses.
  // A well crystal is annular, with two axial faces and two radial ones.  Only
  //  the outward-facing pair is on the path from a source in front of the
  //  detector, so those are the ones that become front/side; the well-facing
  //  pair is kept in its own fields and reported, not silently folded in.
  const bool is_well = (type == AngleDetectorType::Well) || (type == AngleDetectorType::NaIWell);
  if( is_well || (contents.deadLayerFront <= 0.0 && contents.deadLayerSide <= 0.0) )
  {
    if( contents.deadLayerTopUpper > 0.0 )
      contents.deadLayerFront = std::max( contents.deadLayerFront, contents.deadLayerTopUpper );
    if( contents.deadLayerSideOuter > 0.0 )
      contents.deadLayerSide = std::max( contents.deadLayerSide, contents.deadLayerSideOuter );

    if( (contents.deadLayerTopLower > 0.0) || (contents.deadLayerSideInner > 0.0) )
      contents.parseNotes.push_back( "The dead layer lining this detector's well is recorded but"
                                     " not modeled; only the outward-facing dead layer is used." );
  }//if( well-style dead layer )

  // A back dead layer is only defined for planar detectors.
  if( (contents.deadLayerBack > 0.0) && (type != AngleDetectorType::PlanarLEPD)
      && (type != AngleDetectorType::Unknown) )
    contents.parseNotes.push_back( "The file gives an <inactiveGe bottomThickness> for a"
                                   " non-planar detector; it was kept as a back dead layer." );

  // Can this become a CeeLo geometry?  See the comments on each refusal.
  contents.modeASupported = contents.hasGeometry;
  contents.modeAObstruction.clear();

  if( !contents.hasGeometry )
  {
    contents.modeAObstruction = "The file does not give the detector crystal's dimensions.";
  }else if( is_well || contents.well.present )
  {
    contents.modeASupported = false;
    contents.modeAObstruction = "This is a well detector.  InterSpec's response calculation has"
                                " no well-cavity shape, so the physical detector model cannot be"
                                " used; the measured efficiency curve can still be imported.";
  }else if( (type == AngleDetectorType::TrueCoaxHPGe)
            || (type == AngleDetectorType::OpenEndCoaxGeLi) )
  {
    contents.modeASupported = false;
    contents.modeAObstruction = "This is a true-coaxial detector, whose bore goes all the way"
                                " through the crystal.  InterSpec's response calculation only"
                                " models a blind bore, and the difference is largest on-axis;"
                                " the measured efficiency curve can still be imported.";
  }//if( ... ) / else if( ... )
}//applyTypeSemantics(...)


/** `<referenceEfficiencyCurve>`: the measured points, the fitted regions, and
 the geometry the reference measurement was made in. */
void read_reference_curve( const xml_node<char> * const angle_node, const double unit,
                           AngleOutxContents &contents )
{
  const xml_node<char> * const ref_node = child( angle_node, "referenceEfficiencyCurve" );
  if( !ref_node )
    return;

  contents.referenceCurveName = attr_str( ref_node, "name" );
  contents.referenceCurveDescription = attr_str( ref_node, "description" );
  contents.referenceDetectorName = attr_str( child(ref_node, "detector"), "name" );

  // Measured experimental points (energy keV, absolute efficiency).
  const xml_node<char> * const points_node = child( ref_node, "experimentalPoints" );
  XML_FOREACH_CHILD( point, points_node, "point" )
  {
    const Attr energy = attr_num( point, "energy" );
    const Attr eff = attr_num( point, "efficiency" );
    if( !energy.present || !eff.present )
      continue;

    if( (energy.value > 0.0) && (eff.value > 0.0) && !IsNan(eff.value) && !IsInf(eff.value) )
      contents.referencePoints.emplace_back( static_cast<float>(energy.value),
                                             static_cast<float>(eff.value) );
  }//XML_FOREACH_CHILD( point, points_node, "point" )

  std::sort( begin(contents.referencePoints), end(contents.referencePoints),
             []( const pair<float,float> &l, const pair<float,float> &r ) -> bool {
               return l.first < r.first;
             } );

  contents.referenceContainer = read_container( child(ref_node, "container"), unit );
  contents.referenceGeometry = read_geometry( child(ref_node, "geometry"), unit );
  contents.referenceSource = read_source( child(ref_node, "source"), unit );

  // Reference on-axis distance from the reference curve's own holder height
  //  (NOT the top-level near-contact holder).
  if( contents.referenceGeometry.holderHeight > 0.0 )
    contents.referenceDistanceCm = contents.referenceGeometry.holderHeight / PhysicalUnits::cm;

  // A reference measured in contact with the end cap has no holder at all, which
  //  is a legitimate - and common - calibration, but leaves us no distance to
  //  anchor a transfer at.  Say that, rather than letting it surface later as
  //  "non-positive reference distance".
  if( !contents.referencePoints.empty() && (contents.referenceDistanceCm <= 0.0) )
    contents.parseNotes.push_back( "The reference efficiency curve does not say what distance it"
                                   " was measured at (its geometry has no holder height), so the"
                                   " measurement distance has to be entered by hand." );

  // The holder's cap sits between the holder and whatever rests on it; ANGLE
  //  does not say whether `height` already includes it.
  if( contents.referenceGeometry.capThickness > 0.0 )
    contents.parseNotes.push_back( "The reference geometry's holder has a cap; whether its"
                                   " thickness is already included in the holder height is"
                                   " ambiguous, so check the reference distance." );

  // Fitted-region polynomials: each <region> is followed by an XML comment
  //  holding its equation.  Pair each element with its trailing comment sibling.
  const xml_node<char> * const regions_node = child( ref_node, "regions" );
  if( regions_node )
  {
    double prev_end = 0.0;
    for( const xml_node<char> *ch = regions_node->first_node(); ch; ch = ch->next_sibling() )
    {
      if( (ch->type() != rapidxml::node_element) || (SpecUtils::xml_name_str(ch) != "region") )
        continue;

      AngleOutxContents::FitRegion region;
      const Attr start = attr_num( ch, "start" );
      const Attr end = attr_num( ch, "end" );

      // `start` is only written on the first region; each subsequent region
      //  continues from where the previous one ended - so a missing `end` would
      //  silently give the NEXT region a start of zero.
      region.startKeV = start.present ? start.value : prev_end;
      if( !end.present )
      {
        contents.parseNotes.push_back( "A fitted efficiency region in the reference curve has no"
                                       " end energy; the region list was truncated there." );
        break;
      }
      region.endKeV = end.value;
      region.polynomOrder = attr_int( ch, "polynomOrder", 100 );
      prev_end = region.endKeV;

      // Grab the first comment sibling before the next element.
      for( const xml_node<char> *sib = ch->next_sibling();
           sib && (sib->type() != rapidxml::node_element); sib = sib->next_sibling() )
      {
        if( sib->type() == rapidxml::node_comment )
        {
          region.equation = SpecUtils::xml_value_str( sib );
          break;
        }
      }//for( following siblings )

      contents.fitRegions.push_back( region );
    }//for( children of <regions> )
  }//if( regions_node )

  contents.hasReference = (contents.referencePoints.size() >= 2);
}//read_reference_curve(...)


void read_energies( const xml_node<char> * const angle_node, AngleOutxContents &contents )
{
  const xml_node<char> * const energies = child( angle_node, "energies" );
  if( !energies )
    return;

  contents.energiesName = attr_str( energies, "name" );
  XML_FOREACH_CHILD( energy, energies, "energy" )
  {
    double val = 0.0;
    if( node_num( energy, val ) && (val > 0.0) )
      contents.energies.push_back( static_cast<float>(val) );
  }//XML_FOREACH_CHILD( energy, energies, "energy" )
}//read_energies(...)


void read_results( const xml_node<char> * const angle_node, AngleOutxContents &contents )
{
  const xml_node<char> * const results_node = child( angle_node, "results" );
  if( !results_node )
    return;

  XML_FOREACH_CHILD( result, results_node, "result" )
  {
    const Attr energy = attr_num( result, "energy" );
    if( !energy.present )
      continue;

    if( energy.value <= 1.0 )
      throw runtime_error( "ANGLE: <result> energy <= 1 keV (" + to_string(energy.value) + ")" );
    if( energy.value > 14000.0 )
      throw runtime_error( "ANGLE: <result> energy > 14 MeV (" + to_string(energy.value) + ")" );

    AngleResultRow row;
    row.energyKeV = energy.value;
    row.solidAngle = attr_num( result, "solidAngle" ).value;
    row.solidAnglePrecision = attr_num( result, "solidAnglePrecision" ).value;

    // `efficiency` is only written when the calculation had a reference curve.
    const Attr eff = attr_num( result, "efficiency" );
    if( eff.present )
    {
      if( (eff.value < 0.0) || IsNan(eff.value) || IsInf(eff.value) )
        throw runtime_error( "ANGLE: invalid <result> efficiency (" + to_string(eff.value)
                             + " at " + to_string(energy.value) + " keV)" );
      row.efficiency = eff.value;
      row.efficiencyPrecision = attr_num( result, "efficiencyPrecision" ).value;
      row.hasEfficiency = true;
    }//if( eff.present )

    contents.results.push_back( row );
  }//XML_FOREACH_CHILD( result, results_node, "result" )

  std::sort( begin(contents.results), end(contents.results),
             []( const AngleResultRow &l, const AngleResultRow &r ) -> bool {
               return l.energyKeV < r.energyKeV;
             } );
}//read_results(...)


void read_cascade( const xml_node<char> * const angle_node, AngleOutxContents &contents )
{
  const xml_node<char> * const cascade_node = child( angle_node, "cascadeSummingCorrections" );
  if( !cascade_node )
    return;

  XML_FOREACH_CHILD( nuc_node, cascade_node, "nuclide" )
  {
    AngleCascadeNuclide nuc;
    nuc.nuclide = attr_str( nuc_node, "name" );

    XML_FOREACH_CHILD( corr_node, nuc_node, "correction" )
    {
      AngleCascadeCorrection corr;
      corr.energyKeV = attr_num( corr_node, "energy" ).value;
      corr.value = attr_num( corr_node, "value" ).value;
      corr.branchingRatio = attr_num( corr_node, "branchingRatio" ).value;
      corr.correctedBranchingRatio = attr_num( corr_node, "correctedBranchingRatio" ).value;
      if( corr.energyKeV > 0.0 )
        nuc.corrections.push_back( corr );
    }//XML_FOREACH_CHILD( corr_node, nuc_node, "correction" )

    if( !nuc.nuclide.empty() && !nuc.corrections.empty() )
      contents.cascadeCorrections.push_back( nuc );
  }//XML_FOREACH_CHILD( nuc_node, cascade_node, "nuclide" )
}//read_cascade(...)

}//anonymous namespace


const char *to_string( const AngleDetectorType type )
{
  switch( type )
  {
    case AngleDetectorType::ClosedEndCoaxHPGe: return "Closed-end coaxial HPGe";
    case AngleDetectorType::TrueCoaxHPGe:      return "True coaxial HPGe";
    case AngleDetectorType::ClosedEndCoaxGeLi: return "Closed-end coaxial Ge(Li)";
    case AngleDetectorType::OpenEndCoaxGeLi:   return "Open-end coaxial Ge(Li)";
    case AngleDetectorType::PlanarLEPD:        return "Planar LEPD";
    case AngleDetectorType::Well:              return "Well";
    case AngleDetectorType::NaI:               return "NaI";
    case AngleDetectorType::NaIWell:           return "NaI Well";
    case AngleDetectorType::Unknown:           break;
  }//switch( type )

  return "Unknown";
}//to_string( AngleDetectorType )


AngleDetectorType angleDetectorTypeFromString( const std::string &type )
{
  const string trimmed = SpecUtils::trim_copy( type );

  // "NaI Well" must be tested before "NaI", and the coaxial spellings differ
  //  only by their prefix, so compare the whole (case-insensitive) string.
  static const AngleDetectorType all[] = {
    AngleDetectorType::NaIWell,           AngleDetectorType::NaI,
    AngleDetectorType::ClosedEndCoaxHPGe, AngleDetectorType::TrueCoaxHPGe,
    AngleDetectorType::ClosedEndCoaxGeLi, AngleDetectorType::OpenEndCoaxGeLi,
    AngleDetectorType::PlanarLEPD,        AngleDetectorType::Well
  };

  for( const AngleDetectorType t : all )
  {
    if( SpecUtils::iequals_ascii( trimmed, to_string(t) ) )
      return t;
  }

  return AngleDetectorType::Unknown;
}//angleDetectorTypeFromString(...)


const char *to_string( const AngleLayerRole role )
{
  switch( role )
  {
    case AngleLayerRole::Contact:               return "contact";
    case AngleLayerRole::ReflectingLayer:       return "reflecting layer";
    case AngleLayerRole::HousingInner:          return "inner housing";
    case AngleLayerRole::HousingOuter:          return "outer housing";
    case AngleLayerRole::VacuumInner:           return "inner vacuum gap";
    case AngleLayerRole::AntimicrophonicShield: return "antimicrophonic shield";
    case AngleLayerRole::VacuumOuter:           return "vacuum gap";
    case AngleLayerRole::EndCap:                return "end cap";
    case AngleLayerRole::EndCapWindow:          return "end cap window";
    case AngleLayerRole::EndCapCoating:         return "end cap coating";
    case AngleLayerRole::OpticalCoupling:       return "optical coupling";
    case AngleLayerRole::PmtWall:               return "photomultiplier wall";
    case AngleLayerRole::PmtWindow:             return "photomultiplier window";
    case AngleLayerRole::Other:                 break;
  }//switch( role )

  return "layer";
}//to_string( AngleLayerRole )


bool AngleMaterial::empty() const
{
  return name.empty() && elements.empty() && compounds.empty();
}


bool AngleMaterial::isBareReference() const
{
  return elements.empty() && compounds.empty();
}


AngleOutxContents AngleOutx::parse( std::istream &input, const size_t max_bytes )
{
  // rapidxml needs one mutable, null-terminated buffer.  Read directly into
  // that buffer so the untrusted XML is not duplicated in memory, and enforce
  // the limit while reading so non-seekable streams are bounded as well.
  vector<char> xml_buf;
  xml_buf.reserve( 64u * 1024u );
  std::array<char,16u * 1024u> chunk;

  while( input )
  {
    input.read( chunk.data(), static_cast<std::streamsize>(chunk.size()) );
    const std::streamsize count = input.gcount();
    if( count <= 0 )
      break;

    const size_t nread = static_cast<size_t>(count);
    if( nread > (max_bytes - xml_buf.size()) )
      throw runtime_error( "parseAngleOutxFile: input exceeds "
                           + to_string(max_bytes / (1024u*1024u)) + " MiB limit." );

    xml_buf.insert( xml_buf.end(), chunk.data(), chunk.data() + nread );
  }//while( input )

  if( input.bad() )
    throw runtime_error( "parseAngleOutxFile: failed reading input." );

  if( xml_buf.size() < 50 )
    throw runtime_error( "parseAngleOutxFile: input too small." );

  xml_buf.push_back( '\0' );

  rapidxml::xml_document<char> doc;
  try
  {
    // parse_comment_nodes: ANGLE stores each fitted-region efficiency
    //  polynomial as an XML comment after its <region> element; we read those
    //  for the developer-checks self-consistency test.
    doc.parse<rapidxml::parse_trim_whitespace | rapidxml::parse_comment_nodes>( xml_buf.data() );
  }catch( const rapidxml::parse_error &e )
  {
    throw runtime_error( "parseAngleOutxFile: XML parse error: " + string(e.what()) );
  }

  const xml_node<char> * const angle_node = child( &doc, "angle" );
  if( !angle_node )
    throw runtime_error( "parseAngleOutxFile: no <angle> root element - this is not an"
                         " ANGLE file." );

  AngleOutxContents contents;
  contents.generator = attr_str( angle_node, "generator" );
  contents.version = attr_str( angle_node, "version" );
  contents.build = attr_str( angle_node, "build" );
  contents.unitsStr = attr_str( angle_node, "units" );

  // ANGLE allows mm, cm and inches.  Guessing wrong is a silent factor-of-25.4,
  //  so only the spec's own default (an absent attribute) falls back to mm.
  if( contents.unitsStr.empty() )
  {
    contents.lengthUnitPU = PhysicalUnits::mm;
    contents.parseNotes.push_back( "The file does not say what length units it uses;"
                                   " millimeters were assumed." );
  }else if( SpecUtils::iequals_ascii( contents.unitsStr, "mm" ) )
  {
    contents.lengthUnitPU = PhysicalUnits::mm;
  }else if( SpecUtils::iequals_ascii( contents.unitsStr, "cm" ) )
  {
    contents.lengthUnitPU = PhysicalUnits::cm;
  }else if( SpecUtils::iequals_ascii( contents.unitsStr, "in" )
            || SpecUtils::iequals_ascii( contents.unitsStr, "inch" )
            || SpecUtils::iequals_ascii( contents.unitsStr, "inches" ) )
  {
    contents.lengthUnitPU = 25.4 * PhysicalUnits::mm;
  }else
  {
    throw runtime_error( "parseAngleOutxFile: unrecognized length units '"
                         + contents.unitsStr + "' (expected \"mm\", \"cm\" or \"in\")." );
  }//if( units ) / else

  const double unit = contents.lengthUnitPU;

  const xml_node<char> * const det_node = child( angle_node, "detector" );
  if( det_node )
  {
    read_detector( det_node, unit, contents );
    applyTypeSemantics( contents );
  }else
  {
    contents.modeAObstruction = "The file does not contain a detector description.";
  }

  // The top-level container/geometry/source describe the SAMPLE the <results>
  //  are for; the reference curve carries its own copies of all three.
  contents.sampleContainer = read_container( child(angle_node, "container"), unit );
  contents.sampleGeometry = read_geometry( child(angle_node, "geometry"), unit );
  contents.sampleSource = read_source( child(angle_node, "source"), unit );

  // The sample rests on top of the holder, so the holder height is the standoff
  //  from the detector face to the sample's near face.
  if( contents.sampleGeometry.holderHeight > 0.0 )
    contents.sampleDistanceCm = contents.sampleGeometry.holderHeight / PhysicalUnits::cm;

  // The source itself sits inside the container, whose foot, bottom wall and
  //  coatings all stand it further off.  (A Marinelli's source wraps around the
  //  detector rather than sitting in front of it, so its "bottom" is not a
  //  standoff; leave that case at the holder height.)
  contents.sampleSourceDistanceCm = contents.sampleDistanceCm;
  if( contents.sampleContainer.present
      && !SpecUtils::iequals_ascii( contents.sampleContainer.type, "Marinelli" ) )
  {
    double below = contents.sampleContainer.footHeight + contents.sampleContainer.bottomThickness;
    for( const AngleContainer::Coating &c : contents.sampleContainer.coatings )
      below += c.bottom;
    contents.sampleSourceDistanceCm += below / PhysicalUnits::cm;
  }//if( a non-Marinelli container )

  // Any material whose percentages had to be rescaled is worth reporting once.
  {
    const AngleMaterial *rescaled = nullptr;
    for( const AngleOutxContents::Layer &layer : contents.layers )
      if( layer.material.fractionsRescaled )
        rescaled = &layer.material;
    if( contents.sampleSource.material.fractionsRescaled )
      rescaled = &contents.sampleSource.material;
    if( rescaled )
      contents.parseNotes.push_back( "Material '" + rescaled->name + "' has mass fractions that do"
                                     " not add up to 100%; they were rescaled so they do." );
  }

  read_energies( angle_node, contents );
  node_num( child(angle_node, "precision"), contents.precision );
  node_num( child(angle_node, "elapsedTime"), contents.elapsedTimeSec );

  read_reference_curve( angle_node, unit, contents );
  read_results( angle_node, contents );
  read_cascade( angle_node, contents );

  return contents;
}//AngleOutxContents AngleOutx::parse( std::istream &input, size_t max_bytes )
