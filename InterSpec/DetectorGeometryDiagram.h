#ifndef DetectorGeometryDiagram_h
#define DetectorGeometryDiagram_h
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

#include <string>
#include <vector>

#include <Wt/WContainerWidget.h>

namespace ceelo{ struct GeometryDescriptor; }

/** A side-elevation cross-section (the r-z half section, mirrored about the axis) of a
 `ceelo::GeometryDescriptor`: the active crystal, its dead layer and bore, each endcap/housing
 layer, and the collimator - drawn with d3 and redrawn as the geometry form changes.

 The geometry is worked out here, in C++ (#buildModel): each region becomes a polycone profile
 in the crystal frame (front face z = 0, back face z = L, the source toward negative z), with its
 volume, mass and tooltip text.  The JavaScript only scales, draws and shows tooltips: a single
 colour (the chart theme variables), every layer at least a couple of pixels wide and everything
 else to scale, orientation picked from the space it is given.
 */
class DetectorGeometryDiagram : public Wt::WContainerWidget
{
public:
  DetectorGeometryDiagram();
  virtual ~DetectorGeometryDiagram() override;

  /** Redraws from a valid descriptor (one that passed `GeometryDescriptor::problems()`). */
  void setGeometry( const ceelo::GeometryDescriptor &gd );

  /** The form is currently invalid: keep the last drawing but dim it (class `DgdStale`). */
  void setStale( const bool stale );

  /** One plane of a polycone profile: at `z`, the region spans `rmin <= r <= rmax`.  Planes are
   listed in ascending z; a repeated z encodes a step. */
  struct Plane
  {
    double z, rmin, rmax;
  };//struct Plane

  struct Region
  {
    std::string id;     //"crystal", "dead", "bore", "layer0".., "collimator"
    std::string kind;   //"crystal" | "dead" | "void" | "layer" | "collimator" - the CSS class
    std::vector<Plane> profile;
    double volume_cm3 = 0.0, mass_g = 0.0;
    std::string tip_title;              //localized, plain text
    std::vector<std::string> tip_lines; //localized, plain text
  };//struct Region

  struct Model
  {
    bool box = false;
    double z_min = 0.0, z_max = 0.0, r_max = 0.0;
    /** Every region edge, sorted and unique - what the on-screen minimum-width rule keys on, so
     adjacent regions stay adjacent after thin gaps are widened. */
    std::vector<double> z_knots, r_knots;
    std::vector<Region> regions;   //draw order: outermost first
  };//struct Model

  /** The pure-geometry step (no Wt needed): regions, volumes, masses and tooltip text. */
  static Model buildModel( const ceelo::GeometryDescriptor &gd );

  /** The JSON the JavaScript draws from. */
  static std::string toJson( const Model &model );

protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags ) override;
  void defineJavaScript();

  /** Runs `<this element>.xs.<call>`, but only if both the element and the client-side object are
   there.  Wt sends a widget on a tab that has never been shown as a client-side *stub*, so the
   element can be absent when JavaScript members are applied - and a bare call then throws a
   TypeError that aborts the REST of Wt's update block, taking unrelated widgets' event wiring with
   it.  (DrfChart carries the same guard, for the same reason.)
   */
  void doXsJs( const std::string &method_call );

  /** The JavaScript object: `jsRef() + ".xs"`. */

  /** The latest setData(...) call made before the widget was rendered (only the newest matters). */
  /** The last geometry sent to the client, kept so #defineJavaScript can send it AGAIN whenever the
   client-side object is (re)built - Wt re-emits every JavaScript member on a full DOM recreation,
   and the object that gets built is empty.  A one-shot pending-JS queue (what this used to be) left
   the drawing permanently blank after any second full render, until the user next edited a field.
   */
  std::string m_geometryJson;

  bool m_stale;

  /** Whether the client-side object has been built.  A first render does not always carry the Full
   flag, and without this the object would never be created. */
  bool m_jsDefined;
};//class DetectorGeometryDiagram

#endif //DetectorGeometryDiagram_h
