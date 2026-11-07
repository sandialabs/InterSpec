#ifndef DrfAppearance_h
#define DrfAppearance_h
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
#include <memory>
#include <string>

namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml


/** Optional pictures of the instrument a DRF is for: a side-view SVG icon and an exterior 3D model
 (glTF 2.0).  Purely cosmetic - neither is used in any calculation, so neither is part of the DRF
 hash.  A DRF holds one as `std::shared_ptr<const DrfAppearance>`; see
 `DetectorPeakResponse::appearance()`.

 Both files share one frame, the DRF's distance reference: the origin is the point on the detector
 axis that is placed at the user's source-to-detector distance (the "detector face").  For a DRF
 with a CeeLo geometry that is the endcap front, `GeometryDescriptor::endcap_front_offset_cm()` ahead
 of the crystal face, so the DRF's crystal sits inside the model with no offset.  A model therefore
 belongs to a DRF's distance reference rather than to the hardware: a DRF of the same instrument
 that measures distance to another surface needs a copy with its origin moved.  Both files show
 only the outside; the crystal and its housing are drawn from the DRF's geometry.

 - Icon (SVG): a side view.  SVG user units are millimetres; (0,0) is the origin above; the face is
   on the left (the instrument looks along -x); its top is up (-y, as SVG's y points down).  The
   viewBox is the drawn extent and `width`/`height` are ignored, so the same file is drawn as a
   16 px icon or to scale.  Black on transparent - the app tints the silhouette, so colours other
   than opacity are lost.  Only ever shown as an image (never inlined), so it can not run scripts.
 - Exterior model (glTF 2.0 JSON, buffers embedded as data URIs): metres, +Y up, the instrument
   looks along +Z (glTF's "front"), so its body is at z < 0.  Triangle meshes with POSITION
   (optionally NORMAL), node transforms, and `baseColorFactor` are used; textures are ignored.
   Drawn by `InterSpec_resources/GltfMini.js`.

 Serialized as the optional `<Appearance>` element of a DRF's XML (and of its database extras):
 \code{.xml}
 <Appearance version="0">
   <Icon format="svg">...escaped SVG...</Icon>
   <ExteriorModel format="gltf">...glTF JSON...</ExteriorModel>
 </Appearance>
 \endcode
 */
struct DrfAppearance
{
  /** The icon, as returned by #validated_svg; empty if none. */
  std::string icon_svg;

  /** The exterior model, as returned by #validated_gltf; empty if none. */
  std::string model_gltf;

  bool empty() const;
  bool operator==( const DrfAppearance &rhs ) const;

  /** The icon as a `data:image/svg+xml;base64,...` URI; empty if there is no icon. */
  std::string icon_data_uri() const;

  /** Markup showing the icon as a silhouette in the surrounding text colour (so in any theme): an
   inline-block `<span class="DrfIcon">` (see InterSpec.css) masked by the icon, `height_px` tall and
   as wide as the viewBox gives, up to `max_width_px`.  Built only from numbers and base64, so it is
   safe as `Wt::TextFormat::UnsafeXHTML` - which it needs, since Wt's XSS filter removes `data:` URLs.
   Empty if there is no icon.
   */
  std::string icon_html( const double height_px, const double max_width_px ) const;

  /** The icon's viewBox, {x, y, width, height}, in mm; throws if there is no icon. */
  std::array<double,4> icon_view_box() const;

  /** The model's axis-aligned bounds {xmin, ymin, zmin, xmax, ymax, zmax} in metres, from its
   POSITION accessors' min/max.  Returns false if there is no model, or its bounds can not be found.
   */
  bool model_bounds_m( std::array<double,6> &bounds ) const;

  /** Checks an SVG for use as an icon and returns the form to store: re-written from its parsed
   elements, without comments or the XML declaration.  Throws std::exception, with a message for the
   user, if it has no usable viewBox, is over 64 kB, or has any script, event handler, foreign
   object, or link outside itself.
   */
  static std::string validated_svg( const std::string &svg );

  /** Checks a glTF 2.0 model (`.gltf` JSON, or a binary `.glb`, which is converted) and returns
   compact JSON with its buffers embedded.  Throws std::exception, with a message for the user, if
   it needs external files or an extension, has no triangle mesh, is over 2 MB, or is inconsistent.
   */
  static std::string validated_gltf( const std::string &gltf_or_glb );

  /** Appends an `<Appearance>` element to `parent` (nothing if #empty). */
  void toXml( ::rapidxml::xml_node<char> *parent, ::rapidxml::xml_document<char> *doc ) const;

  /** Reads an `<Appearance>` element; returns nullptr if it holds nothing usable.  The icon is
   re-checked with #validated_svg (an icon that fails is dropped); the model is taken as is.
   */
  static std::shared_ptr<const DrfAppearance> fromXml( const ::rapidxml::xml_node<char> *node );
};//struct DrfAppearance

#endif //DrfAppearance_h
