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


#pragma once

/** The scenario matrix for the near-field volumetric-efficiency validation, plus the geometry
 helpers shared by the (developer-only) Monte-Carlo truth generator and the fast committable test
 that checks InterSpec's integration against the recorded truth.

 The point of the matrix is to bracket where the flat-disk approximation stops being adequate:
   - small-near / large-near : a source at contact, narrower and then wider than the crystal.  The
                               wide one is the hard case - broad cos(theta) spread and rays leaving
                               through the side.
   - small-far / large-far   : the same sources far away, where flat-disk must agree.
   - wide-angle-far          : far on-axis but subtending a large angle, which the four corners above
                               all miss.  eta(E,theta) is not a near-field effect and does not fall
                               off with distance, so this is where a distance-only far-field test
                               would wrongly step down to flat-disk.
   - shielded                : an external shield around the source, which is the case the
                               survival-removal mu is meant to fix (plain mu_total is measured
                               -8..-16% at 60 keV there).

 Each runs with a dense and a light matrix, and includes explicit low-energy rows: 60-120 keV is
 where self-attenuation and the in-window Compton credit matter most, and where the depth-aware
 credit g(tau) is least trustworthy.
 */

#include <string>
#include <vector>

#include <Eigen/Core>

#include "io/DetectorResponse.h"
#include "efficiency/EfficiencyCalculator.h"

namespace VolNearField
{

/** Source shape.  Cylinders are end-on (InterSpec CylinderEndOn); boxes are axis-aligned with the
 detector (InterSpec Rectangular), which is CeeLo's identity rotation. */
enum class Shape
{
  Cylinder,
  Box
};//enum class Shape


/** One source geometry, in the CeeLo crystal-face frame convention (z = 0 at the crystal face,
 source in front at negative z).  Distances in cm, matching CeeLo. */
struct Scenario
{
  std::string name;
  double radius_cm = 0.0;        ///< Cylinder only: source outer radius.  Unused for a Box.
  double half_length_cm = 0.0;   ///< Half-extent ALONG THE DETECTOR AXIS: the cylinder half-length,
                                 ///<  or the box half-depth.  Shared between the shapes on purpose -
                                 ///<  it is the only dimension scenario_center() and
                                 ///<  scenario_center_distance_cm() need, so they stay shape-free.
  double standoff_cm = 0.0;      ///< endcap front to the NEAR face of the source
  bool   dense = false;          ///< dense (steel) vs light (water) source matrix
  double shield_cm = 0.0;        ///< external Fe shield around the source; 0 = bare
  // Appended after shield_cm so the positional initializers of the cylinder rows below keep
  //  working unchanged; the shapes' defaults make an un-annotated row a cylinder.
  Shape  shape = Shape::Cylinder;
  double half_width_cm = 0.0;    ///< Box only: x half-extent (transverse)
  double half_height_cm = 0.0;   ///< Box only: y half-extent (transverse)
};//struct Scenario


/** Box half-dimensions in CeeLo's (hx, hy, hz) order, hz being along the detector axis. */
inline Eigen::Vector3d scenario_box_half_dims( const Scenario &s )
{
  return Eigen::Vector3d( s.half_width_cm, s.half_height_cm, s.half_length_cm );
}


/** On-axis distance (cm) the Monte-Carlo anchor curve is recorded at.
 
 Far enough that a point source is an excellent approximation and the transfer is well conditioned,
 and close to the ANGLE file's own 25 cm reference so the MC-anchored and measured-anchored
 transfers are otherwise like-for-like.
 */
inline constexpr double kMcAnchorDistanceCm = 25.0;


/** Energies the Monte-Carlo anchor curve is sampled at: log-spaced across the validated range, with
 extra density at the low end where the efficiency curve turns over hardest. */
inline std::vector<double> anchor_energies()
{
  return { 40.0, 50.0, 60.0, 75.0, 88.0, 105.0, 122.0, 150.0, 200.0, 280.0,
           344.0, 500.0, 661.7, 900.0, 1332.5, 1800.0, 2614.0 };
}


/** Energies the truth bank is recorded at.  Weighted toward 60-120 keV deliberately - see the file
 comment - with a mid and a high line to keep the far-field limit honest. */
inline std::vector<double> scenario_energies()
{
  return { 60.0, 88.0, 122.0, 344.0, 661.7, 1332.5 };
}


/** The matrix itself.  Radii are relative to a ~2.9 cm crystal radius (the ANGLE GEM35-70 in
 test_data), so "small" is well inside it and "large" comfortably exceeds it. */
inline std::vector<Scenario> scenarios()
{
  // Named factory rather than designated initializers: the box fields sit at the end of Scenario
  //  (so the cylinder rows keep their positional braces) and spelling them positionally would mean
  //  writing out the unused cylinder radius every time.
  const auto box = []( std::string nm, const double hx, const double hy, const double hz,
                       const double standoff, const bool d, const double shield ) -> Scenario {
    Scenario s;
    s.name = std::move( nm );
    s.half_width_cm = hx;
    s.half_height_cm = hy;
    s.half_length_cm = hz;
    s.standoff_cm = standoff;
    s.dense = d;
    s.shield_cm = shield;
    s.shape = Shape::Box;
    return s;
  };

  std::vector<Scenario> v;
  for( const bool dense : { false, true } )
  {
    const char * const tag = dense ? "dense" : "light";

    v.push_back( { std::string("small-near-") + tag,      1.0,  0.5,   1.0, dense, 0.0 } );
    v.push_back( { std::string("large-near-") + tag,      4.0,  2.0,   1.0, dense, 0.0 } );
    v.push_back( { std::string("small-far-") + tag,       1.0,  0.5,  50.0, dense, 0.0 } );
    v.push_back( { std::string("large-far-") + tag,       4.0,  2.0,  50.0, dense, 0.0 } );
    // Far on-axis, but a wide flat disc still subtends ~40 degrees at its rim.
    v.push_back( { std::string("wide-angle-far-") + tag, 12.0,  0.5,  15.0, dense, 0.0 } );
    v.push_back( { std::string("shielded-near-") + tag,   2.0,  1.0,   1.0, dense, 0.5 } );

    // Boxes.  Transverse half-extent = the twin cylinder's radius, axial half-extent = its
    //  half-length, so each box row reads directly against its cylinder twin (the box holds 4/pi
    //  the volume; scenario_volume_cm3 handles that).
    //
    //  Boxes exist here because eval_rect had no per-ray aperture kernel and nothing caught it -
    //  every pre-existing row is a cylinder.  Deliberately only the CONTACT corners: a box is
    //  always a 3D integration and each element now costs ~500 ray traces, so a box row is 20-35 s
    //  against ~1 s for its cylinder twin even after the on-axis quarter-box reduction.  The
    //  far-field corners would mostly re-measure what the cylinders already cover, at that price.
    v.push_back( box( std::string("box-large-near-") + tag,      4.0,  4.0, 2.0,  1.0, dense, 0.0 ) );
    v.push_back( box( std::string("box-shielded-near-") + tag,   2.0,  2.0, 1.0,  1.0, dense, 0.5 ) );
    // Deliberately ANISOTROPIC - wide in x, narrow in y.  A cylinder cannot produce this, and it is
    //  the case where WHICH of the three exit planes wins genuinely varies with ray direction, so
    //  it is what would catch an x/y transposition or a centre-ray leftover in eval_rect.
    v.push_back( box( std::string("box-slab-near-") + tag,       4.0,  1.0, 0.5,  1.0, dense, 0.0 ) );
  }//for( dense / light )

  return v;
}//scenarios()


/** Source matrix and shield, named as InterSpec materials.
 
 Both sides of the comparison resolve these through InterSpec's MaterialDB and the MC side converts
 with CeeLoUtils::to_ceelo_material, so the model and the truth are attenuating through the SAME
 composition and density by construction - not through two hand-written definitions that could
 quietly differ.  Dense matters because that is where self-attenuation dominates and the choice of
 removal mu actually shows.
 */
inline const char *scenario_matrix_material( const bool dense )
{
  return dense ? "Stainless steel SS-304" : "Water";
}

inline const char *scenario_shield_material()
{
  return "Fe (iron)";
}


/** Volume (cm^3) of a scenario's source - the MC reports an efficiency per emitted gamma averaged
 over the source, so the model's volume integral has to be divided by this to compare. */
inline double scenario_volume_cm3( const Scenario &s )
{
  if( s.shape == Shape::Box )
    return 8.0 * s.half_width_cm * s.half_height_cm * s.half_length_cm;

  return 3.14159265358979323846 * s.radius_cm * s.radius_cm * (2.0 * s.half_length_cm);
}


/** Centre of the source cylinder in the crystal-face frame, given the endcap-front offset.
 `standoff_cm` is to the NEAR face, so the centre is one half-length further out. */
inline Eigen::Vector3d scenario_center( const Scenario &s, const double endcap_front_offset_cm )
{
  return Eigen::Vector3d( 0.0, 0.0,
                          -( endcap_front_offset_cm + s.standoff_cm + s.half_length_cm ) );
}


/** Distance from the detector's reference point to the source CENTRE, in cm - what the InterSpec
 side uses as its source-to-detector distance. */
inline double scenario_center_distance_cm( const Scenario &s )
{
  return s.standoff_cm + s.half_length_cm;
}

}//namespace VolNearField
