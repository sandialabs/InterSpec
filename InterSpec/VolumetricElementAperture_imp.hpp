/* InterSpec: an application to analyze spectral gamma radiation data.

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
#ifndef VolumetricElementAperture_imp_hpp
#define VolumetricElementAperture_imp_hpp

/** ===================================================================================================
 VALIDATION REFERENCE - NOT THE PRODUCTION PATH.

 The per-ELEMENT aperture fan: for one source element, CeeLo's quadrature of ~128 rays into the
 crystal, each traced through the source and shielding by InterSpec's own ray-tracer.  This was how
 a volumetric source with a Monte-Carlo / transfer response was integrated until 2026-09-03; the
 production path is now the detector-side LINE integration (VolumetricLineIntegration_imp.hpp),
 which evaluates the same integral with the order of integration reversed, at a small fraction of
 the cost and with exact chords in the fit parameters.

 The fan is kept as the independent second quadrature of that integral, reachable only through
 `sm_volumetric_integrator_override == VolumetricIntegrator::Element` (and the tests' direct calls
 to `self_shielding_integration_imp` on a calculator with a response attached).  Two tests keep the
 two in step, and are the reason this code still exists:
   - LineVsElementScenarioMatrix (test_VolumetricLadder.cpp): the integral itself, on a quick
     subset every run and the whole scenario matrix on request;
   - LineVsElementSphericalSource (test_VolumetricLinePath.cpp): spheres, which the scenario matrix
     has none of - and which is where this reference was found to be wrong by up to a factor of two
     until eval_spherical got its own per-ray branch (2026-09-03);
   - ShellWalkMatchesElementCentreRay (test_VolumetricLinePath.cpp): the per-point shell walk
     against eval_cylinder / eval_rect / eval_spherical.
 Everything here is scalar and per element, which is exactly what made it slow (5.8e-4 s per element
 at 512 rays, 11k-67k elements per energy for a box, nothing shared across energies or evaluations).

 The `eval_*` functions in GammaInteractionCalc_imp.hpp that call this are NOT part of the reference:
 they are the adaptive volume integrator that also runs the production flat-disk model (a response
 without a Monte-Carlo characterization).  Only their `if( m_effResponse )` per-ray branches belong
 with this file; they are marked there.

 This is an implementation fragment of GammaInteractionCalc_imp.hpp and is included from it right
 where this block used to live; it must not be included on its own.
 =================================================================================================== */

/** The detector aperture as seen from one source element, expressed in the ASSEMBLY frame.

 CeeLo builds its quadrature in the crystal-face frame (z = 0 at the crystal face, source in front at
 negative z), so its ray directions have to be rotated before InterSpec's ray-tracer can follow them
 through the source and shielding.  For an axially symmetric detector the azimuth about the
 element->detector axis is a free rotation, so the minimal rotation taking CeeLo's element->detector
 direction onto InterSpec's is the whole mapping.

 `weights` and `dirs` come back parallel; `prefactor` is everything in eps_fep except the kernel, so
     eps_fep(element) == prefactor * sum_i weights[i] * t_src(dirs[i])
 with t_src == 1 for a bare element.  Scalar throughout: these depend on the element POSITION, which
 the caller supplies as a scalar - see the note on the frozen Jacobian in #eff_response_factor.
 */
struct ElementAperture
{
  std::vector<double> weights;
  std::vector<std::array<double,3>> dirs;   //unit, assembly frame, element -> detector hemisphere
  double prefactor = 0.0;
  DetectorPeakResponse::EffFlag flag = DetectorPeakResponse::EffFlag::Ok;
  double frac_sigma = 0.0;

  /** The crystal-face-frame -> assembly-frame rotation the directions were mapped through.  Kept so
   a geometry audit can ask where CeeLo's detector axis landed (test_VolumetricLadder.cpp). */
  std::array<std::array<double,3>,3> rotation = { { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} } };
};//struct ElementAperture


/** TEST HOOK - selects how #build_element_aperture defines CeeLo's element->detector direction.

 `false` (production): the direction from the element to the CeeLo-frame image of the assembly's
 detector reference point (`reference_point_position()`, i.e. the endcap front for an EndcapFront
 response) - the same physical point InterSpec's `to_det` aims at, so the two frames agree on the
 element's polar angle and the rotation reduces to the azimuthal one it should be.

 `true` (legacy, 2026-09-02 and earlier): the direction to CeeLo's ORIGIN (the crystal-face centre),
 which is `endcap_front_offset_cm` behind the point `to_det` aims at.  Off axis the two directions
 differ, and the full-frame rotation then tilts the whole fan in the plane of incidence by that
 parallax angle - several degrees for a contact element - toward more oblique source chords.  Kept
 only so the before/after can be measured (VolumetricLadder rung 0d); never set it in production.
 */
inline bool sm_aperture_frame_legacy_origin = false;


/** Builds #ElementAperture for one element.  `to_det` is the (unit) element->detector direction in
 the assembly frame, `dist` their separation in PhysicalUnits, `cos_theta` the incidence cosine. */
inline ElementAperture build_element_aperture( const ceelo::DetectorResponse &resp,
                                               const double energy_keV,
                                               const double dist,
                                               const double cos_theta,
                                               const double to_det[3],
                                               const double det_axis[3],
                                               const int n_rays )
{
  ElementAperture answer;

  // Position in the crystal-face frame: `dist` is from the DETECTOR FACE (InterSpec's one distance
  //  convention), formed through CeeLoUtils::sourcePositionFromFace exactly as the point-source
  //  query and the line cache form theirs - never through the descriptor's own reference point.
  const double theta = std::acos( std::max(-1.0, std::min(1.0, cos_theta)) );
  const Eigen::Vector3d pos = CeeLoUtils::sourcePositionFromFace( resp.descriptor, theta, 0.0,
                                                                  dist / PhysicalUnits::cm );

  const ceelo::ApertureQuadrature q = ceelo::make_aperture_quadrature( resp.geometry(), pos, n_rays );

  std::vector<Eigen::Vector3d> ceelo_dirs;
  resp.fep_ray_weights( energy_keV, q, answer.weights, ceelo_dirs );

  const ceelo::EffResult pre = resp.fep_prefactor( energy_keV, pos, q );
  answer.prefactor = pre.value;
  answer.frac_sigma = (pre.value > 0.0) ? (pre.sigma / pre.value) : 0.0;

  // CeeLo's element->detector direction, toward the SAME physical point InterSpec's `to_det` aims
  //  at: the assembly's detector reference point, whose CeeLo-frame image is
  //  CeeLoUtils::detectorFacePosition( resp.descriptor ), i.e.
  //  (0,0,-endcap_front_offset) - the detector FACE, InterSpec's one distance convention.  Aiming at
  //  CeeLo's origin instead - the legacy behaviour, see #sm_aperture_frame_legacy_origin - disagrees
  //  with `to_det` on the element's polar angle by the parallax of the endcap offset, and the frame
  //  rotation below then tilts the whole fan by that angle in the plane of incidence.
  //  Plain arithmetic rather than Eigen::Geometry - cross()/Rodrigues are a few lines, and this
  //  header is included in every translation unit that fits.
  const Eigen::Vector3d ref = sm_aperture_frame_legacy_origin ? Eigen::Vector3d( 0.0, 0.0, 0.0 )
                                                              : CeeLoUtils::detectorFacePosition( resp.descriptor );
  const Eigen::Vector3d elem_to_ref = ref - pos;
  const double pn = elem_to_ref.norm();
  double u_c[3] = { 0.0, 0.0, 1.0 };
  if( pn > 0.0 )
  {
    u_c[0] = elem_to_ref.x()/pn;  u_c[1] = elem_to_ref.y()/pn;  u_c[2] = elem_to_ref.z()/pn;
  }

#if( PERFORM_DEVELOPER_CHECKS )
  if( !sm_aperture_frame_legacy_origin && (pn > 0.0) )
  {
    // Both frames must now agree on the element's polar angle: -(u_c . a_c) == -(to_det . det_axis)
    //  == cos_theta.  a_c is (0,0,-1), so -(u_c . a_c) is just u_c[2].
    const double ceelo_cos = u_c[2];
    if( std::fabs( ceelo_cos - cos_theta ) > 1.0e-6 )
    {
      char msg[256];
      snprintf( msg, sizeof(msg), "build_element_aperture: CeeLo/assembly polar-angle mismatch"
                " (cos %.8f vs %.8f)", ceelo_cos, cos_theta );
      log_developer_error( __func__, msg );
    }
  }
#endif

  // FULL-FRAME rotation, not a minimal one-vector rotation.
  //
  //  A minimal rotation pins only u_c -> to_det and lets the residual azimuth about that axis fall
  //  out of the arithmetic.  That azimuth was justified as "free by axial symmetry", and it IS free
  //  for the detector - `weights` and `prefactor` are computed in CeeLo's own frame before any
  //  rotation, so no amount of spinning the fan changes eps_fep.  It is NOT free for the SOURCE.
  //  From an off-axis element the crystal's projection is not circular, so the fan is flattened,
  //  and the flattening plane is the plane of incidence.  Spin the fan about to_det and each ray
  //  points somewhere else through the source, which changes its chord and hence exp(-tau_i).
  //
  //  MEASURED consequence of getting this wrong (ApertureFanOrientation): four elements at one
  //  cos_theta = 0.7815, i.e. physically identical up to a rotation about the detector axis, got
  //  fan orientations of -51.8, +52.0, +74.2 and -74.0 degrees relative to the plane of incidence.
  //  A y-mirrored pair happened to map onto exact negatives, so a y-mirror preserved the chord
  //  distribution, but an x-mirrored pair did not - which made eval_rect mirror-symmetric in y to
  //  0.00% and asymmetric in x by 13-97% (BoxFoldSymmetryProbe), silently invalidating the
  //  quarter-box fold in self_shielding_integration_imp and inflating box integrals by 10-24%
  //  against an independent tensor rule (BoxOuterQuadratureVsTensorGL).
  //
  //  So map the whole frame: (element->detector, plane of incidence) in CeeLo's coordinates onto
  //  the same pair in the assembly's.  CeeLo builds `pos` at azimuth 0 with the crystal axis along
  //  z and the source in front at negative z, so its detector axis is -z under the same sign
  //  convention InterSpec uses for `det.axis` (cos_theta == -(to_det . axis) in both).
  const double a_c[3] = { 0.0, 0.0, -1.0 };

  // Orthonormal frame from a primary direction p and a secondary a: (p, in-plane, normal).
  const auto make_frame = []( const double p[3], const double a[3], double e[3][3] ) -> bool {
    e[0][0] = p[0];  e[0][1] = p[1];  e[0][2] = p[2];
    const double ap = a[0]*p[0] + a[1]*p[1] + a[2]*p[2];
    double t[3] = { a[0] - ap*p[0], a[1] - ap*p[1], a[2] - ap*p[2] };
    const double tn = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
    if( !(tn > 1.0e-9) )
      return false;   //element on the detector axis: the plane of incidence is degenerate
    for( int i = 0; i < 3; ++i )
      e[1][i] = t[i]/tn;
    e[2][0] = e[0][1]*e[1][2] - e[0][2]*e[1][1];
    e[2][1] = e[0][2]*e[1][0] - e[0][0]*e[1][2];
    e[2][2] = e[0][0]*e[1][1] - e[0][1]*e[1][0];
    return true;
  };

  double ec[3][3], ea[3][3];
  const bool have_frames = make_frame( u_c, a_c, ec ) && make_frame( to_det, det_axis, ea );

  double R[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };
  if( have_frames )
  {
    // R = sum_k  ea[k] (outer) ec[k]  - the unique rotation carrying one frame onto the other.
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        double v = 0.0;
        for( int k = 0; k < 3; ++k )
          v += ea[k][i] * ec[k][j];
        R[i][j] = v;
      }
    }
  }else
  {
    // On-axis element (either frame degenerate): the crystal projection really is circular here, so
    //  the azimuth genuinely is free and the minimal rotation is correct and sufficient.
    const double v[3] = { u_c[1]*to_det[2] - u_c[2]*to_det[1],
                          u_c[2]*to_det[0] - u_c[0]*to_det[2],
                          u_c[0]*to_det[1] - u_c[1]*to_det[0] };
    const double c = u_c[0]*to_det[0] + u_c[1]*to_det[1] + u_c[2]*to_det[2];
    const double vlen2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

    if( vlen2 > 1.0e-24 )
    {
      const double k = 1.0 / (1.0 + c);
      const double K[3][3] = { {   0.0, -v[2],  v[1] },
                               {  v[2],   0.0, -v[0] },
                               { -v[1],  v[0],   0.0 } };
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          double k2 = 0.0;
          for( int m = 0; m < 3; ++m )
            k2 += K[i][m] * K[m][j];
          R[i][j] += K[i][j] + k2*k;
        }
      }
    }else if( c < 0.0 )
    {
      R[0][0] = R[1][1] = R[2][2] = -1.0;
    }
  }//if( have_frames ) / else

  answer.dirs.reserve( ceelo_dirs.size() );
  for( const Eigen::Vector3d &d : ceelo_dirs )
  {
    answer.dirs.push_back( { R[0][0]*d.x() + R[0][1]*d.y() + R[0][2]*d.z(),
                             R[1][0]*d.x() + R[1][1]*d.y() + R[1][2]*d.z(),
                             R[2][0]*d.x() + R[2][1]*d.y() + R[2][2]*d.z() } );
  }

  for( int i = 0; i < 3; ++i )
    for( int j = 0; j < 3; ++j )
      answer.rotation[i][j] = R[i][j];

  return answer;
}//build_element_aperture(...)

#endif //VolumetricElementAperture_imp_hpp
