#ifndef VolumetricReferenceIntegrator_h
#define VolumetricReferenceIntegrator_h
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

/** INDEPENDENT QUADRATURE REFERENCES for the volumetric-source efficiency integral - test-only.

 The line path (VolumetricLineIntegration_imp.hpp) and the element path (eval_* in
 GammaInteractionCalc_imp.hpp) integrate the same integrand:

     eps = Int_V dV rho(r) P(r) (1/4pi) Int dOmega k(line) T(r, w)

 with rho the emission density (per the TraceActivityType contract), P the response's prefactor
 (ceelo::DetectorResponse::fep_prefactor), k the per-line full-energy interaction probability in the
 crystal (fep_line_probabilities), T the transmission through the source and shields.  Every
 difference between the two production paths is quadrature - PLUS one approximation each: the line
 path interpolates P from a grid, the element path discretises the crystal kernel with a frozen
 128-ray fan.  Monte Carlo (CeeLo) cannot arbitrate quadrature: a model-vs-MC difference contains the
 transfer's own physics error.  So here are two more integrations of the SAME integrand, sharing
 NOTHING with either production quadrature but the CeeLo response object, the shell coefficients and
 the emission contract:

   (A) `reference_tensor_gl` - a PER-VOXEL "EFFTRAN-style" calculation: Gauss-Legendre voxels in
       source-adapted coordinates, each voxel evaluated with the response's own point kernel
       (`eps_fep_element`: P times a Fibonacci aperture fan with per-ray transmission), summed with
       the emission density.  Deterministic; refine the node counts until it stops moving.  Its only
       approximation is the per-voxel fan, which the caller sweeps.
   (B) `reference_random` - plain Monte Carlo over BOTH the emission point and the emission
       direction, ONE ray per sample, with a standard error from independent streams: no fan, no
       grid, no low-discrepancy sequence, no analytic chord.  Importance sampling (1/d^2 about the
       detector's foot point) keeps a 20 m disk tractable.

 Both evaluate the emission density and the normalisation from the contract at
 GammaInteractionCalc::TraceActivityType, written here independently of the production code; the
 shell walk is the line-side `shell_path_to_point_imp`, which `ShellWalkMatchesElementCentreRay`
 pins against the element walkers to 1e-10 (a brute-force ray march would be the fully independent
 alternative; `PerRayKernelIdentityVsRayMarch` already does that check on the kernel).

 Units: the calculator is in PhysicalUnits, CeeLo in cm; every conversion is explicit below.
 The result has the same meaning and units as `DistributedSrcCalcT::integral` on either path.
 */

#include <array>
#include <cmath>
#include <mutex>
#include <random>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <functional>

#include <Eigen/Dense>

#include "SpecUtils/SpecUtilsAsync.h"

#include "io/DetectorResponse.h"
#include "io/ResponseKernel.h"
#include "geometry/Geometry.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/GammaInteractionCalc_imp.hpp"

namespace VolRef
{
using GammaInteractionCalc::ShellType;
using GammaInteractionCalc::GeometryType;
using GammaInteractionCalc::DistributedSrcCalcT;

/** Gauss-Legendre nodes and weights on [-1,1] (Newton iteration on the Legendre polynomial). */
inline void gauss_legendre( const int n, std::vector<double> &x, std::vector<double> &w )
{
  x.assign( n, 0.0 );
  w.assign( n, 0.0 );
  const double pi = PhysicalUnits::pi;
  for( int i = 0; i < (n + 1)/2; ++i )
  {
    double z = std::cos( pi*(i + 0.75)/(n + 0.5) );
    double pp = 0.0;
    for( int iter = 0; iter < 100; ++iter )
    {
      double p1 = 1.0, p2 = 0.0;
      for( int j = 0; j < n; ++j )
      {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0*j + 1.0)*z*p2 - j*p3)/(j + 1.0);
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      const double z1 = z;
      z = z1 - p1/pp;
      if( std::fabs(z - z1) < 1.0e-15 )
        break;
    }
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = 2.0/((1.0 - z*z)*pp*pp);
    w[n - 1 - i] = w[i];
  }
}//gauss_legendre(...)


/** The frame, response and emission contract of one calculator, as the references need them.
 `calc.finalize_shell_coefficients()` must already have been called (the harness's
 `build_scenario_calc` leaves the FEP coefficients at a negative sentinel until then). */
struct ReferenceModel
{
  const DistributedSrcCalcT<double> &calc;
  const ceelo::DetectorResponse &resp;
  const ceelo::Geometry &geom;
  double M[3][3];                //crystal -> assembly rotation (GammaInteractionCalc::detector_frame_rotation)
  Eigen::Vector3d r0;            //crystal-frame position (cm) of the assembly's detector point (the face)
  Eigen::Vector3d crystal_centre_a;   //assembly frame (PU)
  double crystal_bound_radius = 0.0;  //PU: sphere about crystal_centre_a containing the active crystal
  double walk_origin_dist = 0.0;      //PU: how far behind a point the walk's detector-side origin is put
  double emission_scale = 1.0;        //1/V, A/norm, or 1 (see emission())

  explicit ReferenceModel( const DistributedSrcCalcT<double> &c )
    : calc( c ), resp( *c.m_effResponse ), geom( c.m_effResponse->geometry() )
  {
    if( !c.m_effResponse )
      throw std::runtime_error( "ReferenceModel: the calculator has no CeeLo response" );
    if( resp.descriptor.collimator )
      throw std::runtime_error( "ReferenceModel: collimated responses are not supported" );
    for( size_t l = 0; l < c.m_shells.size(); ++l )
      if( c.m_shells[l].fep_trans_len_coef < 0.0 )
        throw std::runtime_error( "ReferenceModel: call finalize_shell_coefficients() first" );

    GammaInteractionCalc::detector_frame_rotation( c.m_detector.axis, 0.0, M );
    r0 = CeeLoUtils::detectorFacePosition( resp.descriptor );

    const double cm = PhysicalUnits::cm;
    const double L = geom.detector_length();
    const Eigen::Vector3d centre_c( 0.0, 0.0, 0.5*L );
    double centre_a[3];
    dir_to_assembly( (centre_c - r0)*cm, centre_a );
    for( int i = 0; i < 3; ++i )
      crystal_centre_a[i] = c.m_detector.position[i] + centre_a[i];
    const double transverse = (geom.shape() == ceelo::DetectorShape::Cylinder)
                                ? geom.detector_radius()
                                : std::hypot( geom.detector_half_x(), geom.detector_half_y() );
    crystal_bound_radius = std::hypot( transverse, 0.5*L ) * cm;

    double det_dist = 0.0;
    for( int i = 0; i < 3; ++i )
      det_dist += c.m_detector.position[i]*c.m_detector.position[i];
    det_dist = std::sqrt( det_dist );
    const std::array<double,3> &od = c.m_shells.back().dims;
    walk_origin_dist = 4.0*(det_dist + od[0] + od[1] + od[2]) + 100.0*cm;

    emission_scale = compute_emission_scale();
  }

  /** Crystal-frame position (cm) of an assembly-frame point (PU). */
  Eigen::Vector3d to_crystal( const double p[3] ) const
  {
    const double cm = PhysicalUnits::cm;
    Eigen::Vector3d pc;
    for( int i = 0; i < 3; ++i )
    {
      double v = 0.0;
      for( int j = 0; j < 3; ++j )
        v += M[j][i]*(p[j] - calc.m_detector.position[j]);   //M^T (p - det)
      pc[i] = v/cm + r0[i];
    }
    return pc;
  }

  Eigen::Vector3d dir_to_crystal( const double d[3] ) const
  {
    Eigen::Vector3d dc;
    for( int i = 0; i < 3; ++i )
      dc[i] = M[0][i]*d[0] + M[1][i]*d[1] + M[2][i]*d[2];
    return dc;
  }

  void dir_to_assembly( const Eigen::Vector3d &dc, double d[3] ) const
  {
    for( int i = 0; i < 3; ++i )
      d[i] = M[i][0]*dc.x() + M[i][1]*dc.y() + M[i][2]*dc.z();
  }

  GeometryType geometry() const { return calc.m_geometry; }
  const std::array<double,3> &outer_dims() const { return calc.m_shells[calc.m_materialIndex].dims; }
  bool hollow() const { return calc.m_materialIndex > 0; }
  const std::array<double,3> &inner_dims() const { return calc.m_shells[calc.m_materialIndex - 1].dims; }

  static bool inside_solid( const GeometryType g, const std::array<double,3> &d, const double p[3] )
  {
    switch( g )
    {
      case GeometryType::Spherical:
        return (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) < d[0]*d[0];
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        return ((p[0]*p[0] + p[1]*p[1]) < d[0]*d[0]) && (std::fabs(p[2]) < d[1]);
      case GeometryType::Rectangular:
        return (std::fabs(p[0]) < d[0]) && (std::fabs(p[1]) < d[1]) && (std::fabs(p[2]) < d[2]);
      case GeometryType::NumGeometryType:
        break;
    }
    return false;
  }

  static double solid_volume( const GeometryType g, const std::array<double,3> &d )
  {
    const double pi = PhysicalUnits::pi;
    switch( g )
    {
      case GeometryType::Spherical:      return (4.0/3.0)*pi*d[0]*d[0]*d[0];
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn: return 2.0*pi*d[0]*d[0]*d[1];
      case GeometryType::Rectangular:    return 8.0*d[0]*d[1]*d[2];
      case GeometryType::NumGeometryType: break;
    }
    return 0.0;
  }

  /** Inside the SOURCE shell: the outer solid, minus the core of a hollow source. */
  bool in_shell( const double p[3] ) const
  {
    if( !inside_solid( geometry(), outer_dims(), p ) )
      return false;
    if( hollow() && inside_solid( geometry(), inner_dims(), p ) )
      return false;
    return true;
  }

  /** Depth below the emitting face, per the in-situ contract (the face toward +z, or the outer
   radius, is the emitting surface). */
  double depth( const double p[3] ) const
  {
    const std::array<double,3> &d = outer_dims();
    switch( geometry() )
    {
      case GeometryType::Spherical:      return d[0] - std::sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
      case GeometryType::CylinderEndOn:  return d[1] - p[2];
      case GeometryType::CylinderSideOn: return d[0] - std::sqrt( p[0]*p[0] + p[1]*p[1] );
      case GeometryType::Rectangular:    return d[2] - p[2];
      case GeometryType::NumGeometryType: break;
    }
    return 0.0;
  }

  /** Emitting-surface area and depth-integral normalisation of the in-situ contract, from the
   derivation (scratch DERIVATION.md section 6) rather than from the production helpers. */
  double compute_emission_scale() const
  {
    const double pi = PhysicalUnits::pi;
    const std::array<double,3> &d = outer_dims();
    if( !calc.m_normalizeByVolume )
      return 1.0;

    if( calc.m_isInSituExponential )
    {
      const double L = calc.m_inSituRelaxationLength;
      if( hollow() )
        throw std::runtime_error( "ReferenceModel: an in-situ source must be the innermost layer" );
      // norm = Int_V exp(-depth/L) dV, computed as a 1D integral by 200-point Gauss-Legendre in the
      //  depth coordinate (exact to round-off for the smooth exponential x polynomial integrands).
      std::vector<double> x, w;
      gauss_legendre( 200, x, w );
      double norm = 0.0, area = 0.0;
      switch( geometry() )
      {
        case GeometryType::Spherical:
          area = 4.0*pi*d[0]*d[0];
          for( size_t i = 0; i < x.size(); ++i )
          {
            const double r = 0.5*d[0]*(x[i] + 1.0);
            norm += 0.5*d[0]*w[i] * 4.0*pi*r*r*std::exp( -(d[0] - r)/L );
          }
          break;
        case GeometryType::CylinderSideOn:
          area = 2.0*pi*d[0]*(2.0*d[1]);
          for( size_t i = 0; i < x.size(); ++i )
          {
            const double r = 0.5*d[0]*(x[i] + 1.0);
            norm += 0.5*d[0]*w[i] * 2.0*pi*r*(2.0*d[1])*std::exp( -(d[0] - r)/L );
          }
          break;
        case GeometryType::CylinderEndOn:
          area = pi*d[0]*d[0];
          for( size_t i = 0; i < x.size(); ++i )
          {
            const double z = d[1]*(x[i] + 1.0);   //depth in [0, 2H]
            norm += d[1]*w[i] * area * std::exp( -z/L );
          }
          break;
        case GeometryType::Rectangular:
          area = (2.0*d[0])*(2.0*d[1]);
          for( size_t i = 0; i < x.size(); ++i )
          {
            const double z = d[2]*(x[i] + 1.0);
            norm += d[2]*w[i] * area * std::exp( -z/L );
          }
          break;
        case GeometryType::NumGeometryType:
          break;
      }
      return area / norm;
    }//in-situ

    double vol = solid_volume( geometry(), d );
    if( hollow() )
      vol -= solid_volume( geometry(), inner_dims() );
    return 1.0 / vol;
  }//compute_emission_scale()

  /** Emission density at p: zero outside the source shell. */
  double emission( const double p[3] ) const
  {
    if( !in_shell( p ) )
      return 0.0;
    if( calc.m_isInSituExponential )
      return emission_scale * std::exp( -depth( p )/calc.m_inSituRelaxationLength );
    return emission_scale;
  }

  /** exp(-tau) from p along the unit assembly direction u to the detector (through every shell's
   FEP removal coefficient, then air from the outermost shell's exit to the endcap plane). */
  double transmission( const double p[3], const double u[3] ) const
  {
    const double S = walk_origin_dist;
    const double o[3] = { p[0] + S*u[0], p[1] + S*u[1], p[2] + S*u[2] };
    GammaInteractionCalc::ShellPathT<double> path;
    GammaInteractionCalc::shell_path_to_point_imp<double>( geometry(), calc.m_shells, o, u, S, path );

    double tau = 0.0;
    for( size_t l = 0; l < calc.m_shells.size(); ++l )
    {
      if( !path.crossed[l] )
        continue;
      const auto &shell = calc.m_shells[l];
      tau += (shell.type == ShellType::Material) ? shell.fep_trans_len_coef * path.own_len[l]
                                                 : shell.fep_trans_len_coef;
    }

    if( calc.m_attenuateForAir && (calc.m_airTransLenCoef > 0.0) )
    {
      // Air from the outermost shell's exit to the endcap-face plane through the detector point.
      const double *det = calc.m_detector.position;
      const double *axis = calc.m_detector.axis;
      // The axis points from the detector INTO the assembly (detector_geom_from_config), so a photon
      //  heading to the detector has u.axis < 0.
      const double u_ax = u[0]*axis[0] + u[1]*axis[1] + u[2]*axis[2];
      if( u_ax < 0.0 )
      {
        const double t_face = ((det[0] - p[0])*axis[0] + (det[1] - p[1])*axis[1] + (det[2] - p[2])*axis[2]) / u_ax;
        const double t_exit = S - path.air;
        tau += calc.m_airTransLenCoef * std::max( 0.0, t_face - t_exit );
      }
    }
    return std::exp( -tau );
  }

  double prefactor( const Eigen::Vector3d &pc ) const
  {
    static const ceelo::ApertureQuadrature no_quadrature;
    return resp.fep_prefactor( calc.m_energy, pc, no_quadrature ).value;
  }

  /** The per-voxel EFFTRAN-style kernel: the response's own point efficiency at p with the
   per-ray source/shield transmission folded in (P times the aperture fan). */
  double point_efficiency( const double p[3], const int n_rays ) const
  {
    const Eigen::Vector3d pc = to_crystal( p );
    const ceelo::ApertureQuadrature q = ceelo::make_aperture_quadrature( geom, pc, n_rays );
    const std::function<double(const Eigen::Vector3d &)> t_src = [this,p]( const Eigen::Vector3d &dir_c ){
      double u[3];
      dir_to_assembly( dir_c, u );
      return transmission( p, u );
    };
    return resp.eps_fep_element( calc.m_energy, pc, q, t_src ).value;
  }

  /** A direction cone from p that contains the whole active crystal (its bounding sphere): the
   sampling measure of the random reference.  Full sphere when p is inside the bound. */
  void direction_cone( const double p[3], double axis[3], double &cos_alpha ) const
  {
    double d = 0.0;
    for( int i = 0; i < 3; ++i )
    {
      axis[i] = crystal_centre_a[i] - p[i];
      d += axis[i]*axis[i];
    }
    d = std::sqrt( d );
    if( d <= crystal_bound_radius )
    {
      axis[0] = axis[1] = 0.0;
      axis[2] = 1.0;
      cos_alpha = -1.0;
      return;
    }
    for( int i = 0; i < 3; ++i )
      axis[i] /= d;
    const double s = crystal_bound_radius/d;
    cos_alpha = std::sqrt( std::max( 0.0, 1.0 - s*s ) );
  }
};//struct ReferenceModel


// ---------------------------------------------------------------------------------------------
// (A) Deterministic per-voxel reference
// ---------------------------------------------------------------------------------------------

struct GlRefResult
{
  double value = 0.0;
  size_t num_evals = 0;
};

/** Options of the tensor Gauss-Legendre voxel reference.  `n_a/n_b/n_c` are the node counts on
 the three source coordinates - (r, cos, phi) for a sphere, (rho, phi, z) for cylinders, (x, y, z)
 for a box.  On-axis symmetry collapses the azimuth to one node (2 pi) / a quarter box.
 `log_radial` maps a cylinder's rho through u = ln(1 + rho^2/h^2) about the axis (h = the standoff
 to the emitting face), which makes a wide disk's 1/d^2 integrand smooth in u.  `exp_depth` maps an
 in-situ source's depth through y = exp(-depth/L), which makes the profile a constant factor. */
struct GlRefOptions
{
  int n_a = 32, n_b = 32, n_c = 32;
  int n_rays = 1024;
  bool log_radial = true;
  bool exp_depth = true;
  bool multithread = true;
};

/** One integration sub-domain: a box in the source coordinates with the emission indicator
 handling the rest (hollow sources are tiled so no voxel straddles the core). */
struct GlDomain
{
  double lo[3], hi[3];
};

inline std::vector<GlDomain> gl_domains( const ReferenceModel &m )
{
  const std::array<double,3> &o = m.outer_dims();
  std::vector<GlDomain> out;
  const auto add = [&out]( const double l0, const double h0, const double l1, const double h1,
                           const double l2, const double h2 ){
    if( (h0 <= l0) || (h1 <= l1) || (h2 <= l2) )
      return;
    out.push_back( GlDomain{ { l0, l1, l2 }, { h0, h1, h2 } } );
  };
  const double pi = PhysicalUnits::pi;
  switch( m.geometry() )
  {
    case GeometryType::Spherical:
    {
      const double ri = m.hollow() ? m.inner_dims()[0] : 0.0;
      add( ri, o[0], -1.0, 1.0, 0.0, 2.0*pi );
      break;
    }
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      if( m.hollow() )
      {
        const std::array<double,3> &in = m.inner_dims();
        add( 0.0, o[0], 0.0, 2.0*pi, in[1], o[1] );
        add( 0.0, o[0], 0.0, 2.0*pi, -o[1], -in[1] );
        add( in[0], o[0], 0.0, 2.0*pi, -in[1], in[1] );
      }else
      {
        add( 0.0, o[0], 0.0, 2.0*pi, -o[1], o[1] );
      }
      break;
    }
    case GeometryType::Rectangular:
    {
      if( m.hollow() )
      {
        const std::array<double,3> &in = m.inner_dims();
        add( -o[0], o[0], -o[1], o[1], in[2], o[2] );
        add( -o[0], o[0], -o[1], o[1], -o[2], -in[2] );
        add( -o[0], o[0], in[1], o[1], -in[2], in[2] );
        add( -o[0], o[0], -o[1], -in[1], -in[2], in[2] );
        add( in[0], o[0], -in[1], in[1], -in[2], in[2] );
        add( -o[0], -in[0], -in[1], in[1], -in[2], in[2] );
      }else
      {
        add( -o[0], o[0], -o[1], o[1], -o[2], o[2] );
      }
      break;
    }
    case GeometryType::NumGeometryType:
      break;
  }
  return out;
}//gl_domains(...)


inline GlRefResult reference_tensor_gl( const DistributedSrcCalcT<double> &calc, const GlRefOptions &opt )
{
  const ReferenceModel m( calc );
  const GeometryType g = m.geometry();
  const bool on_axis = (calc.m_detector.position[0] == 0.0) && (calc.m_detector.position[1] == 0.0);
  const bool cylinder = (g == GeometryType::CylinderEndOn) || (g == GeometryType::CylinderSideOn);
  const double pi = PhysicalUnits::pi;
  const double L = calc.m_inSituRelaxationLength;
  const bool in_situ = calc.m_isInSituExponential;

  // Axial standoff of the detector point above the emitting face, for the log-radial map.
  const double h_face = calc.m_detector.position[2] - ((g == GeometryType::Rectangular) ? m.outer_dims()[2]
                                                                                        : m.outer_dims()[1]);
  const bool use_log_radial = opt.log_radial && cylinder && on_axis && (h_face > 0.0);
  const bool use_exp_depth = opt.exp_depth && in_situ
                             && ((g == GeometryType::CylinderEndOn) || (g == GeometryType::Rectangular));

  std::vector<double> xa, wa, xb, wb, xc, wc;
  gauss_legendre( opt.n_a, xa, wa );
  gauss_legendre( opt.n_b, xb, wb );
  gauss_legendre( opt.n_c, xc, wc );

  const std::vector<GlDomain> domains = gl_domains( m );
  GlRefResult result;

  for( const GlDomain &dom : domains )
  {
    // Coordinate ranges and per-axis node counts after the symmetry reductions and substitutions.
    double lo[3] = { dom.lo[0], dom.lo[1], dom.lo[2] };
    double hi[3] = { dom.hi[0], dom.hi[1], dom.hi[2] };
    int n[3] = { opt.n_a, opt.n_b, opt.n_c };
    const std::vector<double> *nx[3] = { &xa, &xb, &xc };
    const std::vector<double> *nw[3] = { &wa, &wb, &wc };
    double sym_factor = 1.0;

    // Azimuthal symmetry: one node in phi, weight 2 pi.
    int phi_axis = -1;
    if( on_axis && (g == GeometryType::Spherical) )
      phi_axis = 2;
    if( on_axis && cylinder )
      phi_axis = 1;
    if( phi_axis >= 0 )
    {
      n[phi_axis] = 1;
      sym_factor *= 2.0*pi;
    }
    // Quarter box.
    if( on_axis && (g == GeometryType::Rectangular) && !m.hollow() )
    {
      lo[0] = 0.0;
      lo[1] = 0.0;
      sym_factor *= 4.0;
    }

    // The cylinder's radial coordinate: rho, or u = ln(1 + rho^2/h^2).
    int radial_axis = -1;
    if( cylinder )
      radial_axis = 0;
    // Depth substitution: z -> y = exp(-depth/L) on the axial coordinate of end-on/rect sources.
    const int depth_axis = (g == GeometryType::Rectangular) ? 2 : (cylinder ? 2 : -1);

    for( int i0 = 0; i0 < n[0]; ++i0 )
    {
      // Threads over the outermost coordinate would be natural; the per-voxel fan is the cost, so
      //  thread over voxels of this i0 slab instead (fixed-order reduction per slab).
      const size_t slab = static_cast<size_t>(n[1]) * static_cast<size_t>(n[2]);
      std::vector<double> vals( slab, 0.0 );

      const auto eval_voxel = [&]( const size_t idx ){
        const int i1 = static_cast<int>( idx / n[2] );
        const int i2 = static_cast<int>( idx % n[2] );
        const int ii[3] = { i0, i1, i2 };
        double c[3], jac = sym_factor;
        for( int ax = 0; ax < 3; ++ax )
        {
          if( n[ax] == 1 )
          {
            c[ax] = (ax == phi_axis) ? 0.0 : 0.5*(lo[ax] + hi[ax]);
            continue;
          }
          const double t = (*nx[ax])[ii[ax]];
          const double wt = (*nw[ax])[ii[ax]];
          if( (ax == radial_axis) && use_log_radial )
          {
            // rho = h sqrt(e^u - 1), rho drho = (h^2/2) e^u du
            const double u_lo = std::log( 1.0 + lo[ax]*lo[ax]/(h_face*h_face) );
            const double u_hi = std::log( 1.0 + hi[ax]*hi[ax]/(h_face*h_face) );
            const double u = 0.5*(u_hi - u_lo)*t + 0.5*(u_hi + u_lo);
            c[ax] = h_face*std::sqrt( std::max( 0.0, std::exp(u) - 1.0 ) );
            jac *= wt * 0.5*(u_hi - u_lo) * 0.5*h_face*h_face*std::exp( u );   //includes the rho of rho drho
          }else if( (ax == depth_axis) && use_exp_depth )
          {
            // depth = top - z in [d_lo, d_hi]; y = exp(-depth/L): dz = L dy / y, and the profile
            //  exp(-depth/L) = y, so the profile times dz is L dy - the emission profile is divided
            //  back out of the integrand below (it multiplies through emission()).
            const double top = m.outer_dims()[(g == GeometryType::Rectangular) ? 2 : 1];   //the emitting face
            const double d_lo = top - hi[ax], d_hi = top - lo[ax];
            const double y_lo = std::exp( -d_hi/L ), y_hi = std::exp( -d_lo/L );
            const double y = 0.5*(y_hi - y_lo)*t + 0.5*(y_hi + y_lo);
            const double depth = -L*std::log( y );
            c[ax] = top - depth;
            jac *= wt * 0.5*(y_hi - y_lo) * L / y;
          }else
          {
            c[ax] = 0.5*(hi[ax] - lo[ax])*t + 0.5*(hi[ax] + lo[ax]);
            jac *= wt * 0.5*(hi[ax] - lo[ax]);
            if( ax == radial_axis )
              jac *= c[ax];   //rho drho
          }
        }//for( axes )

        double p[3];
        switch( g )
        {
          case GeometryType::Spherical:
          {
            const double r = c[0], ct = c[1], ph = c[2];
            const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
            p[0] = r*st*std::cos( ph );
            p[1] = r*st*std::sin( ph );
            p[2] = r*ct;
            jac *= r*r;
            break;
          }
          case GeometryType::CylinderEndOn:
          case GeometryType::CylinderSideOn:
            p[0] = c[0]*std::cos( c[1] );
            p[1] = c[0]*std::sin( c[1] );
            p[2] = c[2];
            break;
          case GeometryType::Rectangular:
            p[0] = c[0];
            p[1] = c[1];
            p[2] = c[2];
            break;
          case GeometryType::NumGeometryType:
            break;
        }

        const double rho = m.emission( p );
        if( !(rho > 0.0) )
          return;
        vals[idx] = jac * rho * m.point_efficiency( p, opt.n_rays );
      };//eval_voxel

      if( opt.multithread && (slab > 8) )
      {
        SpecUtilsAsync::ThreadPool pool;
        const size_t nchunk = 32;
        for( size_t ch = 0; ch < nchunk; ++ch )
        {
          pool.post( [&,ch](){
            for( size_t idx = (slab*ch)/nchunk; idx < (slab*(ch + 1))/nchunk; ++idx )
              eval_voxel( idx );
          } );
        }
        pool.join();
      }else
      {
        for( size_t idx = 0; idx < slab; ++idx )
          eval_voxel( idx );
      }

      for( const double v : vals )
        result.value += v;
      result.num_evals += slab;
    }//for( i0 )
  }//for( domains )

  return result;
}//reference_tensor_gl(...)


// ---------------------------------------------------------------------------------------------
// (B) Random reference: one ray per sample
// ---------------------------------------------------------------------------------------------

struct McRefResult
{
  double value = 0.0;
  double std_err = 0.0;     //from the per-sample variance (independent draws)
  size_t n = 0;
  double ess = 0.0;         //effective sample size (sum w)^2 / sum w^2
};

enum class RefSampling
{
  /** Emission point uniform in the outer solid (zero weight in a hollow core). */
  Uniform,
  /** End-on cylinders and boxes: transverse density about the detector's foot point on the emitting
   face proportional to 1/(rho^2 + h^2), depth from the in-situ profile when there is one - the
   1/d^2 of the flux cancels in the weight.  Falls back to Uniform elsewhere. */
  InverseSquare
};

struct McRefOptions
{
  size_t n_samples = 1 << 20;
  uint64_t seed = 1;
  RefSampling sampling = RefSampling::InverseSquare;
  bool multithread = true;
  size_t num_chunks = 64;   //independent streams; fixed so the result is thread-count independent
};

inline McRefResult reference_random( const DistributedSrcCalcT<double> &calc, const McRefOptions &opt )
{
  const ReferenceModel m( calc );
  const GeometryType g = m.geometry();
  const std::array<double,3> od = m.outer_dims();
  const double pi = PhysicalUnits::pi;
  const bool in_situ = calc.m_isInSituExponential;
  const double L = calc.m_inSituRelaxationLength;

  // Inverse-square transverse sampling is available for the two geometries whose emitting face is
  //  the +z face and the detector's foot point lies inside that face.
  const double *det = calc.m_detector.position;
  const bool facing = (g == GeometryType::CylinderEndOn) || (g == GeometryType::Rectangular);
  const double top = facing ? od[(g == GeometryType::Rectangular) ? 2 : 1] : 0.0;
  const double h0 = det[2] - top;
  bool foot_inside = false;
  if( facing )
  {
    if( g == GeometryType::CylinderEndOn )
      foot_inside = (std::hypot( det[0], det[1] ) < 0.999*od[0]);
    else
      foot_inside = (std::fabs(det[0]) < 0.999*od[0]) && (std::fabs(det[1]) < 0.999*od[1]);
  }
  const bool inv_sq = (opt.sampling == RefSampling::InverseSquare) && facing && foot_inside && (h0 > 0.0);
  const double vol_outer = ReferenceModel::solid_volume( g, od );

  // Distance from the foot point to the face boundary along the in-plane direction (cos, sin).
  const auto boundary_radius = [&]( const double cph, const double sph ) -> double {
    if( g == GeometryType::CylinderEndOn )
    {
      const double en = det[0]*cph + det[1]*sph;
      const double disc = en*en + od[0]*od[0] - (det[0]*det[0] + det[1]*det[1]);
      return -en + std::sqrt( std::max( 0.0, disc ) );
    }
    double r = 1.0e300;
    if( cph > 0.0 ) r = std::min( r, (od[0] - det[0])/cph );
    if( cph < 0.0 ) r = std::min( r, (-od[0] - det[0])/cph );
    if( sph > 0.0 ) r = std::min( r, (od[1] - det[1])/sph );
    if( sph < 0.0 ) r = std::min( r, (-od[1] - det[1])/sph );
    return r;
  };

  const size_t num_chunks = std::max<size_t>( 1, opt.num_chunks );
  std::vector<double> chunk_sum( num_chunks, 0.0 ), chunk_sq( num_chunks, 0.0 );
  std::vector<size_t> chunk_n( num_chunks, 0 );

  const auto do_chunk = [&]( const size_t chunk )
  {
    std::mt19937_64 rng( opt.seed * 0x9E3779B97F4A7C15ull + 0x100000001B3ull * (chunk + 1) );
    std::uniform_real_distribution<double> U( 0.0, 1.0 );
    const size_t n_lo = (opt.n_samples*chunk)/num_chunks, n_hi = (opt.n_samples*(chunk + 1))/num_chunks;

    const size_t batch = 4096;
    ceelo::ApertureQuadrature q;
    std::vector<double> pre( batch ), p_out;
    std::vector<ceelo::PathSegment> segs;
    double sum = 0.0, sq = 0.0;

    size_t done = n_lo;
    while( done < n_hi )
    {
      const size_t nb = std::min( batch, n_hi - done );
      q.rays.assign( nb, ceelo::KernelRay() );
      q.n_rays_total = static_cast<int>( nb );

      for( size_t i = 0; i < nb; ++i )
      {
        q.rays[i].active_len = 0.0f;   //a sample with no emission leaves the ray empty
        q.rays[i].omega_w = 1.0f;

        // --- emission point and its sampling density q(p)
        double p[3], qdens;
        if( inv_sq )
        {
          const double ph = 2.0*pi*U( rng );
          const double cph = std::cos( ph ), sph = std::sin( ph );
          const double R = boundary_radius( cph, sph );
          const double lam = std::log( 1.0 + R*R/(h0*h0) );
          const double rho = h0*std::sqrt( std::max( 0.0, std::pow( 1.0 + R*R/(h0*h0), U( rng ) ) - 1.0 ) );
          double z, qz;
          const double half = od[(g == GeometryType::Rectangular) ? 2 : 1];
          if( in_situ )
          {
            // Depth from the truncated exponential profile over the column [0, 2 half].
            const double c1 = -std::expm1( -2.0*half/L );
            const double dep = -L*std::log1p( -U( rng )*c1 );
            z = top - dep;
            qz = std::exp( -dep/L ) / (L*c1);
          }else
          {
            z = half*(2.0*U( rng ) - 1.0);
            qz = 1.0/(2.0*half);
          }
          p[0] = det[0] + rho*cph;
          p[1] = det[1] + rho*sph;
          p[2] = z;
          qdens = qz / (pi*(rho*rho + h0*h0)*lam);
        }else
        {
          // Uniform in the outer solid.
          switch( g )
          {
            case GeometryType::Spherical:
            {
              const double r = od[0]*std::cbrt( U( rng ) );
              const double ct = 1.0 - 2.0*U( rng );
              const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
              const double ph = 2.0*pi*U( rng );
              p[0] = r*st*std::cos( ph );
              p[1] = r*st*std::sin( ph );
              p[2] = r*ct;
              break;
            }
            case GeometryType::CylinderEndOn:
            case GeometryType::CylinderSideOn:
            {
              const double r = od[0]*std::sqrt( U( rng ) );
              const double ph = 2.0*pi*U( rng );
              p[0] = r*std::cos( ph );
              p[1] = r*std::sin( ph );
              p[2] = od[1]*(2.0*U( rng ) - 1.0);
              break;
            }
            case GeometryType::Rectangular:
              p[0] = od[0]*(2.0*U( rng ) - 1.0);
              p[1] = od[1]*(2.0*U( rng ) - 1.0);
              p[2] = od[2]*(2.0*U( rng ) - 1.0);
              break;
            case GeometryType::NumGeometryType:
              break;
          }
          qdens = 1.0/vol_outer;
        }

        // --- direction uniform in the cone containing the crystal
        double axis[3], cos_alpha;
        m.direction_cone( p, axis, cos_alpha );
        const double omega_frac = 0.5*(1.0 - cos_alpha);
        const double ct = cos_alpha + (1.0 - cos_alpha)*U( rng );
        const double st = std::sqrt( std::max( 0.0, 1.0 - ct*ct ) );
        const double psi = 2.0*pi*U( rng );
        // Tangent frame about the cone axis (reference vector chosen away from the axis).
        double e1[3], e2[3];
        {
          const double ax = (std::fabs(axis[2]) < 0.9) ? 0.0 : 1.0;
          const double rv[3] = { ax, 0.0, 1.0 - ax };
          e1[0] = axis[1]*rv[2] - axis[2]*rv[1];
          e1[1] = axis[2]*rv[0] - axis[0]*rv[2];
          e1[2] = axis[0]*rv[1] - axis[1]*rv[0];
          const double n1 = std::sqrt( e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2] );
          for( int k = 0; k < 3; ++k )
            e1[k] /= n1;
          e2[0] = axis[1]*e1[2] - axis[2]*e1[1];
          e2[1] = axis[2]*e1[0] - axis[0]*e1[2];
          e2[2] = axis[0]*e1[1] - axis[1]*e1[0];
        }
        double u[3];
        for( int k = 0; k < 3; ++k )
          u[k] = ct*axis[k] + st*(std::cos(psi)*e1[k] + std::sin(psi)*e2[k]);

        // --- per-sample weight without the crystal kernel
        const double rho = m.emission( p );
        double w = 0.0;
        if( rho > 0.0 )
        {
          const Eigen::Vector3d pc = m.to_crystal( p );
          w = rho/qdens * omega_frac * m.transmission( p, u ) * m.prefactor( pc );

          // The single ray through the detector.
          ceelo::KernelRay &kr = q.rays[i];
          const Eigen::Vector3d uc = m.dir_to_crystal( u );
          m.geom.trace_ray( pc, uc, segs );
          std::sort( begin(segs), end(segs), []( const ceelo::PathSegment &a, const ceelo::PathSegment &b ){
            return a.t_start < b.t_start; } );
          kr.omega_w = 1.0f;
          kr.cos_incidence = static_cast<float>( std::fabs( uc.z() ) );
          kr.dir = uc.cast<float>();
          double active = 0.0;
          for( const ceelo::PathSegment &s : segs )
          {
            const double len = s.length();
            if( len <= 1.0e-12 )
              continue;
            if( s.is_scoring )
              active += len;
            if( s.material )
              kr.segs.push_back( { s.material, static_cast<float>(len), s.is_scoring } );
          }
          kr.active_len = static_cast<float>( active );
        }
        pre[i] = w;
      }//for( samples in batch )

      m.resp.fep_line_probabilities( calc.m_energy, q, p_out );
      for( size_t i = 0; i < nb; ++i )
      {
        const double v = pre[i] * p_out[i];
        sum += v;
        sq += v*v;
      }
      done += nb;
    }//while( samples )

    chunk_sum[chunk] = sum;
    chunk_sq[chunk] = sq;
    chunk_n[chunk] = n_hi - n_lo;
  };//do_chunk

  if( opt.multithread && (num_chunks > 1) )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( size_t ch = 0; ch < num_chunks; ++ch )
      pool.post( [&do_chunk,ch](){ do_chunk( ch ); } );
    pool.join();
  }else
  {
    for( size_t ch = 0; ch < num_chunks; ++ch )
      do_chunk( ch );
  }

  double sum = 0.0, sq = 0.0;
  size_t n = 0;
  for( size_t ch = 0; ch < num_chunks; ++ch )
  {
    sum += chunk_sum[ch];
    sq += chunk_sq[ch];
    n += chunk_n[ch];
  }

  McRefResult r;
  r.n = n;
  if( n == 0 )
    return r;
  const double nd = static_cast<double>( n );
  r.value = sum/nd;
  const double var = std::max( 0.0, (sq/nd - r.value*r.value) * nd/(nd - 1.0) );
  r.std_err = std::sqrt( var/nd );
  r.ess = (sq > 0.0) ? (sum*sum/sq) : 0.0;
  return r;
}//reference_random(...)

}//namespace VolRef

#endif //VolumetricReferenceIntegrator_h
