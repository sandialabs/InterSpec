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

// Per-bin unit tests for the GADRAS peak-shape implementation integrated into
// InterSpec/PeakDists_imp.hpp (namespace PeakDists).  Ported from the standalone
// doctest suite scratch/gadras_peak_shape_as_dist/Tests/gadras_peak_dists_tests.cpp.
//
// Reference values (gadras_peak_dists_refdata.inc) are the FORTRAN gold-standard,
// produced by compiling the GADRASw Fortran PeakShapeComputer directly and running
// Distribute() per detector/energy.  The integrated code reproduces these to
// <3e-3 relative / <5e-5 absolute per bin.

#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>

#define BOOST_TEST_MODULE GadrasPeakDists_suite
#include <boost/test/included/unit_test.hpp>

#include "ceres/jet.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakDists_imp.hpp"

using namespace std;

namespace
{
  static const double kRelTol = 3.0e-3;
  static const double kAbsTol = 5.0e-5;   // of the unit-normalized distribution


  // ==========================================================================
  //  Legacy discrete GADRAS peak shape (the 128-point Gaussian-mixture form).
  // --------------------------------------------------------------------------
  //  This is the implementation that used to live in InterSpec/PeakDists_imp.hpp
  //  before it was replaced by the analytic-EMG form.  It reproduces the GADRASw
  //  Fortran shape by summing 128 shifted Gaussians on a fixed zeta grid (a
  //  right-endpoint rectangle-rule quadrature of the exponential tail density).
  //  It is retained here ONLY so we can keep validating against the Fortran gold
  //  standard (MatchesFortranPerBin) and regression-test the new analytic form
  //  against it (AnalyticCloseToDiscrete).  Do not use it in production; the
  //  analytic form (PeakDists::gadras_*) is the maintained implementation.
  //
  //  It also has the PVT ("low photopeak probability") high-tail branch, verbatim
  //  (including its tiny 3E-6 widening and log-extrapolation fudge, which the
  //  analytic skew-normal form does not reproduce).
  //
  //  The grid size and quadrature rule are parameters, so the tests can also refine
  //  it to show the analytic form is its continuous limit.
  // ==========================================================================
  namespace legacy
  {
    static constexpr int kNumGridPoints = 128;   // == SPLINE_POINT_COUNT in the GADRAS code

    inline double snorm_cdf( const double z )
    {
      return 0.5 * (1.0 + std::erf( z * 0.70710678118654752440 ));
    }

    struct GadrasPeakShape
    {
      double sum_skew = 0.0;
      double low_zeta_factor = 1.0;
      double high_zeta_factor = 1.0;
      std::vector<double> gz;   // impulse positions (zeta)
      std::vector<double> w;    // impulse weights
    };

    /** The GADRAS zeta grid (ZetaRange), refined by an integer factor `refine`.

     GADRAS's grid is zeta = sign(u)*(u/Divider)^2 for u = -63, -62, ..., 64 (128 points; Divider
     differing below/above u=0).  Refining keeps the same end-points by stepping u by 1/refine, so
     there are 127*refine + 1 points, with u = 0 at index 63*refine; refine == 1 is exactly GADRAS's.
     */
    inline std::vector<double> build_zeta_grid( const double base_low_skew, const double base_high_skew,
                                                const int refine )
    {
      double low = std::abs( base_low_skew );
      double high = std::abs( base_high_skew );
      if( high > 0.0 )
        low = std::max( low, 0.1 );

      const double divider_lo = 12.0 / std::pow( std::max( 1.0, low ), 0.2 );
      const double divider_hi = 18.0 / std::pow( std::max( 1.0, high ), 0.2 );

      const int npoints = 127*refine + 1;
      std::vector<double> gz( npoints );
      for( int j = 0; j < npoints; ++j )
      {
        const double u = (j - 63.0*refine) / refine;
        const double zeta = u / ((u < 0.0) ? divider_lo : divider_hi);
        gz[j] = (zeta >= 0.0) ? (zeta * zeta) : -(zeta * zeta);
      }
      return gz;
    }

    /** Builds GADRAS's discrete shape.  With `refine == 1` and `midpoint == false` this is exactly
     GADRAS's construction (right-endpoint rectangle rule: each grid interval's mass is the density
     at its upper zeta, placed there); `refine > 1` subdivides GADRAS's grid over the same span, and
     `midpoint == true` evaluates and places the mass at interval midpoints instead, which converges
     to the continuous limit as 1/N^2, not 1/N. */
    inline GadrasPeakShape build_peak_shape( const double energy,
                                             const double low_skew, const double high_skew,
                                             const double low_skew_power, const double high_skew_power,
                                             const double low_skew_extent, const double high_skew_extent,
                                             const PeakDists::GadrasMaterial material,
                                             const int refine = 1,
                                             const bool midpoint = false )
    {
      const bool low_photopeak_probability = (material == PeakDists::GadrasMaterial::LowPhotopeakProbability);
      GadrasPeakShape s;

      // sum_skew (GetSumSkew): raw magnitudes, energy-scaled only when power > 0.
      double skew_p = std::abs( high_skew );
      if( (skew_p > 0.0) && (high_skew_power > 0.0) )
        skew_p *= std::pow( energy / 661.0, high_skew_power );
      double skew_n = std::abs( low_skew );
      if( (skew_n > 0.0) && (low_skew_power > 0.0) )
        skew_n *= std::pow( energy / 661.0, low_skew_power );
      s.sum_skew = std::min( 1.0, (skew_p + skew_n) / 100.0 );

      if( s.sum_skew <= 0.0 )
        return s;

      double low_val  = std::max( 0.0, low_skew );
      double high_val = std::max( 0.0, high_skew );
      if( high_val > 0.0 )
        low_val = std::max( low_val, 0.1 );

      if( (low_val <= 0.0) && (high_val <= 0.0) )
      {
        s.sum_skew = 0.0;
        return s;
      }

      const std::vector<double> gz_edges = build_zeta_grid( low_skew, high_skew, refine );
      const int N = static_cast<int>( gz_edges.size() );
      // Where each interval's density is evaluated and its mass placed.
      std::vector<double> gz = gz_edges;
      if( midpoint )
        for( int i = 1; i < N; ++i )
          gz[i] = 0.5*(gz_edges[i] + gz_edges[i-1]);
      const double *gze = gz_edges.data();

      // First index of the high side (GADRAS's halfPointIdx+1, 1-based); the low side's last
      //  interval ends at zeta = 0, index half-1.
      const int half = 63*refine + 1;
      std::vector<double> gs( N, 0.0 );

      auto slope_scale = []( const double extent ) -> double {
        return (extent >= 0.0) ? (1.0 + extent / 3.0) : std::exp( extent / 3.0 );
      };

      // --- low (left) tail density ---
      if( low_val > 0.0 )
      {
        const double lss = slope_scale( low_skew_extent );
        double sum = 0.0;
        if( material == PeakDists::GadrasMaterial::CZT_CdTe )
        {
          const double sn = 0.1 * lss * low_val;
          for( int i = 1; i < half; ++i )
          {
            gs[i] = (gze[i] - gze[i-1]) * (0.8 * std::exp( gz[i] / sn ) + 0.2 * std::exp( 0.8 * gz[i] / sn ));
            sum += gs[i];
          }
        }
        else
        {
          const double sn = 0.2 * lss * low_val;
          const double fr = 0.04 * low_val;
          for( int i = 1; i < half; ++i )
          {
            gs[i] = (gze[i] - gze[i-1]) * ((1.0 - fr) * std::exp( gz[i] / sn ) + fr * std::exp( 0.4 * gz[i] / sn ));
            sum += gs[i];
          }
        }
        if( !(sum > 0.0) )
        {
          // Tail scale so small every midpoint density underflowed: it is a delta at zero (as the
          //  right-endpoint rule, which samples zeta = 0 itself, would give).
          gs[half-1] = 1.0;
          gz[half-1] = 0.0;
          sum = 1.0;
        }
        for( int i = 0; i < half; ++i )
          gs[i] /= sum;
      }

      // --- high (right) tail density (starts at `half`, matching the Fortran) ---
      if( high_val > 0.0 )
      {
        const double hss = slope_scale( high_skew_extent );
        double sum = 0.0;
        if( material == PeakDists::GadrasMaterial::CZT_CdTe )
        {
          const double sp = 0.1 * hss * high_val;
          for( int i = half; i < N; ++i )
          {
            gs[i] = (gze[i] - gze[i-1]) * (0.8 * std::exp( -gz[i] / sp ) + 0.2 * std::exp( -0.65 * gz[i] / sp ));
            sum += gs[i];
          }
        }
        else if( low_photopeak_probability )
        {
          // PVT-like high tail: Gaussian-CDF differences with a small linear widening,
          // switching to logarithmic extrapolation once that widening term dominates.
          const double divider = high_val / 9.0;
          double gl = 0.5, gintl = 0.5, zeta_last = 0.0;
          for( int i = half; i < N; ++i )
          {
            const double zeta = gze[i] / divider;
            const double gintn = snorm_cdf( zeta );
            const double gn = gintn * (1.0 + 3.0e-6 * hss * zeta);
            if( std::fabs(gintn - gintl)
                < std::fabs( (gintn*zeta - gintl*zeta_last) * 3.0e-6 * hss ) )
            {
              const double gs_intercept = std::log( gs[i-2] );
              const double gs_slope = (std::log(gs[i-1]) - std::log(gs[i-2]))
                                      / (gze[i-1]/divider - gze[i-2]/divider);
              const double zeta_intercept = gze[i-2] / divider;
              for( int j = i; j < N; ++j )
                gs[j] = std::exp( gs_intercept + gs_slope * (gze[j]/divider - zeta_intercept) );
              break;
            }
            gs[i] = gn - gl;
            gl = gn; gintl = gintn; zeta_last = zeta;
          }
          for( int i = half; i < N; ++i )
            sum += gs[i];
        }
        else
        {
          const double sp = 0.2 * hss * high_val;
          for( int i = half; i < N; ++i )
          {
            gs[i] = (gze[i] - gze[i-1]) * std::exp( -gz[i] / sp );
            sum += gs[i];
          }
        }
        if( !(sum > 0.0) )
        {
          gs[half] = 1.0;
          gz[half] = 0.0;
          sum = 1.0;
        }
        for( int i = half; i < N; ++i )
          gs[i] /= sum;
      }

      // Weight the two halves by their relative skew fractions (SCN / SCP).
      const double denom = low_val + high_val;
      const double scn = (denom > 0.0) ? (low_val / denom) : 0.0;
      const double scp = (denom > 0.0) ? (high_val / denom) : 0.0;
      for( int i = 0; i < half; ++i )
        gs[i] *= scn;
      for( int i = half; i < N; ++i )
        gs[i] *= scp;

      // Compact the mixture, dropping negligible weights (scaled so the total dropped stays tiny).
      const double weight_cutoff = 1.0e-9 / refine;
      for( int i = 0; i < N; ++i )
      {
        if( std::fabs( gs[i] ) > weight_cutoff )
        {
          // The PVT tail is built from CDF differences between grid *edges*, so its mass belongs
          //  at the interval it spans; GADRAS places it at the upper edge (as for the others).
          const bool pvt_side = low_photopeak_probability && (i >= half);
          s.gz.push_back( (pvt_side && midpoint && (i > 0)) ? 0.5*(gze[i] + gze[i-1]) : gz[i] );
          s.w.push_back( gs[i] );
        }
      }

      s.low_zeta_factor  = std::pow( 661.0 / energy, std::max( 0.0, low_skew_power ) );
      s.high_zeta_factor = std::pow( 661.0 / energy, std::max( 0.0, high_skew_power ) );

      return s;
    }//build_peak_shape(...)

    inline double peak_shape_cdf( const double zeta, const GadrasPeakShape &s )
    {
      const double gauss = snorm_cdf( zeta );
      if( s.sum_skew <= 0.0 )
        return gauss;

      const double factor = (zeta > 0.0) ? s.high_zeta_factor : s.low_zeta_factor;
      const double z_shape = zeta * factor;

      double shape = 0.0;
      for( size_t j = 0; j < s.w.size(); ++j )
        shape += s.w[j] * snorm_cdf( z_shape - s.gz[j] );

      const double result = (1.0 - s.sum_skew) * gauss + s.sum_skew * shape;
      return std::min( 1.0, std::max( 0.0, result ) );
    }//peak_shape_cdf(...)
  }//namespace legacy

  // A local mirror of the GADRAS detector parameters used to build a shape (the InterSpec
  //  PeakDists API takes the six skew params + material; resolution is used only to derive sigma).
  struct DetParams
  {
    double resolution_offset = 0.0;
    double resolution_661    = 0.0;
    double resolution_power  = 0.0;
    double fwhm_adjustment   = 1.0;

    double low_skew = 0.0, high_skew = 0.0;
    double low_skew_power = 0.0, high_skew_power = 0.0;
    double low_skew_extent = 0.0, high_skew_extent = 0.0;

    PeakDists::GadrasMaterial material = PeakDists::GadrasMaterial::Generic;
    bool low_photopeak_probability = false;
  };

  // The reference-data struct; must match gadras_peak_dists_refdata.inc.
  struct RefCase
  {
    const char *name;
    const char *det_key;
    double      energy;
    double      lo;
    double      hi;
    int         nbins;
    bool        has_high_tail;
    std::vector<double> truth_norm;
  };

#include "gadras_peak_dists_refdata.inc"   // defines kRefCases[]

  DetParams make_params( const std::string &key )
  {
    DetParams p;
    if( key == "identifinder" )            // NaI
    {
      p.resolution_offset = -6.0;  p.resolution_661 = 7.8;     p.resolution_power = 0.55;
      p.low_skew = 0.0;   p.high_skew = 18.0;
      p.low_skew_power = -0.2;   p.high_skew_power = 0.0;
      p.low_skew_extent = -3.0;  p.high_skew_extent = -3.0;
    }
    else if( key == "detective" )          // HPGe
    {
      p.resolution_offset = 1.6;   p.resolution_661 = 0.282001; p.resolution_power = 0.312;
      p.low_skew = 4.61;  p.high_skew = 0.0;
      p.low_skew_power = -0.141;  p.high_skew_power = 0.0;
      p.low_skew_extent = 0.0;    p.high_skew_extent = 0.0;
    }
    else if( key == "d3s" )                // CsI
    {
      p.resolution_offset = -2.47; p.resolution_661 = 7.07;    p.resolution_power = 0.484;
    }
    else if( key == "sam" )                // LaBr3
    {
      p.resolution_offset = 7.0;   p.resolution_661 = 2.6;     p.resolution_power = 0.52;
    }
    else if( key == "czt" )                // CZT
    {
      p.resolution_661 = 1.30176; p.resolution_power = 0.09654; p.resolution_offset = 0.0;
      p.low_skew = 39.64508; p.high_skew = 9.71002;
      p.low_skew_power = 0.63548; p.high_skew_power = 0.07442;
      p.low_skew_extent = 1.05222; p.high_skew_extent = 2.16786;
      p.material = PeakDists::GadrasMaterial::CZT_CdTe;
    }
    else if( key == "pvt" )                // PVT (low photopeak probability)
    {
      p.resolution_offset = 0.0; p.resolution_661 = 0.264; p.resolution_power = 0.674;
      p.low_skew = 0.0; p.high_skew = 11.5;
      p.low_skew_power = 0.01; p.high_skew_power = -0.383;
      p.low_skew_extent = 0.0; p.high_skew_extent = -9.57;
      p.material = PeakDists::GadrasMaterial::LowPhotopeakProbability;
      p.low_photopeak_probability = true;
    }
    else if( key == "rapiscan" )           // RapiScan Metor 6S RPM (PVT, but used here as Generic):
    {                                      //  a long high-E extent, where GADRAS's grid truncates the tail
      p.resolution_offset = -0.0737; p.resolution_661 = 33.4; p.resolution_power = 0.98;
      p.low_skew = 0.0; p.high_skew = 21.9;
      p.low_skew_power = 0.0; p.high_skew_power = 0.0535;
      p.low_skew_extent = 0.0; p.high_skew_extent = 21.5;
    }
    else if( key == "bege" )               // LANL BEGe: a long low-E extent
    {
      p.resolution_offset = 1.0; p.resolution_661 = 0.25; p.resolution_power = 0.5;
      p.low_skew = 5.18; p.high_skew = 0.0;
      p.low_skew_power = 0.0; p.high_skew_power = 0.0;
      p.low_skew_extent = 13.3; p.high_skew_extent = 0.0;
    }
    else if( key == "nanoraider" )         // Kromek nanoRaider-Z CZT: very large low skew
    {
      p.resolution_offset = 0.0; p.resolution_661 = 2.5; p.resolution_power = 0.5;
      p.low_skew = 65.6; p.high_skew = 0.0;
      p.low_skew_power = 0.2; p.high_skew_power = 0.0;
      p.low_skew_extent = 0.0; p.high_skew_extent = 0.0;
      p.material = PeakDists::GadrasMaterial::CZT_CdTe;
    }
    else if( key == "lds_cebr3" )          // LDS CeBr3: negative low magnitude, very negative extents
    {
      p.resolution_offset = 2.0; p.resolution_661 = 4.0; p.resolution_power = 0.5;
      p.low_skew = -5.98; p.high_skew = 25.5;
      p.low_skew_power = 0.0; p.high_skew_power = 0.227;
      p.low_skew_extent = -94.3; p.high_skew_extent = -10.0;
    }
    return p;
  }//make_params(...)

  double det_sigma( const double energy, const DetParams &p )
  {
    return PeakDists::gadras_sigma( energy, p.resolution_offset, p.resolution_661,
                                    p.resolution_power, p.fwhm_adjustment,
                                    p.low_photopeak_probability );
  }

  // Integrate a unit-area GADRAS peak over the given bin edges, by differencing the production
  //  (analytic) CDF.
  std::vector<double> integrate( const double energy, const double sigma, const DetParams &p,
                                 const std::vector<float> &edges )
  {
    const int nbins = static_cast<int>( edges.size() ) - 1;
    std::vector<double> y( std::max(0,nbins), 0.0 );
    if( (nbins <= 0) || (sigma <= 0.0) )
      return y;

    const PeakDists::GadrasPeakShape<double> shape = PeakDists::gadras_build_peak_shape( energy,
                              p.low_skew, p.high_skew, p.low_skew_power, p.high_skew_power,
                              p.low_skew_extent, p.high_skew_extent, p.material );

    double cdf_low = PeakDists::gadras_peak_shape_cdf<double>( (edges[0] - energy)/sigma, shape );
    for( int i = 0; i < nbins; ++i )
    {
      const double cdf_high = PeakDists::gadras_peak_shape_cdf<double>( (edges[i+1] - energy)/sigma, shape );
      y[i] = cdf_high - cdf_low;
      cdf_low = cdf_high;
    }
    return y;
  }//integrate(...)

  // Integrate a unit-area GADRAS peak using the LEGACY discrete form; with the defaults this is
  //  GADRAS's own 128-point construction, used to reproduce the Fortran gold standard.
  std::vector<double> legacy_integrate( const double energy, const double sigma, const DetParams &p,
                                        const std::vector<float> &edges,
                                        const int refine = 1,
                                        const bool midpoint = false )
  {
    const int nbins = static_cast<int>( edges.size() ) - 1;
    std::vector<double> y( std::max(0,nbins), 0.0 );
    if( (nbins <= 0) || (sigma <= 0.0) )
      return y;

    const legacy::GadrasPeakShape shape = legacy::build_peak_shape( energy,
                              p.low_skew, p.high_skew, p.low_skew_power, p.high_skew_power,
                              p.low_skew_extent, p.high_skew_extent, p.material,
                              refine, midpoint );

    double cdf_low = legacy::peak_shape_cdf( (edges[0] - energy)/sigma, shape );
    for( int i = 0; i < nbins; ++i )
    {
      const double cdf_high = legacy::peak_shape_cdf( (edges[i+1] - energy)/sigma, shape );
      y[i] = cdf_high - cdf_low;
      cdf_low = cdf_high;
    }
    return y;
  }//legacy_integrate(...)

  std::vector<float> make_edges( const RefCase &rc )
  {
    std::vector<float> edges( rc.nbins + 1 );
    for( int i = 0; i <= rc.nbins; ++i )
      edges[i] = static_cast<float>( rc.lo + (rc.hi - rc.lo) * (double(i) / rc.nbins) );
    return edges;
  }

  std::vector<double> normalized( const std::vector<double> &v )
  {
    const double s = std::accumulate( v.begin(), v.end(), 0.0 );
    std::vector<double> out( v.size(), 0.0 );
    if( s > 0.0 )
      for( size_t i = 0; i < v.size(); ++i )
        out[i] = v[i] / s;
    return out;
  }

  int count_perbin_failures( const std::vector<double> &got, const std::vector<double> &ref,
                             int &worst_bin, double &worst_rel )
  {
    int fails = 0;
    worst_bin = -1;
    worst_rel = 0.0;
    for( size_t i = 0; i < ref.size(); ++i )
    {
      const double d = std::fabs( got[i] - ref[i] );
      if( d <= kAbsTol )
        continue;
      const double rel = (ref[i] != 0.0) ? (d / std::fabs(ref[i])) : d;
      if( rel > worst_rel ){ worst_rel = rel; worst_bin = int(i); }
      if( rel > kRelTol )
        ++fails;
    }
    return fails;
  }
}//namespace


// The Fortran gold standard is reproduced by the LEGACY discrete (128-point) form, which
//  follows the same rectangle-rule quadrature the Fortran uses.  The analytic-EMG production
//  form is the continuum limit of this and differs by ~1-2% in the tails (see
//  AnalyticCloseToDiscrete), so it is intentionally NOT compared against the Fortran here.
BOOST_AUTO_TEST_CASE( MatchesFortranPerBin )
{
  for( const RefCase &rc : kRefCases )
  {
    const DetParams p = make_params( rc.det_key );
    const double sigma = det_sigma( rc.energy, p );
    const std::vector<float> edges = make_edges( rc );
    const std::vector<double> yn = normalized( legacy_integrate( rc.energy, sigma, p, edges ) );

    int worst_bin; double worst_rel;
    const int fails = count_perbin_failures( yn, rc.truth_norm, worst_bin, worst_rel );
    BOOST_TEST_INFO( "case=" << rc.name << " fails=" << fails
                     << " worst_bin=" << worst_bin << " worst_rel=" << worst_rel );
    BOOST_CHECK_EQUAL( fails, 0 );
  }
}


// The analytic-EMG production form should be close to the legacy discrete form for the
//  non-PVT cases.  They are NOT identical: the discrete form is a right-endpoint
//  rectangle-rule quadrature of the exponential tail density on a fixed 128-point grid,
//  while the analytic form is the exact continuum (infinite-grid) limit of that same
//  density convolved with the Gaussian.  The two therefore agree closely in the core and
//  drift apart by ~1-2% in the far tails, where the discrete grid is coarsest and the
//  quadrature error is largest (the analytic form is the more accurate of the two).
//  We assert a tight tolerance where the (normalized) reference density is appreciable and
//  a looser relative tolerance in the sparse tails.  This includes PVT (the analytic
//  skew-normal high tail vs GADRAS's discrete half-Gaussian), which also checks the PVT
//  analytic form against the Fortran, since the legacy form matches it (MatchesFortranPerBin).
BOOST_AUTO_TEST_CASE( AnalyticCloseToDiscrete )
{
  const double core_abs_tol = 1.0e-3;   // absolute, where the discrete bin has real area
  const double tail_rel_tol = 0.06;     // relative, in the sparse far tails (a few %)
  const double tail_abs_floor = 1.0e-4; // below this normalized area a bin counts as "tail"

  for( const RefCase &rc : kRefCases )
  {
    const DetParams p = make_params( rc.det_key );
    const double sigma = det_sigma( rc.energy, p );
    const std::vector<float> edges = make_edges( rc );

    const std::vector<double> analytic = normalized( integrate( rc.energy, sigma, p, edges ) );
    const std::vector<double> discrete = normalized( legacy_integrate( rc.energy, sigma, p, edges ) );

    int core_fails = 0, tail_fails = 0, worst_bin = -1;
    double worst = 0.0;
    for( size_t i = 0; i < discrete.size(); ++i )
    {
      const double d = std::fabs( analytic[i] - discrete[i] );
      if( discrete[i] >= tail_abs_floor )
      {
        const double rel = d / discrete[i];
        if( rel > worst ){ worst = rel; worst_bin = int(i); }
        if( (d > core_abs_tol) && (rel > tail_rel_tol) )
          ++core_fails;
      }
      else
      {
        // Sparse tail bin: only flag gross absolute disagreement.
        if( d > tail_abs_floor )
          ++tail_fails;
      }
    }

    BOOST_TEST_INFO( "case=" << rc.name << " core_fails=" << core_fails
                     << " tail_fails=" << tail_fails
                     << " worst_rel=" << worst << " worst_bin=" << worst_bin );
    BOOST_CHECK_EQUAL( core_fails, 0 );
    BOOST_CHECK_EQUAL( tail_fails, 0 );
  }
}


BOOST_AUTO_TEST_CASE( HighSkewCarriesTailArea )
{
  for( const RefCase &rc : kRefCases )
  {
    if( !rc.has_high_tail )
      continue;

    const DetParams p = make_params( rc.det_key );
    const double sigma = det_sigma( rc.energy, p );
    const std::vector<float> edges = make_edges( rc );
    const std::vector<double> yn = normalized( integrate( rc.energy, sigma, p, edges ) );

    const double tail_start = rc.energy + 3.0 * sigma;
    double tail_area = 0.0;
    for( int i = 0; i < rc.nbins; ++i )
    {
      const double bc = rc.lo + (rc.hi - rc.lo) * ((i + 0.5) / rc.nbins);
      if( bc >= tail_start )
        tail_area += yn[i];
    }
    // A pure Gaussian has only ~0.00135 of its area beyond +3 sigma; a real high tail is larger.
    BOOST_TEST_INFO( "case=" << rc.name << " tail_area=" << tail_area );
    BOOST_CHECK( tail_area > 0.004 );
  }
}


BOOST_AUTO_TEST_CASE( NonNegativeUnitArea )
{
  for( const RefCase &rc : kRefCases )
  {
    const DetParams p = make_params( rc.det_key );
    const double sigma = det_sigma( rc.energy, p );

    const int N = 4000;
    const double lo = rc.energy - 60.0 * sigma;
    const double hi = rc.energy + 60.0 * sigma;
    std::vector<float> edges( N + 1 );
    for( int i = 0; i <= N; ++i )
      edges[i] = static_cast<float>( lo + (hi - lo) * (double(i) / N) );

    const std::vector<double> y = integrate( rc.energy, sigma, p, edges );

    bool all_nonneg = true;
    for( int i = 0; i < N; ++i )
      all_nonneg = all_nonneg && (y[i] >= -1.0e-15);
    BOOST_TEST_INFO( "case=" << rc.name );
    BOOST_CHECK( all_nonneg );

    const double area = std::accumulate( y.begin(), y.end(), 0.0 );
    BOOST_TEST_INFO( "case=" << rc.name << " area=" << area );
    BOOST_CHECK_CLOSE( area, 1.0, 1.0e-2 );  // within 0.01%
  }
}


BOOST_AUTO_TEST_CASE( PureGaussianMatchesAnalytic )
{
  // With no skew, each bin must equal the difference of the standard normal CDF.
  const char *keys[] = { "d3s", "sam" };
  const double energies[] = { 600.0, 1173.2, 122.06 };

  for( const char *key : keys )
  {
    for( double energy : energies )
    {
      const DetParams p = make_params( key );
      const double sigma = det_sigma( energy, p );

      const int N = 200;
      const double lo = energy - 8.0 * sigma;
      const double hi = energy + 8.0 * sigma;
      std::vector<float> edges( N + 1 );
      for( int i = 0; i <= N; ++i )
        edges[i] = static_cast<float>( lo + (hi - lo) * (double(i) / N) );

      const std::vector<double> y = integrate( energy, sigma, p, edges );

      for( int i = 0; i < N; ++i )
      {
        const double zlo = (edges[i]   - energy) / sigma;
        const double zhi = (edges[i+1] - energy) / sigma;
        const double expected = 0.5 * (std::erf( zhi / std::sqrt(2.0) )
                                       - std::erf( zlo / std::sqrt(2.0) ));
        BOOST_TEST_INFO( "key=" << key << " energy=" << energy << " bin=" << i );
        BOOST_CHECK( std::fabs( y[i] - expected ) <= (1.0e-9 + 1.0e-9*std::fabs(expected)) );
      }
    }
  }
}


// The production array-filling `PeakDists::gadras_integral` (what the fitters call) gives the same
//  per-bin areas as integrating the shape's CDF directly.
BOOST_AUTO_TEST_CASE( ArrayIntegralMatchesCdf )
{
  for( const RefCase &rc : kRefCases )
  {
    const DetParams p = make_params( rc.det_key );
    const double sigma = det_sigma( rc.energy, p );
    const std::vector<float> edges = make_edges( rc );
    const std::vector<double> expected = integrate( rc.energy, sigma, p, edges );

    const double skew[6] = { p.low_skew, p.high_skew, p.low_skew_power, p.high_skew_power,
                             p.low_skew_extent, p.high_skew_extent };
    std::vector<double> got( expected.size(), 0.0 );
    PeakDists::gadras_integral<double>( rc.energy, sigma, 1.0, skew, p.material,
                                        edges.data(), got.data(), got.size() );

    double max_diff = 0.0;
    for( size_t i = 0; i < got.size(); ++i )
      max_diff = (std::max)( max_diff, std::fabs( got[i] - expected[i] ) );
    BOOST_TEST_INFO( "case=" << rc.name << " max_diff=" << max_diff );
    BOOST_CHECK( max_diff < 1.0e-12 );
  }
}//ArrayIntegralMatchesCdf


// The fitters use Ceres autodiff, so the gradient of the per-bin areas w.r.t. the six skew
//  parameters and the mean (which also sets the energy the shape is resolved at) must be right -
//  the shape used to be built in `double`, which silently gave zero gradient for every skew
//  parameter, so no fit could ever move them.
BOOST_AUTO_TEST_CASE( JetGradientMatchesFiniteDifference )
{
  typedef ceres::Jet<double,7> Jet7;  // [0] = mean, [1..6] = the six skew parameters

  struct GradCase
  {
    const char *name;
    PeakDists::GadrasMaterial material;
    double mean, sigma;
    double skew[6];
    bool one_sided_powers;  // powers at zero: compare to a forward difference
  };

  const GradCase cases[] = {
    { "generic", PeakDists::GadrasMaterial::Generic,  300.37, 8.0, { 6.0, 3.0, 0.4, 0.3, 1.5, -2.0 }, false },
    { "czt",     PeakDists::GadrasMaterial::CZT_CdTe, 1001.3, 5.0, { 30.0, 8.0, 0.6, 0.1, 1.0, 2.2 }, false },
    { "low-only",PeakDists::GadrasMaterial::Generic,  186.21, 1.2, { 4.0, 0.0, 0.2, 0.0, -1.0, 0.0 }, false },
    { "power0",  PeakDists::GadrasMaterial::Generic,  300.37, 8.0, { 6.0, 3.0, 0.0, 0.0, 1.5, -2.0 }, true },
    { "pvt",     PeakDists::GadrasMaterial::LowPhotopeakProbability, 477.3, 30.0, { 2.0, 11.5, 0.3, 0.2, 0.5, -1.0 }, false },
    // Long high extent: the high tail is truncated (when USE_GADRAS_TRUNCATION), so its gradient
    //  w.r.t. the magnitude/extent/power includes the truncation terms.
    { "trunc",   PeakDists::GadrasMaterial::Generic,  300.37, 3.0, { 1.0, 21.9, 0.3, 0.2, 0.0, 21.5 }, false },
  };

  for( const GradCase &gc : cases )
  {
    // Bins within +-7.5 sigma, so always inside the integration window (which spans at least +-8
    //  sigma) - otherwise a perturbation could move the window edge across a bin, giving a jump.
    //  The range is offset so no bin edge falls on the mean: when the low and high powers differ,
    //  the shape's density is discontinuous there (each side's zeta is rescaled differently), so a
    //  central difference across it would not be a derivative.
    // (The "trunc" case uses a wider range, into where its truncated tail ends, ~43 sigma, but
    //  inside its integration window, ~51 sigma.)
    const bool wide = (std::string(gc.name) == "trunc");
    const int nbins = wide ? 900 : 300;
    const double lo = gc.mean - 7.33*gc.sigma, hi = gc.mean + (wide ? 47.0 : 7.5)*gc.sigma;
    std::vector<float> edges( nbins + 1 );
    for( int i = 0; i <= nbins; ++i )
      edges[i] = static_cast<float>( lo + (hi - lo)*(double(i)/nbins) );

    const auto eval = [&]( const double mean, const double skew[6] ) -> std::vector<double> {
      std::vector<double> y( nbins, 0.0 );
      PeakDists::gadras_integral<double>( mean, gc.sigma, 1.0, skew, gc.material, edges.data(), y.data(), y.size() );
      return y;
    };

    Jet7 jet_skew[6];
    for( int k = 0; k < 6; ++k )
      jet_skew[k] = Jet7( gc.skew[k], k + 1 );
    std::vector<Jet7> jet_y( nbins, Jet7(0.0) );
    PeakDists::gadras_integral<Jet7>( Jet7(gc.mean, 0), Jet7(gc.sigma), Jet7(1.0), jet_skew, gc.material,
                                      edges.data(), jet_y.data(), jet_y.size() );

    const std::vector<double> y0 = eval( gc.mean, gc.skew );
    for( int i = 0; i < nbins; ++i )
      BOOST_REQUIRE_CLOSE( jet_y[i].a + 1.0, y0[i] + 1.0, 1.0e-10 );

    for( int k = 0; k < 7; ++k )
    {
      // Skip the high-tail parameters when there is no high tail (they then have no effect), and
      //  the high extent for the PVT-like tail (which does not use it).
      if( (gc.skew[1] <= 0.0) && ((k == 2) || (k == 4) || (k == 6)) )
        continue;
      if( (gc.material == PeakDists::GadrasMaterial::LowPhotopeakProbability) && (k == 6) )
        continue;

      const bool forward = gc.one_sided_powers && ((k == 3) || (k == 4));
      const double x = (k == 0) ? gc.mean : gc.skew[k-1];
      const double h = (forward ? 1.0e-7 : 1.0e-5) * (std::max)( 1.0, std::fabs(x) );

      double skew_up[6], skew_down[6];
      std::copy( gc.skew, gc.skew + 6, skew_up );
      std::copy( gc.skew, gc.skew + 6, skew_down );
      double mean_up = gc.mean, mean_down = gc.mean;
      if( k == 0 )
      {
        mean_up += h;
        mean_down -= (forward ? 0.0 : h);
      }else
      {
        skew_up[k-1] += h;
        skew_down[k-1] -= (forward ? 0.0 : h);
      }

      const std::vector<double> y_up = eval( mean_up, skew_up );
      const std::vector<double> y_down = eval( mean_down, skew_down );

      double max_fd = 0.0, max_diff = 0.0;
      for( int i = 0; i < nbins; ++i )
      {
        const double fd = (y_up[i] - y_down[i]) / (forward ? h : 2.0*h);
        max_fd = (std::max)( max_fd, std::fabs(fd) );
        max_diff = (std::max)( max_diff, std::fabs( fd - jet_y[i].v[k] ) );
      }

      BOOST_TEST_INFO( "case=" << gc.name << " par=" << k << " max|fd|=" << max_fd
                       << " max|fd-jet|=" << max_diff );
      BOOST_CHECK( max_fd > 1.0e-6 );  // the parameter really does change the shape
      BOOST_CHECK( max_diff <= ((forward ? 1.0e-3 : 1.0e-4)*max_fd + 1.0e-9) );
    }//for( int k = 0; k < 7; ++k )
  }//for( const GradCase &gc : cases )
}//JetGradientMatchesFiniteDifference


BOOST_AUTO_TEST_CASE( AvoidStationarySkewStart )
{
  const std::vector<bool> fit_amps = { true, true, false, false, false, false };

  // Both GADRAS amplitudes zero: the fitted ones move to the default start
  std::vector<double> values = { 0.0, 0.0, 0.3, 0.0, 1.0, 0.0 };
  PeakDef::avoid_stationary_skew_start( PeakDef::SkewType::GadrasGeneric, values, fit_amps );
  double lower, upper, starting, step;
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::GadrasGeneric, PeakDef::SkewPar0,
                                                lower, upper, starting, step ) );
  BOOST_CHECK_EQUAL( values[0], starting );
  BOOST_CHECK_EQUAL( values[1], starting );
  BOOST_CHECK_EQUAL( values[2], 0.3 );
  BOOST_CHECK_EQUAL( values[4], 1.0 );

  // Only the fitted amplitude moves
  values = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  PeakDef::avoid_stationary_skew_start( PeakDef::SkewType::GadrasCZT, values,
                                        { true, false, false, false, false, false } );
  BOOST_CHECK_EQUAL( values[0], starting );
  BOOST_CHECK_EQUAL( values[1], 0.0 );

  // A non-zero amplitude is not a stationary point - nothing changes (e.g., HPGe has no high tail)
  values = { 3.8, 0.0, 0.0, 0.0, 0.0, 0.0 };
  PeakDef::avoid_stationary_skew_start( PeakDef::SkewType::GadrasGeneric, values, fit_amps );
  BOOST_CHECK_EQUAL( values[0], 3.8 );
  BOOST_CHECK_EQUAL( values[1], 0.0 );

  // Nor is a negative one, which still adds to the tail fraction
  values = { -2.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  PeakDef::avoid_stationary_skew_start( PeakDef::SkewType::GadrasGeneric, values, fit_amps );
  BOOST_CHECK_EQUAL( values[0], -2.0 );
  BOOST_CHECK_EQUAL( values[1], 0.0 );
}//AvoidStationarySkewStart


namespace
{
  // Detectors exercising the parts of parameter space the Fortran reference cases do not: long
  //  extents (where GADRAS's grid truncates the tails), very large magnitudes, PVT, and negative
  //  values straight from shipped Detector.dat files.
  const char * const kExtraKeys[] = { "identifinder", "detective", "czt", "pvt", "rapiscan",
                                      "bege", "nanoraider", "lds_cebr3" };
  const double kExtraEnergies[] = { 186.2, 661.7, 1460.8 };

  bool truncation_negligible( const PeakDists::GadrasPeakShape<double> &s )
  {
    // Whether every tail component ends well inside GADRAS's grid, so truncation makes no difference
    for( int i = 0; i < s.n_low; ++i )
      if( (s.low_reach > 0.0) && (s.low_reach < 15.0*s.low_scale[i]) )
        return false;
    for( int i = 0; i < s.n_high; ++i )
      if( (s.high_reach > 0.0) && (s.high_reach < 15.0*s.high_scale[i]) )
        return false;
    return true;
  }
}//namespace


// The analytic form is the continuous limit of GADRAS's discretization: refining GADRAS's own
//  grid (same zeta range, more points, midpoint rule so it converges quickly), the discrete CDF
//  approaches the analytic one.  With USE_GADRAS_TRUNCATION this holds for every detector,
//  including long-extent ones whose tails run past GADRAS's grid; without it, only where the
//  grid's truncation is negligible.
BOOST_AUTO_TEST_CASE( AnalyticIsContinuumLimit )
{
  const int refine = 16;   // ~2000 grid points
  for( const char *key : kExtraKeys )
  {
    const DetParams p = make_params( key );
    for( const double energy : kExtraEnergies )
    {
      const PeakDists::GadrasPeakShape<double> shape = PeakDists::gadras_build_peak_shape( energy,
                                p.low_skew, p.high_skew, p.low_skew_power, p.high_skew_power,
                                p.low_skew_extent, p.high_skew_extent, p.material );
#if( !USE_GADRAS_TRUNCATION )
      {
        // Check the untruncated shape against GADRAS's reach directly
        const double skew[6] = { p.low_skew, p.high_skew, 0, 0, 0, 0 };
        const pair<double,double> reach = PeakDists::gadras_truncation_reach( skew );
        PeakDists::GadrasPeakShape<double> probe = shape;
        probe.low_reach = reach.first;
        probe.high_reach = reach.second;
        if( !truncation_negligible( probe ) )
          continue;
      }
#endif

      const legacy::GadrasPeakShape discrete = legacy::build_peak_shape( energy,
                                p.low_skew, p.high_skew, p.low_skew_power, p.high_skew_power,
                                p.low_skew_extent, p.high_skew_extent, p.material, refine, true );

      double max_diff = 0.0, worst_z = 0.0;
      for( double z = -170.0; z <= 120.0; z += 0.1 )
      {
        const double a = PeakDists::gadras_peak_shape_cdf<double>( z, shape );
        const double d = legacy::peak_shape_cdf( z, discrete );
        if( std::fabs(a - d) > max_diff )
        {
          max_diff = std::fabs(a - d);
          worst_z = z;
        }
      }

      BOOST_TEST_INFO( "key=" << key << " energy=" << energy << " max|cdf diff|=" << max_diff
                       << " at z=" << worst_z );
      BOOST_CHECK( max_diff < 1.0e-4 );
    }//for( const double energy : kExtraEnergies )
  }//for( const char *key : kExtraKeys )
}//AnalyticIsContinuumLimit


// GADRAS's grid only cuts the tails off for long tails; for the ordinary detectors truncating
//  changes nothing, so the choice of USE_GADRAS_TRUNCATION only matters for the few long ones.
BOOST_AUTO_TEST_CASE( TruncationOnlyMattersForLongTails )
{
  for( const char *key : kExtraKeys )
  {
    const DetParams p = make_params( key );
    const PeakDists::GadrasPeakShape<double> shape = PeakDists::gadras_build_peak_shape( 661.7,
                              p.low_skew, p.high_skew, p.low_skew_power, p.high_skew_power,
                              p.low_skew_extent, p.high_skew_extent, p.material );
    double max_trunc_wt = 0.0;
    for( int i = 0; i < shape.n_low; ++i )
      max_trunc_wt = std::max( max_trunc_wt, shape.low_trunc_wt[i] );
    for( int i = 0; i < shape.n_high; ++i )
      max_trunc_wt = std::max( max_trunc_wt, shape.high_trunc_wt[i] );

    const string k = key;
    const bool long_tail = (k == "rapiscan") || (k == "bege");
    BOOST_TEST_INFO( "key=" << key << " max e^{-R/s}=" << max_trunc_wt );
#if( USE_GADRAS_TRUNCATION )
    if( long_tail )
      BOOST_CHECK( max_trunc_wt > 1.0e-3 );
    else
      BOOST_CHECK( max_trunc_wt < 1.0e-3 );
#else
    BOOST_CHECK( max_trunc_wt == 0.0 );
#endif
  }//for( const char *key : kExtraKeys )

  // The reach itself is GADRAS's ZetaRange end-points
  const double czt[6] = { 39.64508, 9.71002, 0, 0, 0, 0 };
  const pair<double,double> reach = PeakDists::gadras_truncation_reach( czt );
  const double lo_end = 63.0 * std::pow( 39.64508, 0.2 ) / 12.0;
  const double hi_end = 64.0 * std::pow( 9.71002, 0.2 ) / 18.0;
  BOOST_CHECK_CLOSE( reach.first, lo_end*lo_end, 1.0e-10 );
  BOOST_CHECK_CLOSE( reach.second, hi_end*hi_end, 1.0e-10 );

  // A high skew floors |low| at 0.1, which max(1,|low|) then hides; no skew => the 1.0 floor.
  const double hi_only[6] = { 0.0, 18.0, 0, 0, 0, 0 };
  const pair<double,double> reach2 = PeakDists::gadras_truncation_reach( hi_only );
  BOOST_CHECK_CLOSE( reach2.first, std::pow( 63.0/12.0, 2 ), 1.0e-10 );
}//TruncationOnlyMattersForLongTails


// gadras_coverage_limits must give the quantiles of the (continuous-limit) GADRAS shape - checked
//  against bisecting a finely-discretized GADRAS shape - including tails reaching well past the
//  15 sigma RelActAuto used to cap them at.
BOOST_AUTO_TEST_CASE( CoverageLimitsMatchDiscreteQuantiles )
{
  const double p = 1.0E-3;   // 99.9% coverage
  for( const char *key : kExtraKeys )
  {
    const DetParams dp = make_params( key );
    const double skew[6] = { dp.low_skew, dp.high_skew, dp.low_skew_power, dp.high_skew_power,
                             dp.low_skew_extent, dp.high_skew_extent };
    for( const double energy : kExtraEnergies )
    {
      const double sigma = 1.0;
      const pair<double,double> lim = PeakDists::gadras_coverage_limits( energy, sigma, skew, dp.material, p );

#if( !USE_GADRAS_TRUNCATION )
      if( (string(key) == "rapiscan") || (string(key) == "bege") )
        continue;
#endif
      const legacy::GadrasPeakShape discrete = legacy::build_peak_shape( energy,
                                dp.low_skew, dp.high_skew, dp.low_skew_power, dp.high_skew_power,
                                dp.low_skew_extent, dp.high_skew_extent, dp.material, 16, true );
      const auto quantile = [&]( const double target ) -> double {
        double lo = -400.0, hi = 400.0;
        for( int i = 0; i < 80; ++i )
        {
          const double mid = 0.5*(lo + hi);
          if( legacy::peak_shape_cdf( mid, discrete ) < target )
            lo = mid;
          else
            hi = mid;
        }
        return 0.5*(lo + hi);
      };

      const double zlo = quantile( 0.5*p ), zhi = quantile( 1.0 - 0.5*p );
      BOOST_TEST_INFO( "key=" << key << " energy=" << energy << " limits=[" << (lim.first - energy)
                       << ", " << (lim.second - energy) << "] sigma; discrete=[" << zlo << ", " << zhi << "]" );
      // The quantile moves by dz = dCDF/density, so in the far tails a 1E-5 CDF difference can be
      //  a sizable fraction of a sigma; compare relative to the distance from the mean.
      BOOST_CHECK( std::fabs( (lim.first - energy) - zlo ) < 0.05 + 0.01*std::fabs(zlo) );
      BOOST_CHECK( std::fabs( (lim.second - energy) - zhi ) < 0.05 + 0.01*std::fabs(zhi) );
    }//for( const double energy : kExtraEnergies )
  }//for( const char *key : kExtraKeys )

  // A CZT peak at 1.46 MeV really does have its 99.9% extent far beyond 15 sigma.
  const DetParams czt = make_params( "czt" );
  const double czt_skew[6] = { czt.low_skew, czt.high_skew, czt.low_skew_power, czt.high_skew_power,
                               czt.low_skew_extent, czt.high_skew_extent };
  const pair<double,double> czt_lim = PeakDists::gadras_coverage_limits( 1460.8, 1.0, czt_skew, czt.material, p );
  BOOST_CHECK( (1460.8 - czt_lim.first) > 40.0 );

#if( USE_GADRAS_TRUNCATION )
  // ... while a truncated tail stops near where GADRAS's grid does (~43 sigma), rather than ~200.
  const DetParams rs = make_params( "rapiscan" );
  const double rs_skew[6] = { rs.low_skew, rs.high_skew, rs.low_skew_power, rs.high_skew_power,
                              rs.low_skew_extent, rs.high_skew_extent };
  const pair<double,double> rs_lim = PeakDists::gadras_coverage_limits( 661.7, 1.0, rs_skew, rs.material, p );
  BOOST_CHECK( (rs_lim.second - 661.7) < 50.0 );
#endif
}//CoverageLimitsMatchDiscreteQuantiles


// Values GADRAS accepts, and how it treats them (see scratch/GADRAS_skew_par_audit_20261006.md).
BOOST_AUTO_TEST_CASE( NegativeValuesBehaveLikeGadras )
{
  const double energy = 300.0;
  const auto cdf_at = []( const double energy, const double s[6], const PeakDists::GadrasMaterial m,
                          const double z ) -> double {
    const PeakDists::GadrasPeakShape<double> shape = PeakDists::gadras_build_peak_shape( energy,
                                                          s[0], s[1], s[2], s[3], s[4], s[5], m );
    return PeakDists::gadras_peak_shape_cdf<double>( z, shape );
  };

  // A negative power is exactly a power of zero (GADRAS uses MAX(0,power) everywhere).
  const double neg_pow[6] = { 6.0, 3.0, -0.4, -0.8, 1.5, -2.0 };
  const double zero_pow[6] = { 6.0, 3.0, 0.0, 0.0, 1.5, -2.0 };
  for( double z = -20.0; z <= 20.0; z += 0.25 )
    BOOST_CHECK_EQUAL( cdf_at( energy, neg_pow, PeakDists::GadrasMaterial::Generic, z ),
                       cdf_at( energy, zero_pow, PeakDists::GadrasMaterial::Generic, z ) );

  // A negative magnitude still adds to the tail fraction (as |value|), but builds no tail on its
  //  side - so it is NOT the same as zero, nor as +|value|.
  const double neg_low[6]  = { -6.0, 10.0, 0.0, 0.0, 0.0, 0.0 };
  const double zero_low[6] = {  0.0, 10.0, 0.0, 0.0, 0.0, 0.0 };
  const double pos_low[6]  = {  6.0, 10.0, 0.0, 0.0, 0.0, 0.0 };
  const PeakDists::GadrasPeakShape<double> s_neg = PeakDists::gadras_build_peak_shape( energy,
              neg_low[0], neg_low[1], neg_low[2], neg_low[3], neg_low[4], neg_low[5], PeakDists::GadrasMaterial::Generic );
  const PeakDists::GadrasPeakShape<double> s_zero = PeakDists::gadras_build_peak_shape( energy,
              zero_low[0], zero_low[1], zero_low[2], zero_low[3], zero_low[4], zero_low[5], PeakDists::GadrasMaterial::Generic );
  const PeakDists::GadrasPeakShape<double> s_pos = PeakDists::gadras_build_peak_shape( energy,
              pos_low[0], pos_low[1], pos_low[2], pos_low[3], pos_low[4], pos_low[5], PeakDists::GadrasMaterial::Generic );
  BOOST_CHECK( s_neg.sum_skew > s_zero.sum_skew );
  BOOST_CHECK_CLOSE( s_neg.sum_skew, s_pos.sum_skew, 1.0e-12 );
  BOOST_CHECK( std::fabs( cdf_at( energy, neg_low, PeakDists::GadrasMaterial::Generic, -3.0 )
                          - cdf_at( energy, pos_low, PeakDists::GadrasMaterial::Generic, -3.0 ) ) > 1.0e-3 );
  // (vs zero, the difference is the larger tail fraction, which shows up in the high tail)
  BOOST_CHECK( std::fabs( cdf_at( energy, neg_low, PeakDists::GadrasMaterial::Generic, 3.0 )
                          - cdf_at( energy, zero_low, PeakDists::GadrasMaterial::Generic, 3.0 ) ) > 1.0e-3 );

  // ...and GADRAS's own (Fortran-reproducing) discrete form agrees on all of it.
  {
    DetParams p = make_params( "lds_cebr3" );
    const double sigma = 2.0;
    std::vector<float> edges( 401 );
    for( int i = 0; i <= 400; ++i )
      edges[i] = static_cast<float>( energy - 20.0*sigma + 0.25*sigma*i + 0.0123 );
    const std::vector<double> a = normalized( integrate( energy, sigma, p, edges ) );
    const std::vector<double> d = normalized( legacy_integrate( energy, sigma, p, edges, 16, true ) );
    double max_diff = 0.0;
    for( size_t i = 0; i < a.size(); ++i )
      max_diff = std::max( max_diff, std::fabs( a[i] - d[i] ) );
    BOOST_TEST_INFO( "max per-bin diff=" << max_diff );
    BOOST_CHECK( max_diff < 1.0e-4 );
  }

  // Both magnitudes non-positive: no tail is built, so a pure Gaussian (GADRAS: the same, after its
  //  renormalization).
  const double neg_both[6] = { -4.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  for( double z = -6.0; z <= 6.0; z += 0.5 )
    BOOST_CHECK_CLOSE( cdf_at( energy, neg_both, PeakDists::GadrasMaterial::Generic, z ) + 1.0,
                       PeakDists::gadras_std_normal_cdf( z ) + 1.0, 1.0e-12 );

  // The parameter ranges accept what shipped Detector.dat files contain.
  double lower, upper, starting, step;
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::GadrasGeneric, PeakDef::SkewPar0,
                                                lower, upper, starting, step ) );
  BOOST_CHECK( lower <= -5.98 );
  BOOST_REQUIRE( PeakDef::skew_parameter_range( PeakDef::SkewType::GadrasGeneric, PeakDef::SkewPar4,
                                                lower, upper, starting, step ) );
  BOOST_CHECK( lower <= -94.3 );
  BOOST_CHECK( upper >= 21.5 );
}//NegativeValuesBehaveLikeGadras


// GADRAS's PVT ("low photopeak probability") high tail is a half-Gaussian of width high/9 sigma,
//  convolved with the core: a skew-normal.  Check the closed form against direct numerical
//  convolution, and that the extent does not enter it (GADRAS only uses it in a 3E-6 term).
BOOST_AUTO_TEST_CASE( PvtTailIsSkewNormal )
{
  for( const double d : { 0.3, 1.2778, 2.556, 5.0 } )
  {
    for( double z = -6.0; z <= 6.0 + 6.0*d; z += 0.37 )
    {
      // F(z) = Integral_0^inf (2/d) phi(t/d) Phi(z - t) dt, by Simpson's rule
      const int n = 4000;
      const double tmax = 12.0*d, h = tmax/n;
      double sum = 0.0;
      for( int i = 0; i <= n; ++i )
      {
        const double t = i*h;
        const double f = (2.0/d) * std::exp( -0.5*(t/d)*(t/d) ) * 0.3989422804014327
                         * legacy::snorm_cdf( z - t );
        sum += f * ((i == 0) || (i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0));
      }
      const double numeric = sum * h / 3.0;
      const double analytic = PeakDists::gadras_half_gauss_tail_cdf<double>( z, d );
      BOOST_TEST_INFO( "d=" << d << " z=" << z );
      BOOST_CHECK( std::fabs( numeric - analytic ) < 1.0e-9 );
    }
  }

  const double ext0[6] = { 0.0, 11.5, 0.0, 0.0, 0.0, 0.0 };
  const double ext9[6] = { 0.0, 11.5, 0.0, 0.0, 0.0, -9.57 };
  const PeakDists::GadrasPeakShape<double> a = PeakDists::gadras_build_peak_shape( 600.0,
                  ext0[0], ext0[1], ext0[2], ext0[3], ext0[4], ext0[5], PeakDists::GadrasMaterial::LowPhotopeakProbability );
  const PeakDists::GadrasPeakShape<double> b = PeakDists::gadras_build_peak_shape( 600.0,
                  ext9[0], ext9[1], ext9[2], ext9[3], ext9[4], ext9[5], PeakDists::GadrasMaterial::LowPhotopeakProbability );
  BOOST_CHECK( a.pvt_high_tail && b.pvt_high_tail );
  for( double z = -5.0; z <= 20.0; z += 0.5 )
    BOOST_CHECK_EQUAL( PeakDists::gadras_peak_shape_cdf<double>( z, a ), PeakDists::gadras_peak_shape_cdf<double>( z, b ) );
}//PvtTailIsSkewNormal
