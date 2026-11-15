#ifndef InterSpec_RelActCalcAuto_LineWindow_imp_hpp
#define InterSpec_RelActCalcAuto_LineWindow_imp_hpp

/* InterSpec: which gamma lines outside a RelActCalcAuto ROI that ROI models, and with what weight.

 A ROI models every line whose nominal energy lies inside it.  It also models lines outside it whose
 peaks reach into it - but only those that matter, since each modelled line costs an integration over
 the ROI's channels.  Which lines matter depends on the fit parameters (the peak width and skew, the
 energy calibration, and the line's amplitude), so a hard include/exclude rule makes the objective jump
 whenever a line crosses the rule's boundary.  Instead each outside line gets a weight that goes smoothly
 (C2) from 1 to 0, multiplying the counts it puts into the ROI.  It is the product of:

   - a fraction weight, from the fraction of the line's own area that falls in the ROI: 0 at or below
     `sm_min_fraction_in_roi`, 1 at or above `sm_full_fraction_in_roi` (ten times both for Crystal Ball
     shapes), smooth in between (in the logarithm).  A line's tail is modelled about as far out as its peak
     shape's 1E-4 coverage, as before these weights.
   - a strength weight, from the counts the line puts into the ROI divided by the data uncertainty of
     the ROI's channel nearest it: 0 at or below `sm_min_counts_in_noise`, 1 at or above
     `sm_full_counts_in_noise`, smooth in between (in the logarithm).  Since a line's counts in the ROI
     are at most those, dropping a line for this can move no residual by more than `sm_min_counts_in_noise`.
     (A strong line's tail is not modelled further out than the fraction weight allows: doing so while the
     fit is still far from the data - e.g., unattenuated low-energy lines early in a fit of a shielded
     source - pulled huge tails into neighbouring ROIs and sent solves on long detours.)
   - a cap weight, from the line's distance to the ROI in units of its own peak sigma: 1 up to the
     shape's cap (as many sigma as RelActCalcAuto ever looks for a peak's tail), 0 beyond
     `(1 + sm_cap_taper_fraction)` times that.

 Each factor is a smooth function of quantities that are themselves smooth in the parameters, so the
 objective is continuous, with a continuous gradient, however the parameters move; and a line is
 skipped only where its weight and the weight's derivatives are exactly zero, so skipping is exact.

 This header holds the pure, data-independent parts, so the production code and the unit tests use the
 same functions.
 */

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace RelActCalcAutoImp
{
namespace LineWindow
{

/** A line outside a ROI is not modelled in it when at most this fraction of its area falls in the ROI... */
constexpr double sm_min_fraction_in_roi = 5.0E-6;

/** ...and is modelled in full (if it matters, see below) when at least this does; this matches the 1E-4 (half
 in each tail) coverage RelActCalcAuto used to decide which lines a ROI models before these weights. */
constexpr double sm_full_fraction_in_roi = 5.0E-5;

/** A line outside a ROI is not modelled in it when the counts it would put into the ROI are at most this
 many times the data uncertainty of the ROI's channel nearest the line. */
constexpr double sm_min_counts_in_noise = 0.01;

/** A line outside a ROI is fully modelled in it when the counts it puts into the ROI are at least this many
 times the data uncertainty of the ROI's channel nearest the line; between the two it is tapered. */
constexpr double sm_full_counts_in_noise = 0.1;

/** The cap weight falls from 1, at the cap, to 0 at `(1 + sm_cap_taper_fraction)` times the cap. */
constexpr double sm_cap_taper_fraction = 0.15;

/** A peak whose side facing a ROI is Gaussian puts nothing into the ROI from further than this many sigma: a
 Gaussian has under 1E-15 of its area beyond it, and PeakDists' channel integrals stop there. */
constexpr double sm_gaussian_side_nsigma = 8.0;

/** The number of channels at each ROI edge whose data uncertainty `edge_noise(...)` takes the smallest of. */
constexpr size_t sm_edge_noise_channels = 3;


inline double scalar_part( const double value )
{
  return value;
}

/** The value part of a `ceres::Jet` (or anything with a scalar `.a`). */
template<typename J>
double scalar_part( const J &value )
{
  return value.a;
}


/** The quintic "smootherstep": 0 for `u <= 0`, 1 for `u >= 1`, and `6u^5 - 15u^4 + 10u^3` between, so that
 it, and its first and second derivatives, are continuous everywhere (both derivatives are zero at 0 and 1). */
template<typename T>
T smoothstep( const T &u )
{
  const double value = scalar_part( u );
  if( !(value > 0.0) )
    return T(0.0);
  if( value >= 1.0 )
    return T(1.0);
  return u*u*u*(u*(u*6.0 - 15.0) + 10.0);
}


/** 0 for `value` at or below `lower` (> 0), 1 at or above `upper`, and smooth in the logarithm of `value`
 between. */
template<typename T>
T log_band_weight( const T &value, const double lower, const double upper )
{
  using std::log;
  const double scalar = scalar_part( value );
  if( !(scalar > lower) )
    return T(0.0);
  if( scalar >= upper )
    return T(1.0);

  return smoothstep( T( (log(value) - std::log(lower)) / (std::log(upper) - std::log(lower)) ) );
}


/** The strength weight for a line that puts `counts_over_noise` (counts divided by the data uncertainty
 of the ROI's channel nearest the line) into a ROI: 0 at or below `sm_min_counts_in_noise`, 1 at or above
 `sm_full_counts_in_noise`, and smooth in the logarithm of `counts_over_noise` between. */
template<typename T>
T strength_weight( const T &counts_over_noise )
{
  return log_band_weight( counts_over_noise, sm_min_counts_in_noise, sm_full_counts_in_noise );
}


/** The fraction weight for a line with `fraction` of its area in a ROI: 0 at or below
 `sm_min_fraction_in_roi*scale`, 1 at or above `sm_full_fraction_in_roi*scale`, smooth in the logarithm
 between.  `scale` is 10 for the Crystal Ball shapes (whose 1E-3 coverage was used before), else 1. */
template<typename T>
T fraction_weight( const T &fraction, const double scale = 1.0 )
{
  return log_band_weight( fraction, scale*sm_min_fraction_in_roi, scale*sm_full_fraction_in_roi );
}


/** The cap weight for a line `nsigma` of its own peak sigmas from a ROI: 1 up to `cap_nsigma`, 0 from
 `(1 + sm_cap_taper_fraction)*cap_nsigma`, smooth between.  A non-positive `cap_nsigma` means no cap. */
template<typename T>
T cap_weight( const T &nsigma, const double cap_nsigma )
{
  if( !(cap_nsigma > 0.0) )
    return T(1.0);

  const double value = scalar_part( nsigma );
  const double taper_width = sm_cap_taper_fraction * cap_nsigma;
  if( value <= cap_nsigma )
    return T(1.0);
  if( !(value < (cap_nsigma + taper_width)) )
    return T(0.0);
  return T(1.0) - smoothstep( T( (nsigma - cap_nsigma) / taper_width ) );
}


/** How far beyond a ROI edge a line can be and still have a non-zero cap weight: the distance `D` with
 `D = (1 + sm_cap_taper_fraction) * cap_nsigma * sigma(edge + direction*D)`, found by fixed-point iteration
 from `sigma(edge)`, and using the larger of the sigmas at the edge and at the far end at each step (so it is
 a bound whenever sigma is monotonic over that interval, as detector resolutions are).

 @param sigma_at Gives the peak sigma (keV) at an energy.
 @param edge_energy The ROI edge.
 @param direction -1 for the region below the ROI, +1 for above.
 @param cap_nsigma The shape's cap, in sigma; non-positive for no cap (then `max_reach` is returned).
 @param max_reach The largest distance worth considering (e.g., to the end of the energy span that has
        lines); also returned if the iteration does not settle (a pathological, steep, trial width).
 */
template<typename SigmaFcn>
double reach_beyond_edge( const SigmaFcn &sigma_at, const double edge_energy, const int direction,
                          const double cap_nsigma, const double max_reach )
{
  if( !(cap_nsigma > 0.0) || !(max_reach > 0.0) )
    return std::max( 0.0, max_reach );

  const double k = (1.0 + sm_cap_taper_fraction) * cap_nsigma;
  const double edge_sigma = sigma_at( edge_energy );
  if( !std::isfinite(edge_sigma) || !(edge_sigma > 0.0) )
    return max_reach;

  double reach = k * edge_sigma;
  for( int iter = 0; iter < 25; ++iter )
  {
    if( reach >= max_reach )
      return max_reach;

    const double far_sigma = sigma_at( edge_energy + direction*reach );
    if( !std::isfinite(far_sigma) )
      return max_reach;

    const double next = k * std::max( edge_sigma, far_sigma );
    if( next <= reach*(1.0 + 1.0E-9) )
      return std::min( max_reach, reach*(1.0 + 1.0E-6) );
    reach = next;
  }//for( fixed-point iterations )

  return max_reach;
}//reach_beyond_edge(...)


/** The noise scale a line's counts in a ROI are compared to: the smallest data uncertainty of the
 `sm_edge_noise_channels` channels at the ROI edge nearest the line (`lower_edge` true for lines below
 the ROI), and at least 1 count.

 @param uncerts The data uncertainty of each of the ROI's `nchannel` channels.
 */
inline double edge_noise( const float * const uncerts, const size_t nchannel, const bool lower_edge )
{
  double noise = std::numeric_limits<double>::infinity();
  const size_t n = std::min( nchannel, sm_edge_noise_channels );
  for( size_t i = 0; i < n; ++i )
  {
    const size_t index = lower_edge ? i : (nchannel - 1 - i);
    const double value = uncerts[index];
    if( std::isfinite(value) )
      noise = std::min( noise, value );
  }

  return std::isfinite(noise) ? std::max( 1.0, noise ) : 1.0;
}//edge_noise(...)

}//namespace LineWindow
}//namespace RelActCalcAutoImp

#endif //InterSpec_RelActCalcAuto_LineWindow_imp_hpp
