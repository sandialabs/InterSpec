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

#include <set>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cassert>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/RelActCalcAuto_Roi.h"
#include "InterSpec/RelActCalcAuto_imp.hpp"

using namespace std;

using RelActCalcAuto::RoiEdge;
using RelActCalcAuto::RoiRange;

namespace
{
  /** Long-tailed shapes (e.g., Crystal Ball with a small power) could otherwise ask for enormous ROIs. */
  const double sm_max_edge_nsigma = 20.0;

  /** Where, within its first/last channel, a resolved ROI's lower/upper energy is placed. */
  const double sm_channel_placement = 0.25;


  /** An energy interval on its way to becoming a resolved ROI. */
  struct Window
  {
    /** Bounds, in the energy calibration of the spectrum. */
    double lower = 0.0, upper = 0.0;

    /** The part of the window the peaks cover (i.e., without the continuum sidebands). */
    double core_lower = 0.0, core_upper = 0.0;

    /** True energies of the lowest and highest lines the window covers. */
    double lower_anchor = 0.0, upper_anchor = 0.0;

    /** If the lower bound is where the sideband shared with the window below was divided. */
    bool shares_lower_sideband = false;

    vector<size_t> inputs;
    PeakContinuum::OffsetType continuum = PeakContinuum::OffsetType::Linear;
    bool auto_continuum = false;
  };//struct Window


  /** The channel number (fractional) of `energy`, limited to the spectrum; some calibration types
   (e.g., lower channel energies) throw for energies outside of the spectrum. */
  double channel_within_spectrum( const double energy, const SpecUtils::EnergyCalibration &cal )
  {
    const size_t nchannel = cal.num_channels();
    if( !(energy > cal.energy_for_channel( 0 )) )
      return 0.0;
    if( !(energy < cal.energy_for_channel( static_cast<double>(nchannel) )) )
      return static_cast<double>( nchannel );
    return cal.channel_for_energy( energy );
  }//channel_within_spectrum(...)


  /** The channels RelActCalcAuto integrates for a ROI with the given bounds; this is the same
   nearest-channel rounding as `RoiRangeChannels::channel_range(...)`. */
  pair<size_t,size_t> nearest_channels( const double lower, const double upper,
                                        const SpecUtils::EnergyCalibration &cal )
  {
    const size_t nchannel = cal.num_channels();
    const auto nearest = [&cal,nchannel]( const double energy ) -> size_t {
      const double channel = std::floor( channel_within_spectrum( energy, cal ) + 0.5 );
      return std::min( static_cast<size_t>( std::max( 0.0, channel ) ), nchannel - 1 );
    };
    return { nearest( lower ), nearest( upper ) };
  }//nearest_channels(...)


  string energy_str( const double energy )
  {
    // Up to six significant figures, without exponent notation for any usual energy.
    char buffer[64];
    snprintf( buffer, sizeof(buffer), "%.6g", energy );
    return buffer;
  }


  void merge_into( Window &into, const Window &from )
  {
    into.upper = std::max( into.upper, from.upper );
    into.lower = std::min( into.lower, from.lower );
    into.core_upper = std::max( into.core_upper, from.core_upper );
    into.core_lower = std::min( into.core_lower, from.core_lower );
    into.lower_anchor = std::min( into.lower_anchor, from.lower_anchor );
    into.upper_anchor = std::max( into.upper_anchor, from.upper_anchor );
    for( const size_t index : from.inputs )
    {
      if( std::find( begin(into.inputs), end(into.inputs), index ) == end(into.inputs) )
        into.inputs.push_back( index );
    }
    std::sort( begin(into.inputs), end(into.inputs) );
    into.continuum = RelActCalcAutoRoi::merged_continuum_type( into.continuum, from.continuum );
    into.auto_continuum = (into.auto_continuum || from.auto_continuum);
  }//merge_into(...)


  /** Sorts windows by lower energy, and resolves the ones that overlap (or touch).

   Windows whose peak coverage regions are closer than half of their combined continuum sidebands are
   merged into one window.  Windows further apart than that, but whose sidebands overlap, stay separate
   ROIs, and divide the shared sideband between them in proportion to their sidebands - so each keeps
   at least half of its own.  (Each ROI still accounts for the tails of the other's peaks.)
   */
  vector<Window> merge_overlapping( vector<Window> windows )
  {
    std::sort( begin(windows), end(windows), []( const Window &lhs, const Window &rhs ){
      return (lhs.lower < rhs.lower) || ((lhs.lower == rhs.lower) && (lhs.upper < rhs.upper));
    } );

    vector<Window> merged;
    for( Window w : windows )
    {
      if( merged.empty() || (w.lower > merged.back().upper) )
      {
        merged.push_back( w );
        continue;
      }

      Window &prev = merged.back();
      const double gap = w.core_lower - prev.core_upper;
      const double prev_sideband = std::max( 0.0, prev.upper - prev.core_upper );
      const double this_sideband = std::max( 0.0, w.core_lower - w.lower );
      const double total_sideband = prev_sideband + this_sideband;

      if( !(total_sideband > 0.0) || (gap < 0.5*total_sideband) )
      {
        merge_into( prev, w );
        continue;
      }

      const double boundary = prev.core_upper + gap*prev_sideband/total_sideband;
      prev.upper = boundary;
      w.lower = boundary;
      w.shares_lower_sideband = true;
      merged.push_back( w );
    }//for( Window w : windows )

    return merged;
  }//merge_overlapping(...)
}//namespace


namespace RelActCalcAutoRoi
{

double edge_distance( const double energy, const PeakShape &shape,
                      const RelActCalcAuto::RoiEdge &edge, const bool lower_side )
{
  if( !(shape.fwhm > 0.0) || std::isinf(shape.fwhm) )
    throw runtime_error( "ROI edge: invalid FWHM (" + std::to_string(shape.fwhm) + " keV) at "
                         + energy_str(energy) + " keV." );

  if( !(edge.tail_fraction > 0.0) || !(edge.tail_fraction < 0.5) )
    throw runtime_error( "ROI edge: the tail fraction must be greater than 0 and less than 0.5." );

  if( !(edge.sideband_fwhm >= 0.0) || std::isinf(edge.sideband_fwhm) )
    throw runtime_error( "ROI edge: the continuum sideband must not be negative." );

  RelActCalcAuto::PeakDefImp<double> peak;
  peak.m_mean = energy;
  peak.m_sigma = shape.fwhm / 2.35482;
  peak.m_amplitude = 1.0;
  peak.m_skew_type = shape.skew_type;
  const size_t num_skew = PeakDef::num_skew_parameters( shape.skew_type );
  for( size_t i = 0; i < num_skew; ++i )
    peak.m_skew_pars[i] = shape.skew_pars[i];

  // `peak_coverage_limits(...)` leaves half of the fraction it is given out in each tail, so asking
  //  for twice our one-sided fraction gives the limit on the side we care about.
  const pair<double,double> limits = peak.peak_coverage_limits( 2.0*edge.tail_fraction, sm_max_edge_nsigma );
  const double coverage = lower_side ? (energy - limits.first) : (limits.second - energy);

  return std::max( 0.0, coverage ) + edge.sideband_fwhm*shape.fwhm;
}//edge_distance(...)


PeakContinuum::OffsetType merged_continuum_type( const PeakContinuum::OffsetType lhs,
                                                 const PeakContinuum::OffsetType rhs )
{
  using OT = PeakContinuum::OffsetType;

  if( lhs == rhs )
    return lhs;

  if( (lhs == OT::External) || (rhs == OT::External) )
    return (lhs == OT::External) ? rhs : lhs;

  // Describe each continuum by its polynomial order, and what kind of step (if any) it has.
  enum class Step { None, Data, PeakCdf };
  struct Form { int order; Step step; bool bilinear; };

  const auto form = []( const OT type ) -> Form {
    switch( type )
    {
      case OT::NoOffset:        return { -1, Step::None, false };
      case OT::Constant:        return { 0, Step::None, false };
      case OT::Linear:          return { 1, Step::None, false };
      case OT::Quadratic:       return { 2, Step::None, false };
      case OT::Cubic:           return { 3, Step::None, false };
      case OT::FlatStep:        return { 0, Step::Data, false };
      case OT::LinearStep:      return { 1, Step::Data, false };
      case OT::BiLinearStep:    return { 1, Step::Data, true };
      case OT::FlatStepCDF:     return { 0, Step::PeakCdf, false };
      case OT::LinearStepCDF:   return { 1, Step::PeakCdf, false };
      case OT::BiLinearStepCDF: return { 1, Step::PeakCdf, true };
      case OT::External:        break;
    }
    assert( 0 );
    return { 1, Step::None, false };
  };//form

  const Form a = form( lhs ), b = form( rhs );
  const int order = std::max( a.order, b.order );
  const Step step = std::max( a.step, b.step );
  const bool bilinear = (a.bilinear || b.bilinear);

  if( step == Step::None )
  {
    switch( order )
    {
      case -1: return OT::NoOffset;
      case 0:  return OT::Constant;
      case 1:  return OT::Linear;
      case 2:  return OT::Quadratic;
      default: return OT::Cubic;
    }
  }//if( step == Step::None )

  // There is no quadratic (or higher) step continuum, so the most flexible stepped form stands in.
  const bool cdf = (step == Step::PeakCdf);
  if( bilinear || (order >= 2) )
    return cdf ? OT::BiLinearStepCDF : OT::BiLinearStep;
  if( order == 1 )
    return cdf ? OT::LinearStepCDF : OT::LinearStep;
  return cdf ? OT::FlatStepCDF : OT::FlatStep;
}//merged_continuum_type(...)


std::vector<ResolvedRoi> resolve( const std::vector<RelActCalcAuto::RoiRange> &input_rois,
                                  const RoiEdge &default_lower_edge,
                                  const RoiEdge &default_upper_edge,
                                  const std::function<PeakShape(double)> &peak_shape,
                                  const std::function<double(double)> &line_position,
                                  const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                                  const std::vector<Line> &lines,
                                  const std::vector<SidebandObstacle> &sideband_obstacles,
                                  std::vector<std::string> &warnings )
{
  if( !energy_cal || !energy_cal->valid() || (energy_cal->num_channels() < 16) )
    throw runtime_error( "ROI resolution: invalid energy calibration." );

  const size_t nchannel = energy_cal->num_channels();
  const double spectrum_lower = energy_cal->energy_for_channel( 0 );
  const double spectrum_upper = energy_cal->energy_for_channel( static_cast<double>(nchannel) );

  const auto position = [&line_position]( const double energy ) -> double {
    return line_position ? line_position( energy ) : energy;
  };

  // Validate the input ROIs.
  for( size_t index = 0; index < input_rois.size(); ++index )
  {
    const RoiRange &roi = input_rois[index];
    const bool single_line_ok = (roi.range_limits_type == RoiRange::RangeLimitsType::LineAnchored);

    if( !(roi.lower_energy < roi.upper_energy) && !(single_line_ok && (roi.lower_energy == roi.upper_energy)) )
      throw runtime_error( "Energy range lower value (" + SpecUtils::printCompact(roi.lower_energy, 3)
                           + " keV) is larger or equal to upper value ("
                           + SpecUtils::printCompact(roi.upper_energy, 3) + " keV)" );

    if( roi.lower_energy < 0.0 )
      throw runtime_error( "Energy range, [" + SpecUtils::printCompact(roi.lower_energy, 3)
                           + " - " + SpecUtils::printCompact(roi.upper_energy, 3)
                           + "], extends below zero keV." );

    // Only fixed ROIs are required to not overlap each other; resolved ROIs that overlap are merged.
    if( roi.range_limits_type != RoiRange::RangeLimitsType::Fixed )
      continue;

    for( size_t other_index = index + 1; other_index < input_rois.size(); ++other_index )
    {
      const RoiRange &other = input_rois[other_index];
      if( (other.range_limits_type == RoiRange::RangeLimitsType::Fixed)
          && (roi.lower_energy < other.upper_energy)
          && (other.lower_energy < roi.upper_energy) )
      {
        throw runtime_error( "RelActAutoCostFcn: input energy ranges are overlapping ["
                             + std::to_string(roi.lower_energy) + ", " + std::to_string(roi.upper_energy)
                             + "] and [" + std::to_string(other.lower_energy) + ", "
                             + std::to_string(other.upper_energy) + "]." );
      }
    }//for( loop over other input ROIs )
  }//for( loop over input ROIs )


  // Fixed ROIs are used exactly as given.
  vector<ResolvedRoi> answer;
  for( size_t index = 0; index < input_rois.size(); ++index )
  {
    const RoiRange &roi = input_rois[index];
    if( roi.range_limits_type != RoiRange::RangeLimitsType::Fixed )
      continue;

    const double mid_energy = 0.5*(roi.lower_energy + roi.upper_energy);
    if( !(mid_energy >= spectrum_lower) || !(mid_energy <= spectrum_upper) )
    {
      warnings.push_back( "Not using the [" + SpecUtils::printCompact(roi.lower_energy, 3)
                          + " - " + SpecUtils::printCompact(roi.upper_energy, 3) + "] ROI because over"
                          + " half of it is not within the spectrums energy range ("
                          + energy_str(spectrum_lower) + " - " + energy_str(spectrum_upper) + " keV)." );
      continue;
    }

    ResolvedRoi resolved;
    resolved.roi = roi;
    resolved.roi.auto_continuum = false;
    resolved.roi.lower_edge = resolved.roi.upper_edge = RoiRange::EdgeOverride{};
    resolved.info.input_roi_indices.push_back( index );
    resolved.info.auto_continuum = roi.auto_continuum;
    answer.push_back( std::move(resolved) );
  }//for( loop over input ROIs )


  // Where the peak for a line sits in the spectrum.
  const auto line_pos = [&position]( const Line &line ) -> double {
    return line.observed_energy ? line.energy : position( line.energy );
  };

  // A window for input ROI `index` covering the lines from `lower_line` to `upper_line` (true
  //  energies), whose peaks sit at `lower_pos` and `upper_pos` in the spectrum.  The window covers
  //  both a line's true energy and where the energy-calibration adjustment puts it, so an adjustment
  //  that went astray (e.g., a deviation pair at its limit, where few lines constrain it) can not
  //  leave the line outside its ROI; the next fit then sees the peak, and can correct the adjustment.
  const auto make_window = [&]( const size_t index, const double lower_line, const double lower_pos,
                                const double upper_line, const double upper_pos ) -> Window {
    const RoiRange &roi = input_rois[index];
    const RoiEdge lower_edge = roi.lower_edge.apply_to( default_lower_edge );
    const RoiEdge upper_edge = roi.upper_edge.apply_to( default_upper_edge );
    const PeakShape lower_shape = peak_shape( lower_line ), upper_shape = peak_shape( upper_line );
    const RoiEdge lower_coverage{ lower_edge.tail_fraction, 0.0 }, upper_coverage{ upper_edge.tail_fraction, 0.0 };

    Window w;
    w.lower_anchor = lower_line;
    w.upper_anchor = upper_line;
    w.core_lower = std::min( lower_pos, lower_line ) - edge_distance( lower_line, lower_shape, lower_coverage, true );
    w.core_upper = std::max( upper_pos, upper_line ) + edge_distance( upper_line, upper_shape, upper_coverage, false );
    w.lower = w.core_lower - lower_edge.sideband_fwhm*lower_shape.fwhm;
    w.upper = w.core_upper + upper_edge.sideband_fwhm*upper_shape.fwhm;

    // The sidebands are for the continuum, so they stop short of peaks - e.g., unmodeled x-rays a
    //  little below a line; an obstacle centered within the peak-covering region is one of our peaks.
    for( const SidebandObstacle &obstacle : sideband_obstacles )
    {
      const double center = 0.5*(obstacle.lower + obstacle.upper);
      if( (center >= w.core_lower) && (center <= w.core_upper) )
        continue;

      if( (center < w.core_lower) && (obstacle.upper > w.lower) )
        w.lower = std::min( obstacle.upper, w.core_lower );
      else if( (center > w.core_upper) && (obstacle.lower < w.upper) )
        w.upper = std::max( obstacle.lower, w.core_upper );
    }//for( const SidebandObstacle &obstacle : sideband_obstacles )

    w.inputs.push_back( index );
    w.continuum = roi.continuum_type;
    w.auto_continuum = roi.auto_continuum;
    return w;
  };//make_window


  // Windows around the lines each non-fixed ROI is to cover.
  vector<Window> windows;
  for( size_t index = 0; index < input_rois.size(); ++index )
  {
    const RoiRange &roi = input_rois[index];
    if( roi.range_limits_type == RoiRange::RangeLimitsType::Fixed )
      continue;

    switch( roi.range_limits_type )
    {
      case RoiRange::RangeLimitsType::Fixed:
        assert( 0 );
        break;

      case RoiRange::RangeLimitsType::LineAnchored:
      {
        // Lines partly outside of the spectrum are fit to where the spectrum ends.
        if( !(position(roi.upper_energy) >= spectrum_lower) || !(position(roi.lower_energy) <= spectrum_upper) )
        {
          warnings.push_back( "Not using the ROI for the " + energy_str(roi.lower_energy)
                              + ((roi.lower_energy == roi.upper_energy) ? string()
                                                              : (" to " + energy_str(roi.upper_energy)))
                              + " keV lines, as they are not within the spectrums energy range ("
                              + energy_str(spectrum_lower) + " - " + energy_str(spectrum_upper) + " keV)." );
          break;
        }

        windows.push_back( make_window( index, roi.lower_energy, position(roi.lower_energy),
                                        roi.upper_energy, position(roi.upper_energy) ) );

        // A single ROI spanning much of the spectrum (many FWHM, and a factor of three in energy) is
        //  almost always meant to be split around its lines; its continuum can not follow the spectrum
        //  over the range (and the fit then widens peaks to try to).
        const Window &w = windows.back();
        const double mid_fwhm = peak_shape( 0.5*(roi.lower_energy + roi.upper_energy) ).fwhm;
        if( (mid_fwhm > 0.0) && ((w.upper - w.lower) > 25.0*mid_fwhm) && (w.upper > 3.0*std::max( w.lower, 1.0 )) )
        {
          warnings.push_back( "The ROI for the " + energy_str(roi.lower_energy) + " to "
                              + energy_str(roi.upper_energy) + " keV lines is fit as a single ROI "
                              + energy_str( std::round( w.upper - w.lower ) ) + " keV wide, which its continuum"
                              " can not follow; to instead make ROIs around the lines in it, use \"Split by"
                              " lines\" (CanBeBrokenUp)." );
        }
        break;
      }//case LineAnchored

      case RoiRange::RangeLimitsType::CanBeBrokenUp:
      case RoiRange::RangeLimitsType::CanExpandForFwhm:
      {
        const bool clip_to_range = (roi.range_limits_type == RoiRange::RangeLimitsType::CanBeBrokenUp);

        // The range is in true energy; where its lines sit in the spectrum may be a little outside it.
        const double clip_lower = std::min( roi.lower_energy, position(roi.lower_energy) );
        const double clip_upper = std::max( roi.upper_energy, position(roi.upper_energy) );

        vector<Window> line_windows;
        for( const Line &line : lines )
        {
          if( !line.significant || !(line.energy >= roi.lower_energy) || !(line.energy <= roi.upper_energy) )
            continue;

          const double pos = line_pos( line );
          if( !(pos >= spectrum_lower) || !(pos <= spectrum_upper) )
            continue;

          Window w = make_window( index, line.energy, pos, line.energy, pos );
          if( clip_to_range )
          {
            w.lower = std::max( w.lower, clip_lower );
            w.upper = std::min( w.upper, clip_upper );
            w.core_lower = std::max( w.core_lower, w.lower );
            w.core_upper = std::min( w.core_upper, w.upper );
          }

          if( w.lower < w.upper )
            line_windows.push_back( std::move(w) );
        }//for( const Line &line : lines )

        if( line_windows.empty() )
          warnings.push_back( "No significant gamma lines were found in the "
                              + SpecUtils::printCompact(roi.lower_energy, 4) + " to "
                              + SpecUtils::printCompact(roi.upper_energy, 4) + " keV range." );

        windows.insert( end(windows), begin(line_windows), end(line_windows) );
        break;
      }//case CanBeBrokenUp / CanExpandForFwhm
    }//switch( roi.range_limits_type )
  }//for( loop over input ROIs )


  // A floating peak that a line-anchored ROI's window reaches gets a window of its own (using that
  //  ROI's edges), so the ROI is extended to cover it as fully as its own lines; merging (below) then
  //  treats it like any other neighboring line.  Floating peaks within a fixed ROI, or within a range
  //  their windows were already made for (above), are left as they are.
  vector<Window> floating_windows;
  for( const Line &line : lines )
  {
    if( !line.floating_peak )
      continue;

    const double pos = line_pos( line );
    if( !(pos >= spectrum_lower) || !(pos <= spectrum_upper) )
      continue;

    const bool handled = std::any_of( begin(input_rois), end(input_rois), [&line]( const RoiRange &roi ) {
      return (roi.range_limits_type != RoiRange::RangeLimitsType::LineAnchored)
             && (line.energy >= roi.lower_energy) && (line.energy <= roi.upper_energy);
    } );
    if( handled )
      continue;

    for( const Window &reaching : windows )
    {
      assert( reaching.inputs.size() == 1 );
      const size_t index = reaching.inputs.front();
      if( input_rois[index].range_limits_type != RoiRange::RangeLimitsType::LineAnchored )
        continue;

      Window w = make_window( index, line.energy, pos, line.energy, pos );
      if( (w.lower < reaching.upper) && (reaching.lower < w.upper) )
      {
        floating_windows.push_back( std::move(w) );
        break;
      }
    }//for( const Window &reaching : windows )
  }//for( const Line &line : lines )

  windows.insert( end(windows), begin(floating_windows), end(floating_windows) );


  // Fixed ROIs take precedence over the windows that run into them; a window completely containing a
  //  fixed ROI is split around it, as long as each piece still covers some of its lines.  This is done
  //  before merging windows, so windows on either side of a fixed ROI are not merged through it.
  vector<ResolvedRoi> fixed_rois = answer;
  std::sort( begin(fixed_rois), end(fixed_rois), []( const ResolvedRoi &lhs, const ResolvedRoi &rhs ){
    return lhs.roi.lower_energy < rhs.roi.lower_energy;
  } );

  set<pair<size_t,size_t>> overlap_warned; //{input ROI index, fixed_rois index}
  for( size_t fixed_index = 0; fixed_index < fixed_rois.size(); ++fixed_index )
  {
    const RoiRange &fixed = fixed_rois[fixed_index].roi;

    vector<Window> clipped;
    for( Window w : windows )
    {
      if( (w.upper <= fixed.lower_energy) || (w.lower >= fixed.upper_energy) )
      {
        clipped.push_back( std::move(w) );
        continue;
      }

      // Only worth mentioning if the fixed ROI takes over some of the lines this window was for, not
      //  just some of its continuum sideband.
      assert( w.inputs.size() == 1 );
      const RoiRange &input = input_rois[w.inputs.front()];
      if( (w.upper_anchor > fixed.lower_energy) && (w.lower_anchor < fixed.upper_energy)
         && overlap_warned.insert( {w.inputs.front(), fixed_index} ).second )
      {
        const bool lines = (input.range_limits_type == RoiRange::RangeLimitsType::LineAnchored);
        warnings.push_back( "The ROI for " + string(lines ? "the " : "") + energy_str(input.lower_energy)
                            + ((input.lower_energy == input.upper_energy) ? string() : (" to " + energy_str(input.upper_energy)))
                            + " keV" + string(lines ? " lines" : "") + " overlaps the fixed ["
                            + energy_str(fixed.lower_energy) + " - " + energy_str(fixed.upper_energy)
                            + " keV] ROI, which was used instead where they overlap." );
      }

      // Pieces of the window below and above the fixed ROI, and the lines each still covers.
      if( (w.lower < fixed.lower_energy) && (w.lower_anchor < fixed.lower_energy) )
      {
        Window below = w;
        below.upper = fixed.lower_energy;
        below.core_upper = std::min( below.core_upper, below.upper );
        below.upper_anchor = std::min( w.upper_anchor, fixed.lower_energy );
        clipped.push_back( std::move(below) );
      }

      if( (w.upper > fixed.upper_energy) && (w.upper_anchor > fixed.upper_energy) )
      {
        Window above = w;
        above.lower = fixed.upper_energy;
        above.core_lower = std::max( above.core_lower, above.lower );
        above.lower_anchor = std::max( w.lower_anchor, fixed.upper_energy );
        clipped.push_back( std::move(above) );
      }
    }//for( Window w : windows )

    windows.swap( clipped );
  }//for( loop over fixed_rois )

  windows = merge_overlapping( windows );


  // Limit the windows to the spectrum, and place their bounds a quarter of the way into their
  //  first and last channels.
  vector<pair<size_t,size_t>> window_channels;
  vector<Window> kept_windows;
  for( Window &w : windows )
  {
    w.lower = std::max( w.lower, std::max( spectrum_lower, 0.0 ) );
    w.upper = std::min( w.upper, spectrum_upper );
    if( !(w.lower < w.upper) )
      continue;

    const double lower_channel = std::max( 0.0, std::floor( channel_within_spectrum( w.lower, *energy_cal ) ) );
    const double upper_channel = std::max( 0.0, std::floor( channel_within_spectrum( w.upper, *energy_cal ) ) );
    size_t first = std::min( static_cast<size_t>(lower_channel), nchannel - 1 );
    size_t last = std::min( static_cast<size_t>(upper_channel), nchannel - 1 );
    if( last < first )
      std::swap( first, last );

    // With a negative energy offset, the bound placed in the channel containing zero keV is below zero.
    while( (first < last) && (energy_cal->energy_for_channel( first + sm_channel_placement ) < 0.0) )
      first += 1;

    window_channels.emplace_back( first, last );
    kept_windows.push_back( w );
  }//for( Window &w : windows )

  // Keep windows from sharing channels with fixed ROIs, or with each other.
  vector<pair<size_t,size_t>> fixed_channels;
  for( const ResolvedRoi &fixed : fixed_rois )
    fixed_channels.push_back( nearest_channels( fixed.roi.lower_energy, fixed.roi.upper_energy, *energy_cal ) );

  for( size_t i = 0; i < kept_windows.size(); ++i )
  {
    pair<size_t,size_t> &channels = window_channels[i];
    for( const pair<size_t,size_t> &fixed : fixed_channels )
    {
      if( (channels.second < fixed.first) || (channels.first > fixed.second) )
        continue;

      // The window was already clipped to the fixed ROI in energy, so only the channel(s) at the
      //  boundary can be shared; give them to the fixed ROI.
      if( channels.first >= fixed.first )
        channels.first = fixed.second + 1;
      else
        channels.second = (fixed.first > 0) ? (fixed.first - 1) : 0;
    }
  }//for( loop over kept_windows )

  vector<Window> final_windows;
  vector<pair<size_t,size_t>> final_channels;
  for( size_t i = 0; i < kept_windows.size(); ++i )
  {
    pair<size_t,size_t> &channels = window_channels[i];

    // Windows that divided a shared sideband meet within a channel; that channel goes to the lower one,
    //  unless that would leave too little of the upper one, which is then merged (below).  Other windows
    //  meeting in a channel are merged too - typically their sidebands stopped short of each other's
    //  peaks, leaving no continuum between them.  (So windows that just touch are merged, while ones
    //  overlapping a little more divide their sidebands; the presets were tuned with this behavior.)
    if( kept_windows[i].shares_lower_sideband && !final_windows.empty()
       && (channels.first == final_channels.back().second) && (channels.second > (channels.first + 1)) )
    {
      channels.first += 1;
    }

    if( channels.second <= channels.first )
      continue;

    if( !final_windows.empty() && (channels.first <= final_channels.back().second) )
    {
      merge_into( final_windows.back(), kept_windows[i] );
      final_channels.back().second = std::max( final_channels.back().second, channels.second );
      continue;
    }

    final_windows.push_back( kept_windows[i] );
    final_channels.push_back( channels );
  }//for( loop over kept_windows )

  for( size_t i = 0; i < final_windows.size(); ++i )
  {
    const Window &w = final_windows[i];
    const pair<size_t,size_t> &channels = final_channels[i];

    ResolvedRoi resolved;
    resolved.roi.lower_energy = energy_cal->energy_for_channel( channels.first + sm_channel_placement );
    resolved.roi.upper_energy = energy_cal->energy_for_channel( channels.second + sm_channel_placement );
    resolved.roi.continuum_type = w.continuum;
    resolved.roi.auto_continuum = false;
    resolved.roi.range_limits_type = RoiRange::RangeLimitsType::Fixed;
    resolved.info.input_roi_indices = w.inputs;
    resolved.info.lower_anchor_energy = w.lower_anchor;
    resolved.info.upper_anchor_energy = w.upper_anchor;
    resolved.info.auto_continuum = w.auto_continuum;
    answer.push_back( std::move(resolved) );
  }//for( loop over final windows )

  std::sort( begin(answer), end(answer), []( const ResolvedRoi &lhs, const ResolvedRoi &rhs ){
    return lhs.roi.lower_energy < rhs.roi.lower_energy;
  } );

  return answer;
}//resolve(...)

}//namespace RelActCalcAutoRoi
