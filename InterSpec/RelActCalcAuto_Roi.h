#ifndef RelActCalcAuto_Roi_h
#define RelActCalcAuto_Roi_h
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
#include <vector>
#include <functional>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalcAuto.h"

namespace SpecUtils
{
  struct EnergyCalibration;
}


/** Resolves the ROIs of a `RelActCalcAuto::Options` into the concrete ROIs a fit uses.

 An input `Fixed` ROI is used exactly as given.  The other range types say which gamma lines a ROI
 is to cover, and the ROI's extent comes from the detector's peak shape at those lines (see
 `RelActCalcAuto::RoiEdge`), so one preset adapts to any detector resolution:
 - `LineAnchored`: a single ROI from the lower to the upper line, extended past each of them.
 - `CanBeBrokenUp`: a window around each significant line in the range, limited to the range.
 - `CanExpandForFwhm` (deprecated): like `CanBeBrokenUp`, but not limited to the range.

 A continuum sideband is only useful where there are no peaks, so sidebands stop short of the peaks
 found in the spectrum (see `resolve(...)`).  Where windows overlap: if the regions their peaks
 cover (the lines, plus the peak tails; i.e., without the continuum sidebands) are closer than half
 of their combined sidebands, the windows are merged into one ROI; otherwise they stay separate
 ROIs, and divide the continuum between them in proportion to their sidebands.

 Resolved ROIs are `Fixed`, sorted, and do not overlap - not even by a channel (except that input
 `Fixed` ROIs that touch each other share their boundary channel, as they always have).  They are
 limited to the spectrum, and never overlap an input `Fixed` ROI, which takes precedence.  The bounds of a
 resolved (non-`Fixed`) ROI are placed a quarter of the way into its first and last channels, so the
 nearest-channel rounding RelActCalcAuto uses, and the floor rounding used elsewhere, agree on which
 channels the ROI covers.
 */
namespace RelActCalcAutoRoi
{
  /** The detector peak shape at an energy. */
  struct PeakShape
  {
    double fwhm = 0.0;
    PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;
    std::array<double,6> skew_pars{};
  };//struct PeakShape


  /** A gamma line, or a floating peak, that ROIs are built around. */
  struct Line
  {
    /** True energy of the line; or for a floating peak with `observed_energy`, its energy in the
     spectrum's energy calibration. */
    double energy = 0.0;

    /** Whether the line is worth its own window in a `CanBeBrokenUp` ROI; before a fit is done, the
     lines at a peak found in the spectrum are, and afterwards the lines the fit predicts to be
     statistically significant. */
    bool significant = true;

    /** A floating peak, rather than a source gamma line.  Besides getting a window in any `CanBeBrokenUp`
     (or `CanExpandForFwhm`) range it is within, a floating peak that a `LineAnchored` ROI reaches is
     covered as fully as the ROI's own lines (otherwise whether the peak was within the ROI would depend
     on the detector's resolution): by that ROI, or if the peak is far enough from the ROI's lines that
     two lines would not be merged into one ROI, by a ROI of its own. */
    bool floating_peak = false;

    /** For a floating peak whose energy was read off the spectrum, rather than a known gamma energy:
     it is not moved by the `line_position` mapping. */
    bool observed_energy = false;
  };//struct Line


  /** An energy range, in the energy calibration of the spectrum, the continuum sidebands of ROIs are
   kept out of; e.g., where there is a peak (whether from a modeled source or not). */
  struct SidebandObstacle
  {
    double lower = 0.0, upper = 0.0;
  };//struct SidebandObstacle


  struct ResolvedRoi
  {
    /** Always `Fixed`, with #RelActCalcAuto::RoiRange::auto_continuum false; the continuum type is
     the one to use, or to start from if `info.auto_continuum` is true. */
    RelActCalcAuto::RoiRange roi;

    RelActCalcAuto::RoiResolutionInfo info;
  };//struct ResolvedRoi


  /** How far (keV, never negative) to extend a ROI past a line at `energy`, below it if
   `lower_side`, otherwise above it: to where only `edge.tail_fraction` of the line's peak area lies
   beyond the ROI edge, plus `edge.sideband_fwhm` FWHM.
   */
  double edge_distance( const double energy, const PeakShape &shape,
                        const RelActCalcAuto::RoiEdge &edge, const bool lower_side );


  /** The continuum type for a ROI formed by merging ROIs using `lhs` and `rhs`: the higher polynomial
   order of the two, with a step if either has one (a peak-CDF step, if either uses one).

   `External` is only returned if both are `External`.
   */
  PeakContinuum::OffsetType merged_continuum_type( const PeakContinuum::OffsetType lhs,
                                                   const PeakContinuum::OffsetType rhs );


  /** Resolves `input_rois` into the ROIs to fit.

   @param input_rois The ROIs, as specified by the user.
   @param default_lower_edge, default_upper_edge How far ROIs extend below/above their lines, where
          a ROI does not override it (see `RelActCalcAuto::default_roi_edge(...)`).
   @param peak_shape Returns the peak shape for a line at a true energy.
   @param line_position Maps a true line energy to where its peak sits in the spectrum's energy
          calibration (i.e., applies any fitted energy-calibration adjustment); ROIs cover both
          this position and the true energy, so an adjustment that went astray can not leave a
          line outside its ROI.  May be empty, meaning no adjustment.
   @param energy_cal The spectrum's energy calibration; used to limit ROIs to the spectrum, and to
          place ROI bounds within channels.
   @param lines Candidate lines for `CanBeBrokenUp` (and `CanExpandForFwhm`) ROIs, and the floating
          peaks (see `Line::floating_peak`); need not be sorted or unique.  A floating peak within an
          input `Fixed` ROI is left to it, and one no ROI reaches gets no ROI.
   @param sideband_obstacles Where continuum sidebands should not extend to - typically the peaks
          found in the spectrum.  A ROI's sideband stops short of an obstacle centered outside of
          the region the ROI's peaks cover (an obstacle centered within it is taken to be one of the
          ROI's own peaks); the peak-covering region itself is never reduced.  `Fixed` ROIs are not
          affected.
   @param warnings Descriptions of ROIs that could not be used as specified are appended here.
   @returns The ROIs to fit, sorted by energy.  Input `Fixed` ROIs whose center is outside the
          spectrum are left out (with a warning), as are ROIs whose lines are all outside it.

   Throws if an input ROI is invalid (e.g., its lower energy is above its upper energy), if input
   `Fixed` ROIs overlap each other, or if the energy calibration is invalid.
   */
  std::vector<ResolvedRoi> resolve( const std::vector<RelActCalcAuto::RoiRange> &input_rois,
                                    const RelActCalcAuto::RoiEdge &default_lower_edge,
                                    const RelActCalcAuto::RoiEdge &default_upper_edge,
                                    const std::function<PeakShape(double)> &peak_shape,
                                    const std::function<double(double)> &line_position,
                                    const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                                    const std::vector<Line> &lines,
                                    const std::vector<SidebandObstacle> &sideband_obstacles,
                                    std::vector<std::string> &warnings );
}//namespace RelActCalcAutoRoi

#endif //RelActCalcAuto_Roi_h
