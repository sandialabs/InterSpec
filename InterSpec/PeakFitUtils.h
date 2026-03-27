#ifndef PeakFitUtils_h
#define PeakFitUtils_h
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

#include <deque>
#include <memory>
#include <vector>

// Forward declarations
class PeakDef;
class SpecMeas;
class InterSpec;
namespace SpecUtils
{
  class Measurement;
}


namespace PeakFitUtils
{

/** Gives approximate nominal FWHM for a CsI or NaI detector */
float nai_fwhm_fcn( const float energy );

/** Gives approximate nominal FWHM for a LaBr3 detector */
float labr_fwhm_fcn( const float energy );

/** Gives approximate nominal FWHM for a CZT detector */
float czt_fwhm_fcn( const float energy );

/** Gives approximate nominal FWHM for a HPGe detector */
float hpge_fwhm_fcn( const float energy );


enum class CoarseResolutionType : int
{
  /** NaI, CsI */
  Low,

  /** LaBr3 */
  LaBr,

  /** CZT */
  CZT,

  /** Generic medium resolution (CZT or LaBr3) */
  MedRes,

  /** Not HPGe, but could be NaI, CsI, LaBr, or CZT.
   Returned by the FWHM-based classifier when it can distinguish HPGe from non-HPGe,
   but cannot further resolve which non-HPGe type.
   */
  LowOrMedRes,

  /** HPGe */
  High,

  /** Unknown */
  Unknown
};//enum class CoarseResolutionType

  
CoarseResolutionType coarse_resolution_from_peaks( const std::vector<std::shared_ptr<const PeakDef>> &peaks );

/** Convenience function */
CoarseResolutionType coarse_resolution_from_peaks( const std::deque<std::shared_ptr<const PeakDef>> &peaks );
  
  
/** Tries to guess if the passed in spectrum is from a high-resolution system (i.e., HPGe), or a lower resolution system.

 Primarily intended to determine defaults for peak-fitting.

 Currently looks at number of channels and/or average channel width, but may look at actual spectrum features for ambiguous cases
 in the future.

 For a more detailed classification, see \c classify_det_type or \c coarse_det_type.
 */
bool is_high_res( const std::shared_ptr<const SpecUtils::Measurement> &meas );


/** Tries to determine if the current foreground spectrum is from a high-resolution system (i.e., HPGe), or a lower resolution system.
 Must only be called from the Wt application thread.

 For a more detailed classification, see \c classify_det_type or \c coarse_det_type.
 */
bool is_likely_high_res( InterSpec *viewer );


/** Classifies a spectrum's detector type based on peak FWHMs and optionally spectrum metadata.

 Uses GA-optimized parameters to compare measured peak widths against reference FWHM curves
 for each detector type, with relative distance weighting and confidence thresholds.

 @param peaks Candidate or fitted peaks to classify from
 @param spectrum Optional spectrum for supplementary metadata analysis (channel count, keV/channel).
                 If null, only FWHM-based classification is performed.
 @return The classified detector resolution type
 */
CoarseResolutionType classify_det_type(
  const std::vector<PeakDef> &peaks,
  const std::shared_ptr<const SpecUtils::Measurement> &spectrum = nullptr );


/** Determines the coarse detector resolution type using a three-tier approach:
 1. Check SpecUtils::DetectorType from file parsing (most reliable)
 2. Search instrument metadata strings for detector-type keywords
 3. Use GA-optimized FWHM-based classification on candidate peaks

 @param meas The spectrum measurement to classify
 @param spec The SpecMeas file (for detector_type, instrument strings, etc.)
 @return The classified detector resolution type
 */
CoarseResolutionType coarse_det_type(
  const std::shared_ptr<const SpecUtils::Measurement> &meas,
  const std::shared_ptr<const SpecMeas> &spec );
}//namespace PeakFitUtils

#endif //PeakFitUtils_h
