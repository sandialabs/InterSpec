#ifndef PeakFitLM_h
#define PeakFitLM_h
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

#include <vector>
#include <memory>
#include <utility>

#include <Wt/WFlags>

// Forward declarations
class PeakDef;
class DetectorPeakResponse;
namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils


namespace PeakFitLM
{

/** By default:
 - peak means are only constrained to be in the ROI
 - peak widths are only constrained to be within a reasonable range, for the spectrum
 - If peak means are closer than 1.25 sigma together, there is a punishment
 - If there are multiple peaks, their FWHM is a linear function of thier mean; the FWHM can vary by +-15% over the ROI

 However, these choices can be overridden
 
 Note: this must be an enum, and not a `enum class` because it doesnt look like Wt 3.x WFlags doent play nicely with enum classes.
 */
enum PeakFitLMOptions
{
  /** By default, if peaks are less than 1.25 sigma away from each other, they start getting punished,
   with the punishment growing linearly as they get closer.
   Specifying this option turns this punishment off.
   */
  DoNotPunishForBeingToClose  = 0x01,

  /** Punishes peaks for their areas being less than less than sqrt(data betwee mean +- 1.75*sigma) .
   As of 20250415, totally not tested.
   */
  PunishForPeakBeingStatInsig = 0x02,

  /** Normally when fitting multiple peaks, the peak FWHM is allowed to vary +-15% throughout the
   ROI, with the FWHM of each peak being a linear function of the fraction of energy through the ROI.
   This option allows the FWHM of each peak to independently vary, within the reasonable range
   for the detector type/spectrum.
   */
  AllPeakFwhmIndependent      = 0x04,

  // TODO: add option to just have FWHM vary with sqrt(energy)

  /** By default, peak means and sigmas are allowed to be fit within ROI and detector limits, however
   specifying this option will reduce the sigma to only change up to about 50% (or if `AllPeakFwhmIndependent`
   is not specified, then peak sigmas can range 50% additional from previous range - not applied on a peak-to-peal basis),
   or the mean to move within 50% of a sigma.

   This is somewhat the analigous option to `PeakFitChi2Fcn::kRefitPeakParameters`.
   */
  MediumRefinementOnly        = 0x08,

  /** Similar to `MediumRefinementOnly`, but limits to 15% of a sigma.
   */
  SmallRefinementOnly         = 0x10,
};//enum PeakFitLMOptions

void fit_peak_for_user_click_LM( std::vector< std::shared_ptr<const PeakDef> > &results,
                                const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                const std::vector< std::shared_ptr<const PeakDef> > &coFitPeaks,
                                const double mean0, const double sigma0,
                                const double area0,
                                const float roiLowerEnergy,
                                const float roiUpperEnergy,
                                const bool isHPGe );

/** Analog of `void fitPeaks(...)`, but using the Ceres based L-M fit method.

 Note different order of funciton arguments, as compared to `fitPeaks(...)`.
 All input peaks must be in the same ROI (e.g., share the same PeakContinuum).
 
 Uses default `PeakFitLMOptions` options (e.g., none of them), unless you specify `is_refit = true`,
 then `PeakFitLMOptions::MediumRefinementOnly` option is specified.
 
 May return fewer peaks than passed in if a peak doesnt pass the `stat_threshold` or `hypothesis_threshold`
 (if these are above 0.0), or if two peaks fit within 1 sigma of each other.
 */
void fit_peaks_LM( std::vector<std::shared_ptr<const PeakDef>> &results,
                  const std::vector<std::shared_ptr<const PeakDef>> input_peaks,
                  std::shared_ptr<const SpecUtils::Measurement> data,
                  const double stat_threshold,
                  const double hypothesis_threshold,
                  const bool is_refit,
                  const bool isHPGe ) throw();



std::vector<std::shared_ptr<const PeakDef>> fit_peaks_in_range_LM( const double x0, const double x1,
                                      const double ncausalitysigma,
                                      const double stat_threshold,
                                      const double hypothesis_threshold,
                                      const std::vector<std::shared_ptr<const PeakDef>> all_peaks,
                                      const std::shared_ptr<const SpecUtils::Measurement> data,
                                      const bool isRefit,
                                      const bool isHPGe );

/** Refit peaks that share an ROI.
 * @param data The data to fit
 * @param detector Currently unused
 * @param inpeaks The peaks to refit - the widths and means of these peaks will serve as initial guesses; they
 *        also dictate which quantities are fit for each peak, the skew and continuum types used, etc.
 * @param fit_options The options to use for the fit.
 * @return The refitted peaks, if fitting was sucessful.
 * 
 * On failure, returns an empty vector.
 */
std::vector<std::shared_ptr<const PeakDef>> refitPeaksThatShareROI_LM(
                                   const std::shared_ptr<const SpecUtils::Measurement> &data,
                                   const std::shared_ptr<const DetectorPeakResponse> &detector,
                                   const std::vector<std::shared_ptr<const PeakDef>> &inpeaks,
                                   const Wt::WFlags<PeakFitLMOptions> fit_options = 0 );

// Need to implement the equivalent of `search_for_peaks(...)` which uses Minuit2 based `AutoPeakSearchChi2Fcn` class.
// Also, the `searchForPeakFromUser(...)` 

}//namespace PeakFitLM


#endif //PeakFitLM_h
