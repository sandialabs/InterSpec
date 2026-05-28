#ifndef EnergyCal_h
#define EnergyCal_h
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
#include <utility>

// Forward declarations
class PeakDef;
namespace SpecUtils
{
  enum class EnergyCalType : int;
  struct EnergyCalibration;
}//namespace SpecUtils


namespace EnergyCal
{

struct RecalPeakInfo
{
  double peakMean;
  double peakMeanUncert;
  double peakMeanBinNumber;
  double photopeakEnergy;
};//struct RecalPeakInfo



/** Fits polynomial energy equation based on peaks with assigned
 nuclide-photopeaks.
 
 @param peakinfos Relevant information collected from peaks.
 @param fitfor Whether to fit for each coefficient.  This vector must be sized
        for the number of coefficients in the calabration equation (e.g.,
        calibration order).  A true value for an index indicates fit for the
        coefficient.  If any value is false (i.e. dont fit for that parameter)
        then the 'coefs' parameter must be exactly the same size as this vector
        and the value at the corresponding index will be used for that parameter.
 @param nchannels The number of channels of the spectrum.
 @param dev_pairs The non-linear deviation pairs.
 @param coefs Provides the fit coeficients for the energy calibration. Also,
        if any parameter is not being fit for than this vector also provides
        its (fixed) value, and this vector must be same size as the 'fitfor'
        vector.
 @param coefs_uncert The uncertainties on the fit parameters.
 @returns Chi2 of the found solution.
 
 Throws exception on error.
 */
double fit_energy_cal_poly( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                            const std::vector<bool> &fitfor,
                            const size_t nchannels,
                            const std::vector<std::pair<float,float>> &dev_pairs,
                            std::vector<float> &coefs,
                            std::vector<float> &coefs_uncert );

/** Analogous to #fit_energy_cal_poly, but for full range fraction.
 */
double fit_energy_cal_frf( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                           const std::vector<bool> &fitfor,
                           const size_t nchannels,
                           const std::vector<std::pair<float,float>> &dev_pairs,
                           std::vector<float> &coefs,
                           std::vector<float> &coefs_uncert );

/** Given the lower channel energies, will determine the best polynomial coeffcients to reproduce
 the binning.
 
 @param ncoeffs The number of coefficents the answer will contain
 @param lower_channel_energies The lower channel energies to fit to
 @param coefs The fit coefficients will be placed into this vector
 @returns the average absolute value of error (in keV) of each channel.
 
 Throws exception on error.
 */
double fit_poly_from_channel_energies( const size_t ncoeffs,
                                       const std::vector<float> &lower_channel_energies,
                                       std::vector<float> &coefs );

/** Given the lower channel energies, will determine the best full range fraction coeffcients to
reproduce the binning.

 @param ncoeffs The number of coefficents the answer will contain
 @param lower_channel_energies The lower channel energies to fit to; assumes the number of gamma
        channels is one less than the number of entries in this vector.
 @param coefs The fit coefficients will be placed into this vector
 @returns the average absolute value of error (in keV) of each channel.

Throws exception on error.
*/
double fit_full_range_fraction_from_channel_energies( const size_t ncoeffs,
                                             const std::vector<float> &lower_channel_energies,
                                             std::vector<float> &coefs );


/** Uses Ceres (Levenberg-Marquardt, automatic differentiation) to fit calibration coefficients.

 This is the fallback path when #fit_energy_cal_poly / #fit_energy_cal_frf fail (e.g., due to
 deviation pairs interacting awkwardly with the linear least-squares formulation).  It is also
 the only entry point for #SpecUtils::EnergyCalType::LowerChannelEdge, which has no LLS form:
 LowerChannelEdge fits a *relative* offset/gain pair such that

   new_lower_edge[i] = baseline_lower_edge[0] + offset
                       + gain * (baseline_lower_edge[i] - baseline_lower_edge[0])

 where `offset` is the energy delta applied at channel 0 (0 means no shift) and `gain` is the
 span multiplier (1 means identity). This matches the values displayed to the user. The caller
 is responsible for materializing the new lower-channel-edge vector from the returned
 (offset, gain) and the baseline.

 Throws exception on error.

 @param peakInfo Relevant information collected from peaks.
 @param nchannels The number of gamma channels in the spectrum the calibration is being evaluated
        for.  Must be at least 7.
 @param eqnType The type of energy calibration being fit for; Polynomial, FullRangeFraction,
        UnspecifiedUsingDefaultPolynomial, or LowerChannelEdge.
 @param fitfor Whether to fit for each coefficient.  This vector must be sized
        for the number of coefficients in the calibration equation (for LowerChannelEdge: 2).
        A true value for an index indicates fit for the coefficient; a false value holds the
        starting value fixed.
 @param startingCoefs The starting value of energy coefficients - must be same size as 'fitfor'.
        For LowerChannelEdge: [offset, gain] with size 2.
 @param devpair The non-linear deviation pairs (must be empty for LowerChannelEdge).
 @param coefs Provides the fit coeficients for the energy calibration.
 @param coefs_uncert The uncertainties on the fit parameters (extracted via Ceres covariance with
        DENSE_SVD).  Fixed parameters get uncertainty 0.
 @param warning_msg If warnings were encountered (e.g. no full convergence, covariance not
        computable), they will be placed here.  Will be empty if all looks good.
 @param lower_channel_edges Required iff eqnType == LowerChannelEdge; must be the baseline
        lower-channel-edge vector with size nchannels+1.  Must be empty for the other
        calibration types.
 @returns Chi2 of the found solution (= 2 * Ceres final_cost).
 */
double fit_energy_cal_iterative( const std::vector<EnergyCal::RecalPeakInfo> &peakInfo,
                                 const size_t nchannels,
                                 const SpecUtils::EnergyCalType eqnType,
                                 const std::vector<bool> fitfor,
                                 std::vector<float> &startingCoefs,
                                 const std::vector<std::pair<float,float>> &devpair,
                                 std::vector<float> &coefs,
                                 std::vector<float> &coefs_uncert,
                                 std::string &warning_msg,
                                 const std::vector<float> &lower_channel_edges = {} );



/** Translates the input peaks for a given energy calibration change.
 
 This translation keeps peak features (mean, ROI start and end, FWHM, continuum) at the same channel
 numbers, but their energies will correspondingly change.  This is useful for energy calibration
 changes.  Note that this function accounts for peaks that share a continuum.
 
 Will throw exception on error (an invalid calibration, less than 5 channels, etc), and on success
 will return a deque of peaks in the same (energy) order as the input, and with same number of
 entries.
 */
std::deque< std::shared_ptr<const PeakDef> >
translatePeaksForCalibrationChange( const std::deque<std::shared_ptr<const PeakDef>> &inputPeaks,
                              const std::shared_ptr<const SpecUtils::EnergyCalibration> &old_cal,
                              const std::shared_ptr<const SpecUtils::EnergyCalibration> &new_cal );



/** Propagates the difference of energy calibration between 'orig_cal' and 'new_cal' to
 'other_cal', returning the answer.  E.g., if you map what channels coorispond to what other
 channels for orig_cal<==>other_cal, then the returned answer gives new_cal<==>answer.
 
 This function is useful when performing energy calibrations on files that have multiple detectors
 with different energy calibrations.
 
 The reason this is a non-trivial operation is that when you perform a a change in energy
 calibration to one spectrum, and expect a second spectrum with a different starting energy
 calibration (maybe higher-order, or more/less-quadratic, etc) to behave similarly (e.g., a peak at
 626 keV in both spectra to begin with, should end up at 661 keV after the change), what you need to
 do is ensure the relative channels stay fixed to each other (e.x., channel 5 of spectrum 1 may
 correspond to same energy as channel 19.3 of spectrum two, while channel 1024 of spectrum 1
 corresponds to the same energy of channel 999 of spectrum two, etc).  That is what this function
 does.
 
 Throws exception if any input #SpecUtils::EnergyCalibration is null, invalid, or if 'orig_cal'
 or 'new_cal' are not polynomial or full range fraction, or if applying the difference causes
 'other_cal' to become invalid.
 
 @returns a valid #SpecUtils::EnergyCalibration pointer that is the same calibration type, and same
          number of channels as 'other_cal'.
 */
std::shared_ptr<const SpecUtils::EnergyCalibration>
propogate_energy_cal_change( const std::shared_ptr<const SpecUtils::EnergyCalibration> &orig_cal,
                             const std::shared_ptr<const SpecUtils::EnergyCalibration> &new_cal,
                             const std::shared_ptr<const SpecUtils::EnergyCalibration> &other_cal );


/** Internal function for fitting polynomial coefficients from channel-energy pairs.
 
 @param channels_energies Vector of (channel, energy) pairs to fit to
 @param poly_terms Number of polynomial terms to fit
 @returns Vector of polynomial coefficients
 
 Throws exception on error.
 */
std::vector<float> fit_for_poly_coefs( const std::vector<std::pair<double,double>> &channels_energies,
                                      const int poly_terms );

/** Internal function for fitting full range fraction coefficients from channel-energy pairs.
 
 @param channels_energies Vector of (channel, energy) pairs to fit to  
 @param nchannels Number of channels in the spectrum
 @param nterms Number of terms to fit (max 5, with 5th term being 1/(1+60*x))
 @returns Vector of full range fraction coefficients
 
 Throws exception on error.
 */
std::vector<float> fit_for_fullrangefraction_coefs( const std::vector<std::pair<double,double>> &channels_energies,
                                                    const size_t nchannels, 
                                                    const int nterms );

}//namespace EnergyCal

#endif //EnergyCal_h
