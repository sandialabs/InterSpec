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
#include <string>
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

/** Returns a new lower-channel-energy calibration, with each of the original channel energies
 transformed as:
   E_new[i] = offset + gain * E_orig[i]
 That is, `offset` shifts every channel edge by the same amount (in keV), while `gain` is a
 dimensionless multiplicative scaling of the energy axis with nominal value 1.0 (e.g., gain = 2.0
 doubles the energy of every channel edge).  offset = 0, gain = 1 is the identity.

 @param orig The original lower channel energy calibration; must be a valid
        #SpecUtils::EnergyCalType::LowerChannelEdge calibration.
 @param offset The energy, in keV, to add to every channel edge (applied after the gain scaling).
 @param gain The dimensionless multiplicative scaling of the original channel energies (nominal
        1.0).

 Throws exception if `orig` isnt a valid lower-channel-energy calibration, or if the result would
 not be monotonically increasing (e.g., a non-positive `gain`).
 */
std::shared_ptr<const SpecUtils::EnergyCalibration>
adjust_lower_channel_energy_cal( const std::shared_ptr<const SpecUtils::EnergyCalibration> &orig,
                                 const double offset, const double gain );


/** Fits the {offset, gain} adjustment of a lower-channel-energy calibration (see
 #adjust_lower_channel_energy_cal for their definitions) using peaks with assigned
 nuclide-photopeaks, relative to the ORIGINAL (pre-adjustment) channel energies.

 Note that `peakinfos` channel numbers (RecalPeakInfo::peakMeanBinNumber) should be computed from
 the calibration currently in effect - which shares channel numbers with `orig_cal` - and the
 returned {offset, gain} are then the new *cumulative* adjustments, relative to `orig_cal`.

 @param peakinfos Relevant information collected from peaks.
 @param orig_cal The original lower-channel-energy calibration adjustments are relative to.
 @param fitfor Which of {offset, gain} to fit; must have exactly 2 entries, at least one true.
 @param coefs The fit {offset, gain}; an entry not being fit provides its fixed value on input,
        in which case this vector must have 2 entries on input.
 @param coefs_uncert The uncertainties of the fit {offset, gain}.
 @returns Chi2 of the found solution.

 Throws exception on error.
 */
double fit_energy_cal_lower_channel( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                            const std::shared_ptr<const SpecUtils::EnergyCalibration> &orig_cal,
                                     const std::vector<bool> &fitfor,
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


/** Setup for #fit_energy_cal_ceres.

 The fit parameters are the calibration coefficients (for #SpecUtils::EnergyCalType::LowerChannelEdge
 these are the {offset, gain} adjustments relative to the original channel energies - see
 #adjust_lower_channel_energy_cal), plus, optionally, the offsets of individual deviation pairs.
 */
struct EnergyCalCeresFitSetup
{
  /** Must be Polynomial, UnspecifiedUsingDefaultPolynomial, FullRangeFraction, or
   LowerChannelEdge. */
  SpecUtils::EnergyCalType cal_type;

  size_t num_channels = 0;

  /** Which coefficients to fit; for LowerChannelEdge must have exactly 2 ({offset, gain})
   entries.  At least one entry (of this and/or fit_dev_pair_offsets) must be true. */
  std::vector<bool> fitfor;

  /** Starting (and fixed, where not fit) coefficient values; same size as `fitfor`. */
  std::vector<float> starting_coefs;

  /** The ORIGINAL lower channel energies (num_channels+1 entries); only used, and required, for
   LowerChannelEdge. */
  std::vector<float> lower_channel_energies;

  /** The current deviation pairs; may be empty.  Must be empty for LowerChannelEdge. */
  std::vector<std::pair<float,float>> dev_pairs;

  /** Which deviation pair *offsets* to fit; either empty (fit none), or `dev_pairs.size()`
   entries.  The deviation pair energies are never altered. */
  std::vector<bool> fit_dev_pair_offsets;
};//struct EnergyCalCeresFitSetup


/** Results of #fit_energy_cal_ceres. */
struct EnergyCalCeresFitResult
{
  /** The fit (and passed-through fixed) coefficients, and their uncertainties (zero for fixed
   coefficients, or if the covariance computation failed). */
  std::vector<float> coefs;
  std::vector<float> coef_uncerts;

  /** The deviation pairs: energies unchanged from the input, offsets re-fit where requested.
   Note: if any offsets were fit, and only a single input pair was given, an implicit fixed
   {0,0} pair is prepended (matching the convention SpecUtils applies). */
  std::vector<std::pair<float,float>> dev_pairs;

  /** Uncertainties of the fit deviation pair offsets (zero where not fit); same size as
   `dev_pairs`. */
  std::vector<float> dev_pair_offset_uncerts;

  double chi2 = 0.0;

  /** Non-fatal quality warnings (e.g., the minimizer failed to fully converge, or the fit looks
   under-determined); empty if all looks good. */
  std::string warning_msg;
};//struct EnergyCalCeresFitResult


/** Fits energy calibration coefficients - and optionally deviation pair offsets - to the peaks,
 using a Ceres (Levenberg-Marquardt, auto-differentiated) fit.

 This is the non-linear counterpart of #fit_energy_cal_poly / #fit_energy_cal_frf /
 #fit_energy_cal_lower_channel: those linear fits should be preferred when applicable (they are
 exact), with this function used as the fallback when they fail, or when deviation pair offsets
 are being fit (which the linear fits cant do).

 Throws exception on invalid input, or if the minimization hard-fails; a fit that converges only
 approximately is returned, with `warning_msg` set.
 */
EnergyCalCeresFitResult fit_energy_cal_ceres( const std::vector<RecalPeakInfo> &peakinfos,
                                              const EnergyCalCeresFitSetup &setup );



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
 
 Throws exception if any input #SpecUtils::EnergyCalibration is null or invalid, or if applying
 the difference causes 'other_cal' to become invalid.  'orig_cal' and 'new_cal' should be of the
 same calibration type (polynomial, full range fraction, or lower channel energy).
 
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
