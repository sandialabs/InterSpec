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

/// \TODO: we could probably make a fit_energy_cal_lower_channel_energies that adjusts the offset
///        and gain equivalents for lower channel energy defined calibrations.


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


/** Uses Minuit2 to fit for calibration coefficients.
 
 
 Throws exception on error.
 
 @param peakinfos Relevant information collected from peaks.
 @param fitfor Whether to fit for each coefficient.  This vector must be sized
        for the number of coefficients in the calabration equation (e.g.,
        calibration order).  A true value for an index indicates fit for the
        coefficient.  If any value is false (i.e. dont fit for that parameter)
        then the 'coefs' parameter must be exactly the same size as this vector
        and the value at the corresponding index will be used for that parameter.
 @param nchannels The number of gamma channels in the spectrum the calibration is being evaluated
        for.  Must be at least 7.
 @param eqnType The type of energy calibration being fit for; must be either polynomial or FRF.
 @param startingCoefs The starting value of energy coefficients - must be same size as 'fitfor'.
 @param dev_pairs The non-linear deviation pairs.
 @param coefs Provides the fit coeficients for the energy calibration.
 @param coefs_uncert The uncertainties on the fit parameters.
 @param warning_msg If warnings were encountered (like it isnt a great fit), they will be placed
        here.  Will be empty if all looks good.
 @returns Chi2 of the found solution.
 
 \deprecated Please use #fit_energy_cal_poly.  This function is around while #fit_energy_cal_poly
             continues to be tested.
 */
double fit_energy_cal_iterative( const std::vector<EnergyCal::RecalPeakInfo> &peakInfo,
                                 const size_t nchannels,
                                 const SpecUtils::EnergyCalType eqnType,
                                 const std::vector<bool> fitfor,
                                 std::vector<float> &startingCoefs,
                                 const std::vector<std::pair<float,float>> &devpair,
                                 std::vector<float> &coefs,
                                 std::vector<float> &coefs_uncert,
                                 std::string &warning_msg );



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


/** Reads an input CALp file and returns a valid energy calibration.
 
 @param input Input stream with CALp file information.  If nullptr is returned, this function will seekg the stream back to its tellg position
        but if calibration is successfully parsed, than position will be at the end of information used.
 @param num_channels The number of channels this calibration *might* be applied to; needed to fill in channel energies.  If less than
        two, will return nullptr.  If CALp file is for lower channel energies (e.g., has a "Exact Energies" segment), then if the CALp has
        more channels than the data, the input energy calibration will be truncated to specified \p num_channels.  If CALp has less
        channels than nullptr will be returned.  If energy in CALp is invalid, then nullptr will be returned.
 @param det_name [out] The detector name as given in the CALp file, or empty if not give.  Note that detector name is an InterSpec
        specific extension to CALp files.
 @returns nullptr on error, otherwise a valid energy calibration.
 
 \verbatim
 #PeakEasy CALp File Ver:  4.00
 Offset (keV)           :  1.50000e+00
 Gain (keV / Chan)      :  3.00000e+00
 2nd Order Coef         :  0.00000e+00
 3rd Order Coef         :  0.00000e+00
 4th Order Coef         :  0.00000e+00
 Deviation Pairs        :  5
 7.70000e+01 -1.00000e+00
 1.22000e+02 -5.00000e+00
 2.39000e+02 -5.00000e+00
 6.61000e+02 -2.90000e+01
 2.61400e+03  0.00000e+00
 #END
 \endverbatim
 */
std::shared_ptr<SpecUtils::EnergyCalibration>
energy_cal_from_CALp_file( std::istream &input, const size_t num_channels, std::string &det_name );

/** Writes the given energy calibration object as a CALp file.
 
 @param output The stream to write the output to.
 @param The energy calibration to write.
 @param detector_name The name of the detector - an InterSpec specific extension of the CALp file format.  If blank, wont be written.
 @returns if CALp file was successfully written.
 
 Note, if the energy calibration is Full Range Fraction, then it will be converted to polynomial, and those coefficients written out, but also
 the original FRF coefficients will be written out after the other content - this is a InterSpec specific extension of CALp file format.
 */
bool write_CALp_file( std::ostream &output,
                      const std::shared_ptr<const SpecUtils::EnergyCalibration> &cal,
                      const std::string &detector_name );
}//namespace EnergyCal

#endif //EnergyCal_h
