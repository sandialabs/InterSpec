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

#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakFitUtils.h"
#include "SpecUtils/EnergyCalibration.h"


namespace PeakFitUtils
{

bool is_high_res( const std::shared_ptr<const SpecUtils::Measurement> &meas )
{
  // I dont think I've seen a HPGe spectrum straight from the MCA with less than 4096, but we'll be
  //  conservative and allow HPGe to go down to 2048 channels
  const size_t min_high_res_chan = 2048;
  
  // I have seen some poorly tuned/defined, one-off type MCAs have a ton of channels for a low
  //  resolution scintillator, but we'll be reasonable and say low-res has max of 4096 channels,
  //  which should account for (nearly all?) the CZT and LaBr detectors
  const size_t max_low_res_chan = 4096;
  
  // If channel counts leaves things ambiguous, we'll look at the width of channels to make a choice
  //   9 MeV spectrum with 9k channels: 1.098 kev/chan.
  //   Which, 1024 channels * 1.1 --> 1127 keV, which is less than 1400 keV almost all dets go to.
  const float max_high_res_chann_width = 1.1;
                                               
  
  const size_t nchannel = meas ? meas->num_gamma_channels() : 0;
  
  if( nchannel < min_high_res_chan )
    return false;
  
  if( nchannel > max_low_res_chan )
    return true;
  
  const auto cal = meas->energy_calibration();
  if( !cal || !cal->valid()
     || (cal->type() == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial) )
    return (nchannel > max_low_res_chan);
  
  const float energy_range = cal->upper_energy() - cal->lower_energy();
  const float avrg_channel_width = energy_range / nchannel;
  
  // Note: some (usually low resolution) systems have a sqrt(energy) scale, which invalidates
  //       the average energy assumptions, and havent been looked into yet.
  
  return (avrg_channel_width <= max_high_res_chann_width);
}//bool is_high_res(...)


}//namespace PeakFitUtils
