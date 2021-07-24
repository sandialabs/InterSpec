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
#include <iostream>
#include <stdexcept>

#include <boost/math/distributions/normal.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectionLimitCalc.h"

using namespace std;

namespace DetectionLimitCalc
{

CurieMdaInput::CurieMdaInput()
  : gamma_energy(0.0f), peak_region_nchannel(0.0f), num_lower_side_channels(0),
    num_upper_side_channels(0), detection_probability(0.0f), additional_uncertainty(0.0f)
{
}


CurieMdaResult::CurieMdaResult()
  : first_lower_continuum_channel(0), last_lower_continuum_channel(0), lower_continuum_counts_sum(0.0f),
    first_upper_continuum_channel(0), last_upper_continuum_channel(0), upper_continuum_counts_sum(0.0f),
    first_peak_region_channel(0), last_peak_region_channel(0), peak_region_counts_sum(0.0f),
    continuum_eqn{ 0.0f, 0.0f },
    estimated_peak_continuum_counts(0.0f), estimated_peak_continuum_uncert(0.0f),
    decision_threshold(0.0f), detection_limit(0.0f), source_counts(0.0f),
    lower_limit(0.0f), upper_limit(0.0f)
{
}


std::ostream &print_summary( std::ostream &strm, const CurieMdaResult &result, const float w )
{
  /*
   {
   //If peak looks present, use 2.5 FWHM for peak width
   //If no peak is located, then use 1.2 FWHM for peak width
   const float fwhm = detector->peakResolutionFWHM(result.input.gamma_energy);
   const shared_ptr<const SpecUtils::EnergyCalibration> cal = result.input.spectrum->energy_calibration();
   const float peak_channel = cal->channel_for_energy(result.input.gamma_energy);
   const float lower_energy = cal->energy_for_channel( peak_channel - 0.5*result.input.peak_region_nchannel );
   const float upper_energy = cal->energy_for_channel( peak_channel + 0.5*result.input.peak_region_nchannel );
   
   strm << "Hacky estimate of width is " << (upper_energy - lower_energy)/fwhm << " FWHMs" << endl;
   }
   */
  
  strm << "Activity less than ";
  if( w > 0.0 )
    strm << PhysicalUnits::printToBestActivityUnits(w*result.upper_limit) << " (" << result.upper_limit << " counts)" << endl;
  else
    strm << result.upper_limit << " counts" << endl;
  
  strm << "Activity greater than ";
  if( w > 0.0 )
    strm << PhysicalUnits::printToBestActivityUnits(w*result.lower_limit) << " (" << result.lower_limit << " counts)" << endl;
  else
    strm << result.lower_limit << " counts" << endl;
  
  strm << "Nominal activity estimate: ";
  if( w > 0.0 )
     strm << PhysicalUnits::printToBestActivityUnits(w*result.source_counts) <<" (" << result.source_counts << " counts)" << endl;
  else
    strm << result.source_counts << " counts" << endl;
  
  strm << "Continuum Starting Channel: " << result.first_lower_continuum_channel << endl;
  strm << "Continuum Ending Channel: " << result.last_upper_continuum_channel + 1 << endl;
  strm << "Peak Starting Channel: " << result.first_peak_region_channel << endl;
  strm << "Peak Ending Channel: " << result.last_peak_region_channel + 1 << endl;
  
  
  strm << "Critical Limit (decision threshold) in Peak Region: "<< result.decision_threshold << " counts";
  if( w > 0 )
    strm << " or " << PhysicalUnits::printToBestActivityUnits(w*result.decision_threshold) << endl;
  
  strm << "Detection limit: " << result.detection_limit << " counts";
  if( w > 0 )
    strm << ", or " << PhysicalUnits::printToBestActivityUnits(w*result.detection_limit) << endl;
  
  
  strm << "Gross Counts in Peak Foreground: " << result.source_counts << endl;
  //Gross Counts in Peak Background: 0
  strm << "Gross Counts in Continuum Foreground: " << result.estimated_peak_continuum_counts << endl;
  //Gross Counts in Continuum Background: 0
  strm << "Uncert in Peak Region: " << result.estimated_peak_continuum_uncert << endl;
  //Variance in Continuum Region: 14.17189
  //Detector Efficiency at Energy: 0.1304517
  //Range of Gammas Entering Detector: (0, 0.1656136)
  //Solid Angle: 0.0002390204
  //Range of Gammas from Source: (0, 692.8849)
  //Gamma Emission Rate per uCi of Source: 3.699356E+10
  //Transmission through Shielding: 1
  
  strm << "Continuum eqn (cnts/keV): " << result.continuum_eqn[0]
       << (result.continuum_eqn[1] > 0.0 ? " + " : " - ")
       << fabs(result.continuum_eqn[1]) << "*x"
       << ", with x=energy-" << result.input.gamma_energy << endl;
  
  return strm;
};//print_summary( CurieMdaResult )

CurieMdaResult currie_mda_calc( const CurieMdaInput &input )
{
  using namespace SpecUtils;
  
  const shared_ptr<const Measurement> spec = input.spectrum;
  const shared_ptr<const EnergyCalibration> cal = spec ? spec->energy_calibration() : nullptr;
  const size_t nchannel = cal ? cal->num_channels() : size_t(0);
  
  if( !nchannel || !cal->valid() || !spec->gamma_counts() || !cal->channel_energies() )
    throw std::runtime_error( "mda_counts_calc: invalid spectrum passed in." );
  
  if( input.gamma_energy <= 0.0f )
    throw std::runtime_error( "mda_counts_calc: invalid gamma_energy." );
  
  if( (input.peak_region_nchannel) <= 2.0f || (input.peak_region_nchannel >= nchannel) )
    throw std::runtime_error( "mda_counts_calc: invalid peak_region_nchannel." );
  
  if( input.num_lower_side_channels < 1 || (input.num_lower_side_channels >= nchannel)  )
    throw std::runtime_error( "mda_counts_calc: invalid num_lower_side_channels." );
  
  if( input.num_upper_side_channels < 1 || (input.num_upper_side_channels >= nchannel) )
    throw std::runtime_error( "mda_counts_calc: invalid num_upper_side_channels." );
  
  if( input.detection_probability <= 0.05f || input.detection_probability >= 1.0f )
    throw std::runtime_error( "mda_counts_calc: invalid detection_probability." );
  
  if( input.additional_uncertainty < 0.0f || input.additional_uncertainty >= 1.0f )
    throw std::runtime_error( "mda_counts_calc: invalid additional_uncertainty." );
  
  const vector<float> &gamma_counts = *spec->gamma_counts();
  const vector<float> &gamma_energies = *cal->channel_energies();
  
  assert( gamma_energies.size() == (gamma_counts.size() + 1) );
  
  const float mean_channel = cal->channel_for_energy( input.gamma_energy );
  const float peak_region_lower_ch = mean_channel - 0.5f*input.peak_region_nchannel;
  const float peak_region_upper_ch = mean_channel + 0.5f*input.peak_region_nchannel;
  
  CurieMdaResult result;
  result.input = input;
  
  result.first_peak_region_channel = static_cast<size_t>( std::round(peak_region_lower_ch) );
  
  // We need to makeup for the channel number defining the left side of each channel, so we will
  //  subtract off 0.5 from the channel we are supposed to go up through.
  result.last_peak_region_channel = static_cast<size_t>( std::round(peak_region_upper_ch - 0.5) );
  
  // TODO: need to verify the peak region is actually the correct number of channels, and do a better job of centering the requested number of channels around the input.gamma_energy
  
  if( result.first_peak_region_channel < (input.num_lower_side_channels + 1) )
    throw std::runtime_error( "mda_counts_calc: lower peak region is outside spectrum energy range" );
  
  result.last_lower_continuum_channel = result.first_peak_region_channel - 1;
  result.first_lower_continuum_channel = result.last_lower_continuum_channel - input.num_lower_side_channels + 1;
  
  result.first_upper_continuum_channel = result.last_peak_region_channel + 1;
  result.last_upper_continuum_channel = result.first_upper_continuum_channel + input.num_upper_side_channels - 1;
  
  if( result.last_upper_continuum_channel >= nchannel  )
    throw std::runtime_error( "mda_counts_calc: upper peak region is outside spectrum energy range" );
  
  
  result.lower_continuum_counts_sum = spec->gamma_channels_sum(result.first_lower_continuum_channel, result.last_lower_continuum_channel);
  result.peak_region_counts_sum = spec->gamma_channels_sum(result.first_peak_region_channel, result.last_peak_region_channel);
  result.upper_continuum_counts_sum = spec->gamma_channels_sum(result.first_upper_continuum_channel, result.last_upper_continuum_channel);
  
  /*
   cout << "Lower region:\n\tChan\tEne\tCounts" << endl;
   for( size_t i = result.first_lower_continuum_channel; i <= result.last_lower_continuum_channel; ++i )
   cout << "\t" << i << "\t" << gamma_energies[i] << "\t" << gamma_counts[i] << endl;
   cout << "\tSum: " << result.lower_continuum_counts_sum << endl;
   
   cout << "\nPeak region:\n\tChan\tEne\tCounts" << endl;
   for( size_t i = result.first_peak_region_channel; i <= result.last_peak_region_channel; ++i )
   cout << "\t" << i << "\t" << gamma_energies[i] << "\t" << gamma_counts[i] << endl;
   cout << "\tSum: " << result.peak_region_counts_sum << endl;
   
   cout << "\nUpper region:\n\tChan\tEne\tCounts" << endl;
   for( size_t i = result.first_upper_continuum_channel; i <= result.last_upper_continuum_channel; ++i )
   cout << "\t" << i << "\t" << gamma_energies[i] << "\t" << gamma_counts[i] << endl;
   cout << "\tSum: " << result.upper_continuum_counts_sum << endl;
   */
  
  const double lower_cont_counts = spec->gamma_channels_sum(result.first_lower_continuum_channel, result.last_lower_continuum_channel);
  const double upper_cont_counts = spec->gamma_channels_sum(result.first_upper_continuum_channel, result.last_upper_continuum_channel);
  const double lower_cont_width = spec->gamma_channel_upper(result.last_lower_continuum_channel)
  - spec->gamma_channel_lower(result.first_lower_continuum_channel);
  const double upper_cont_width = spec->gamma_channel_upper(result.last_upper_continuum_channel)
  - spec->gamma_channel_lower(result.first_upper_continuum_channel);
  
  const double lower_cont_density = lower_cont_counts / lower_cont_width;
  const double lower_cont_density_uncert = lower_cont_density / sqrt(lower_cont_counts);
  
  const double upper_cont_density = upper_cont_counts / upper_cont_width;
  const double upper_cont_density_uncert = upper_cont_density / sqrt(upper_cont_counts);
  
  const double peak_cont_density = 0.5*(lower_cont_density + upper_cont_density);
  const double peak_cont_density_uncert = 0.5*sqrt( upper_cont_density_uncert*upper_cont_density_uncert
                                                   + lower_cont_density_uncert*lower_cont_density_uncert );
  const double peak_cont_frac_uncert = peak_cont_density_uncert / peak_cont_density;
  
  
  const double peak_area_width = spec->gamma_channel_upper(result.last_peak_region_channel)
  - spec->gamma_channel_lower(result.first_peak_region_channel);
  const double peak_cont_sum = peak_cont_density * peak_area_width;
  const double peak_cont_sum_uncert = peak_cont_sum * peak_cont_frac_uncert;
  
  result.estimated_peak_continuum_counts = static_cast<float>( peak_cont_sum );
  result.estimated_peak_continuum_uncert = static_cast<float>( peak_cont_sum_uncert );
  
  
  // The equation is centered around the input.gamma_energy with the density of counts at normal
  //  value at that point.  The Slope will be through the midpoints of each continuum.
  // TODO: should do a proper least-squares fit to the continuum; I think this will give us a slightly true-er answer
  
  const double lower_cont_mid_energy = spec->gamma_channel_lower(result.first_lower_continuum_channel) + 0.5*lower_cont_width;
  const double upper_cont_mid_energy = spec->gamma_channel_lower(result.first_upper_continuum_channel) + 0.5*upper_cont_width;
  
  result.continuum_eqn[1] = (upper_cont_density - lower_cont_density) / (upper_cont_mid_energy - lower_cont_mid_energy);
  result.continuum_eqn[0] = lower_cont_density - result.continuum_eqn[1]*(lower_cont_mid_energy - input.gamma_energy);
  
  
  {// begin sanity check on continuum eqn
    const double peak_start_eq = spec->gamma_channel_lower(result.first_peak_region_channel) - input.gamma_energy;
    const double peak_end_eq = spec->gamma_channel_upper(result.last_peak_region_channel) - input.gamma_energy;
    
    const double peak_cont_eq_integral = result.continuum_eqn[0] * (peak_end_eq - peak_start_eq)
    + result.continuum_eqn[1] * 0.5 * (peak_end_eq*peak_end_eq - peak_start_eq*peak_start_eq);
    const double upper_cont_eq = result.continuum_eqn[0] + (upper_cont_mid_energy - input.gamma_energy)*result.continuum_eqn[1];
    
    assert( fabs(upper_cont_eq - upper_cont_density) < 0.1 ); //arbitrary precision test, for development
    assert( fabs(peak_cont_eq_integral - peak_cont_sum) < 0.1 ); //arbitrary precision test, for development
  }// end sanity check on continuum eqn
  
  
  typedef boost::math::policies::policy<boost::math::policies::digits10<6> > my_pol_6;
  const boost::math::normal_distribution<float,my_pol_6> gaus_dist( 0.0f, 1.0f );
  
  //  TODO: If/when we start having k_alpha != k_beta, then we probably need to be more careful
  //        around single vs double sided quantile.
  //   Will map 0.8414->1.00023, 0.95->1.64485, 0.975->1.95996, 0.995->2.57583, ...
  const float k = boost::math::quantile( gaus_dist, input.detection_probability );
  
  
  const double peak_cont_sigma = sqrt( peak_cont_sum_uncert*peak_cont_sum_uncert + peak_cont_sum );
  
  result.decision_threshold = k * peak_cont_sigma; //Note if using non-symmetric coverage, we would use k_alpha here
  
  // TODO: The calculation of result.detection_limit is using the simplified form requiring k_alpha == k_beta.
  //       If this is not the case it is an iterative solution (see eqn 129 on pg 47 of AQ-48)
  const double add_uncert = input.additional_uncertainty;
  if( k*k*add_uncert >= 1.0 )
  {
    result.detection_limit = -999; //TODO: indicate non-applicability better
  }else
  {
    // TODO: Using tbl 16 of AQ-48, I get a slightly high answer of 193.05 vs the expected answer of 191.906.
    //       If in the numerator I replace k*k with just k, then I get 191.983
    //       And if I plug those tables numbers into eqn 129 I get 0.474598 vs their 0.471705,
    //       so I am currently suspecting an error in the table.
    result.detection_limit = (2.0*result.decision_threshold + k*k) / (1.0 - k*k*add_uncert*add_uncert);
  }
  
  
  const float source_counts = result.peak_region_counts_sum - result.estimated_peak_continuum_counts;
  result.source_counts = source_counts;
  
  double region_sigma = peak_cont_sum_uncert*peak_cont_sum_uncert + result.peak_region_counts_sum;
  
  // TODO: I *think* this is right; e.g., use the nominal estimate of signal counts to estimate total uncertainty impact due to the "additional uncertainty" of the measurement, but I need to double check this.
  if( (source_counts > 0) && (add_uncert > 0) )
    region_sigma += source_counts*source_counts * add_uncert*add_uncert;
  region_sigma = sqrt( region_sigma );
  
  result.lower_limit = source_counts - k*region_sigma;
  result.upper_limit = source_counts + k*region_sigma;
  
  /*
   cout << "lower_cont_counts=" << lower_cont_counts << endl;
   cout << "upper_cont_counts=" << upper_cont_counts << endl;
   cout << "peak_cont_density=" << peak_cont_density << endl;
   cout << "peak_cont_density_uncert=" << peak_cont_density_uncert << endl;
   cout << "peak_cont_sum=" << peak_cont_sum << endl;
   cout << "peak_cont_sum_uncert=" << peak_cont_sum_uncert << endl;
   cout << "peak_cont_sigma=" << peak_cont_sigma << endl;
   cout << "result.decision_threshold=" << result.decision_threshold << endl;
   cout << "result.detection_limit=" << result.detection_limit << endl;
   cout << "result.source_counts=" << result.source_counts << endl;
   cout << "result.lower_limit=" << result.lower_limit << endl;
   cout << "result.upper_limit=" << result.upper_limit << endl;
   */
  
  return result;
};//mda_counts_calc



}//namespace DetectionLimitCalc

