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

#include "InterSpec/PeakFit.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace DetectionLimitCalc
{

CurieMdaInput::CurieMdaInput()
  : spectrum(nullptr),
    gamma_energy(0.0f), roi_lower_energy(0.0f), roi_upper_energy(0.0f),
    num_lower_side_channels(0), num_upper_side_channels(0),
    detection_probability(0.0f), additional_uncertainty(0.0f)
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


#if( PERFORM_DEVELOPER_CHECKS )

template<class T>
bool floats_equiv_enough( const T a, const T b )
{
  // This function checks the arguments are EITHER within 'abs_epsilon' of each other,
  //  OR 'rel_epsilon * max(a,b)' of each other.
  const T abs_epsilon = 1.0E-10f; //arbitrary, could pick like std::numeric_limits<float>::epsilon();
  const T rel_epsilon = 1.0E-5f;
  
  const auto diff = fabs(a - b);
  const auto maxval = std::max(fabs(a),fabs(b));
  
  return (diff < (rel_epsilon*maxval) || (diff < abs_epsilon));
};


void CurieMdaResult::equal_enough( const CurieMdaResult &test, const CurieMdaResult &expected )
{
  //CurieMdaInput::equal_enough( test.input, expected.input );
  
  vector<string> errors;
  char buffer[512];
  
  if( test.first_lower_continuum_channel != expected.first_lower_continuum_channel )
    errors.push_back( "Test first_lower_continuum_channel ("
                     + std::to_string(test.first_lower_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.first_lower_continuum_channel) );
  
  if( test.last_lower_continuum_channel != expected.last_lower_continuum_channel )
    errors.push_back( "Test last_lower_continuum_channel ("
                     + std::to_string(test.last_lower_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.last_lower_continuum_channel) );
  
  if( !floats_equiv_enough(test.lower_continuum_counts_sum, expected.lower_continuum_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
              "Test lower_continuum_counts_sum (%.6G) does not equal expected (%.6G)",
             test.lower_continuum_counts_sum, expected.lower_continuum_counts_sum );
    errors.push_back( buffer );
  }
  
  if( test.first_upper_continuum_channel != expected.first_upper_continuum_channel )
    errors.push_back( "Test first_upper_continuum_channel ("
                     + std::to_string(test.first_upper_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.first_upper_continuum_channel) );
  
  
  if( test.last_upper_continuum_channel != expected.last_upper_continuum_channel )
    errors.push_back( "Test last_upper_continuum_channel ("
                     + std::to_string(test.last_upper_continuum_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.last_upper_continuum_channel) );
  
  if( !floats_equiv_enough(test.upper_continuum_counts_sum, expected.upper_continuum_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test lower_continuum_counts_sum (%.6G) does not equal expected (%.6G)",
             test.upper_continuum_counts_sum, expected.upper_continuum_counts_sum );
    errors.push_back( buffer );
  }
  
  
  if( test.first_peak_region_channel != expected.first_peak_region_channel )
    errors.push_back( "Test first_peak_region_channel ("
                     + std::to_string(test.first_peak_region_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.first_peak_region_channel) );
  
  if( test.last_peak_region_channel != expected.last_peak_region_channel )
    errors.push_back( "Test last_peak_region_channel ("
                     + std::to_string(test.last_peak_region_channel)
                     + ") does not equal expected ("
                     + std::to_string(expected.last_peak_region_channel) );


  if( !floats_equiv_enough(test.peak_region_counts_sum, expected.peak_region_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test peak_region_counts_sum (%.6G) does not equal expected (%.6G)",
             test.peak_region_counts_sum, expected.peak_region_counts_sum );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.peak_region_counts_sum, expected.peak_region_counts_sum) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test peak_region_counts_sum (%.6G) does not equal expected (%.6G)",
             test.peak_region_counts_sum, expected.peak_region_counts_sum );
    errors.push_back( buffer );
  }
  
  /*
  for( int i = 0; i < 2; ++i )
  {
    if( !floats_equiv_enough(test.continuum_eqn[i], expected.continuum_eqn[i]) )
    {
      snprintf( buffer, sizeof(buffer),
               "Test continuum_eqn[%i] (%.6G) does not equal expected (%.6G)",
               i, test.continuum_eqn[i], expected.continuum_eqn[i] );
      errors.push_back( buffer );
    }
  }//for( int i = 0; i < 2; ++i )
  */
  
  if( !floats_equiv_enough(test.estimated_peak_continuum_counts, expected.estimated_peak_continuum_counts) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test estimated_peak_continuum_counts (%.6G) does not equal expected (%.6G)",
             test.estimated_peak_continuum_counts, expected.estimated_peak_continuum_counts );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.estimated_peak_continuum_uncert, expected.estimated_peak_continuum_uncert) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test estimated_peak_continuum_uncert (%.6G) does not equal expected (%.6G)",
             test.estimated_peak_continuum_uncert, expected.estimated_peak_continuum_uncert );
    errors.push_back( buffer );
  }
  
  

  if( !floats_equiv_enough(test.decision_threshold, expected.decision_threshold) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test decision_threshold (%.6G) does not equal expected (%.6G)",
             test.decision_threshold, expected.decision_threshold );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.detection_limit, expected.detection_limit) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test detection_limit (%.6G) does not equal expected (%.6G)",
             test.detection_limit, expected.detection_limit );
    errors.push_back( buffer );
  }
  

  if( !floats_equiv_enough(test.source_counts, expected.source_counts) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test source_counts (%.6G) does not equal expected (%.6G)",
             test.source_counts, expected.source_counts );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.lower_limit, expected.lower_limit) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test lower_limit (%.6G) does not equal expected (%.6G)",
             test.lower_limit, expected.lower_limit );
    errors.push_back( buffer );
  }
  
  
  if( !floats_equiv_enough(test.upper_limit, expected.upper_limit) )
  {
    snprintf( buffer, sizeof(buffer),
             "Test upper_limit (%.6G) does not equal expected (%.6G)",
             test.upper_limit, expected.upper_limit );
    errors.push_back( buffer );
  }
  
  if( errors.empty() )
    return;
  
  string err_msg = "CurieMdaResult::equal_enough: test and expected values are not equal\n";
  for( size_t i = 0; i < errors.size(); ++i )
    err_msg += (i ? "\n\t" : "\t") + errors[i];
    
  throw runtime_error( err_msg );
}//CurieMdaResult::equal_enough(...)
#endif //PERFORM_DEVELOPER_CHECKS



std::ostream &print_summary( std::ostream &strm, const CurieMdaResult &result, const float w )
{
  /*
   {
   //If peak looks present, use 2.5 FWHM for peak width
   //If no peak is located, then use 1.2 FWHM for peak width
   const float fwhm = detector->peakResolutionFWHM(result.input.gamma_energy);
   const shared_ptr<const SpecUtils::EnergyCalibration> cal = result.input.spectrum->energy_calibration();
   const float peak_channel = cal->channel_for_energy(result.input.gamma_energy);
   const float lower_energy = result.input.roi_lower_energy;
   const float upper_energy = result.input.roi_upper_energy;
   
   strm << "Width is " << (upper_energy - lower_energy)/fwhm << " FWHMs" << endl;
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
    strm << " or " << PhysicalUnits::printToBestActivityUnits(w*result.decision_threshold);
  
  strm << "\nDetection limit: " << result.detection_limit << " counts";
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
    throw runtime_error( "mda_counts_calc: invalid spectrum passed in." );
  
  if( IsInf(input.roi_lower_energy) || IsNan(input.roi_lower_energy)
     ||  IsInf(input.roi_upper_energy) || IsNan(input.roi_upper_energy) )
    throw runtime_error( "mda_counts_calc: invalid ROI specified." );
  
  if( input.roi_lower_energy >= input.roi_upper_energy )
    throw runtime_error( "mda_counts_calc: upper ROI energy must be greater than lower energy." );
  
  if( (input.gamma_energy < input.roi_lower_energy)
      || (input.gamma_energy > input.roi_upper_energy) )
    throw runtime_error( "mda_counts_calc: gamma energy must be between lower and upper ROI." );
  
  if( input.num_lower_side_channels < 1 || (input.num_lower_side_channels >= nchannel)  )
    throw runtime_error( "mda_counts_calc: invalid num_lower_side_channels." );
  
  if( input.num_upper_side_channels < 1 || (input.num_upper_side_channels >= nchannel) )
    throw runtime_error( "mda_counts_calc: invalid num_upper_side_channels." );
  
  if( input.detection_probability <= 0.05f || input.detection_probability >= 1.0f )
    throw runtime_error( "mda_counts_calc: invalid detection_probability." );
  
  if( input.additional_uncertainty < 0.0f || input.additional_uncertainty >= 1.0f )
    throw runtime_error( "mda_counts_calc: invalid additional_uncertainty." );
  
  const vector<float> &gamma_counts = *spec->gamma_counts();
  const vector<float> &gamma_energies = *cal->channel_energies();
  
  assert( gamma_energies.size() == (gamma_counts.size() + 1) );
  
  const float peak_region_lower_ch = cal->channel_for_energy( input.roi_lower_energy );
  const float peak_region_upper_ch = cal->channel_for_energy( input.roi_upper_energy );
  
  if( (peak_region_lower_ch - input.num_lower_side_channels) < 0.0 )
    throw runtime_error( "mda_counts_calc: lower energy goes off spectrum" );
  
  if( (peak_region_upper_ch + input.num_upper_side_channels) > cal->num_channels() )
    throw runtime_error( "mda_counts_calc: upper energy goes off spectrum" );
  
  CurieMdaResult result;
  result.input = input;
  
  
  result.first_peak_region_channel = static_cast<size_t>( std::round(peak_region_lower_ch) );
  
  // If we pass in exactly the channel boundary, or really close to it, we want to round in the
  //  reasonable way, otherwise we need to makeup for the channel number defining the left side of
  //  each channel, so we will subtract off 0.5 from the channel we are supposed to go up through.
  if( fabs(peak_region_upper_ch - std::floor(peak_region_upper_ch)) < 0.01 )
    result.last_peak_region_channel = static_cast<size_t>( std::floor(peak_region_upper_ch) - 1 );
  else
    result.last_peak_region_channel = static_cast<size_t>( std::round(peak_region_upper_ch - 0.5) );
  
  
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
  const double lower_cont_density_uncert = ((lower_cont_counts <= 0.0) ? 0.0 : (lower_cont_density / sqrt(lower_cont_counts)));
  
  const double upper_cont_density = upper_cont_counts / upper_cont_width;
  const double upper_cont_density_uncert = ((upper_cont_counts <= 0.0) ? 0.0 : (upper_cont_density / sqrt(upper_cont_counts)));
  
  const double peak_cont_density = 0.5*(lower_cont_density + upper_cont_density);
  const double peak_cont_density_uncert = 0.5*sqrt( upper_cont_density_uncert*upper_cont_density_uncert
                                                   + lower_cont_density_uncert*lower_cont_density_uncert );
  const double peak_cont_frac_uncert = ((peak_cont_density > 0.0) ? (peak_cont_density_uncert / peak_cont_density) : 1.0);
  
  
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
  
#if( PERFORM_DEVELOPER_CHECKS )
  {// begin sanity check on continuum eqn
    const double peak_start_eq = spec->gamma_channel_lower(result.first_peak_region_channel) - input.gamma_energy;
    const double peak_end_eq = spec->gamma_channel_upper(result.last_peak_region_channel) - input.gamma_energy;
    
    const double peak_cont_eq_integral = result.continuum_eqn[0] * (peak_end_eq - peak_start_eq)
          + result.continuum_eqn[1] * 0.5 * (peak_end_eq*peak_end_eq - peak_start_eq*peak_start_eq);
    const double upper_cont_eq = result.continuum_eqn[0] + (upper_cont_mid_energy - input.gamma_energy)*result.continuum_eqn[1];
    
    // Arbitrary chosen precision tests, for development
    const double eq_dens = fabs(peak_cont_eq_integral - peak_cont_sum);
    assert( (eq_dens < 0.1)
           || (eq_dens < 1.0E-5*std::max(peak_cont_eq_integral, peak_cont_sum)) );
    
    const double eq_diff = fabs(peak_cont_eq_integral - peak_cont_sum);
    assert( eq_diff < 0.1 || eq_diff < 1.0E-5*std::max(peak_cont_eq_integral, peak_cont_sum) );
  }// end sanity check on continuum eqn
#endif //PERFORM_DEVELOPER_CHECKS
  
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
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  assert( !IsInf(result.lower_limit) && !IsNan(result.lower_limit) );
  assert( !IsInf(result.upper_limit) && !IsNan(result.upper_limit) );
#endif
  
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


DeconRoiInfo::DeconRoiInfo()
: roi_start( 0.0f ),
  roi_end( 0.0f ),
  continuum_type( PeakContinuum::OffsetType::NoOffset ),
  fix_continuum_to_edges( false ),
  num_lower_side_channels( 0 ),
  num_upper_side_channels( 0 ),
  peak_infos()
{
}


DeconRoiInfo::PeakInfo::PeakInfo()
 : energy( 0.0f ),
   fwhm( 0.0f ),
   counts_per_bq_into_4pi( 0.0 )
{
}



DeconComputeInput::DeconComputeInput()
 : distance( 0.0 ),
   activity( 0.0 ),
   include_air_attenuation( false ),
   shielding_thickness( 0.0 ),
   drf( nullptr ),
   measurement( nullptr ),
   roi_info()
{
}


DeconComputeResults decon_compute_peaks( const DeconComputeInput &input )
{
  DeconComputeResults result;
  result.input = input;
  result.chi2 = 0.0;
  result.num_degree_of_freedom = 0;
  
  // Lets sanity check input
  if( (input.distance <= 0.0) || IsNan(input.distance) || IsInf(input.distance) )
    throw runtime_error( "decon_compute_peaks: invalid input distance" );
  
  // Lets sanity check input
  if( (input.activity <= 0.0) || IsNan(input.activity) || IsInf(input.activity) )
    throw runtime_error( "decon_compute_peaks: invalid input activity" );
  
  if( input.include_air_attenuation
      && ((input.shielding_thickness < 0.0)
          || IsNan(input.shielding_thickness)
          || IsInf(input.shielding_thickness)
          || (input.shielding_thickness >= input.distance)) )
    throw runtime_error( "decon_compute_peaks: invalid input shielding thickness" );
  
  if( !input.drf || !input.drf->isValid() || !input.drf->hasResolutionInfo() )
    throw runtime_error( "decon_compute_peaks: invalid DRF input" );
  
  
  if( !input.measurement || (input.measurement->num_gamma_channels() < 2)
     || !input.measurement->energy_calibration() //ptr should always be valid anyway
     || !input.measurement->energy_calibration()->valid() )
    throw runtime_error( "decon_compute_peaks: invalid spectrum input" );
  
  // Check if there is anything to do
  if( input.roi_info.empty() )
    return result;
    
  // We should be good to go,
  vector<PeakDef> inputPeaks, fitPeaks;
  
  for( const DeconRoiInfo &roi : input.roi_info )
  {
    const float  &roi_start = roi.roi_start; //This _should_ already be rounded to nearest bin edge; TODO: check that this is rounded
    const float  &roi_end = roi.roi_end; //This _should_ already be rounded to nearest bin edge; TODO: check that this is rounded
    const bool   &fix_continuum = roi.fix_continuum_to_edges;
    const PeakContinuum::OffsetType &continuum_type = roi.continuum_type;
    const size_t &num_lower_side_channels = roi.num_lower_side_channels;
    const size_t &num_upper_side_channels = roi.num_upper_side_channels;
    
    
    shared_ptr<PeakContinuum> peak_continuum;
    std::shared_ptr<const SpecUtils::Measurement> computed_global_cont;
    
    // Find the largest peak in the ROIs, energy to use as the continuum "reference energy"
    float reference_energy = 0.0f;
    for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
      reference_energy = std::max( reference_energy, peak_info.energy );
    
    assert( reference_energy != 0.0f );
    
    for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
    {
      const float &energy = peak_info.energy;
      const float fwhm = input.drf->peakResolutionFWHM( energy );
      const float sigma = fwhm / 2.634;
      const double det_eff = input.drf->efficiency( energy, input.distance );
      const double counts_4pi = peak_info.counts_per_bq_into_4pi;
      double air_atten = 1.0;
      
      if( input.include_air_attenuation )
      {
        const double air_len = input.distance - input.shielding_thickness;
        const double mu_air = GammaInteractionCalc::transmission_coefficient_air( energy, air_len );
        air_atten = exp( -1.0 * mu_air );
      }
      
      const float amplitude = input.activity * counts_4pi * det_eff * air_atten;
    
      PeakDef peak( energy, sigma, amplitude );
      peak.setFitFor( PeakDef::CoefficientType::Mean, false );
      peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
      peak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
      
      
      if( peak_continuum )
      {
        peak.setContinuum( peak_continuum );
      }else
      {
        peak_continuum = peak.continuum();
        peak_continuum->setType( continuum_type );
        peak_continuum->setRange( roi_start, roi_end );
        
        // First, we'll find a linear continuum as the starting point, and then go through
        //  and modify it how we need
        peak_continuum->calc_linear_continuum_eqn( input.measurement, reference_energy,
                                                  roi_start, roi_end,
                                                  num_lower_side_channels,
                                                  num_upper_side_channels );
        
        for( size_t order = 0; order < 2; ++order )  //peak_continuum->parameters().size()
          peak_continuum->setPolynomialCoefFitFor( order, !fix_continuum );
        
        switch( continuum_type )
        {
          case PeakContinuum::NoOffset:
            break;
            
          case PeakContinuum::Constant:
          
            
            break;
            
          case PeakContinuum::Linear:
          case PeakContinuum::FlatStep:
          {
            
            
            
            
            
            
            break;
          }//
            
          case PeakContinuum::Quadratic:
          case PeakContinuum::Cubic:
          case PeakContinuum::LinearStep:
          case PeakContinuum::BiLinearStep:
            break;
            
            
          case PeakContinuum::External:
            if( !computed_global_cont )
              computed_global_cont = estimateContinuum( input.measurement );
            peak_continuum->setExternalContinuum( computed_global_cont );
            break;
        }//switch( assign continuum )
      }//if( peak_continuum ) / else
      
      inputPeaks.push_back( std::move(peak) );
    }//for( const DeconRoiInfo::PeakInfo &peak_info : roi.peak_infos )
  }//for( size_t i = 0; i < energies.size(); ++i )
  
  if( inputPeaks.empty() )
    throw runtime_error( "decon_compute_peaks: No peaks given in ROI(s)" );
    
#warning "Totally not computing things yet"
  throw runtime_error( "decon_compute_peaks: computations not implemented" );
  
  /*
  ROOT::Minuit2::MnUserParameters inputFitPars;
  PeakFitChi2Fcn::addPeaksToFitter( inputFitPars, inputPeaks, spec, PeakFitChi2Fcn::kFitForPeakParameters );
  
  const int npeaks = static_cast<int>( inputPeaks.size() );
  PeakFitChi2Fcn chi2Fcn( npeaks, spec, nullptr );
  chi2Fcn.useReducedChi2( false );
  
  assert( inputFitPars.VariableParameters() != 0 );
  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputFitPars );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2Fcn, inputParamState, strategy );
  
  unsigned int maxFcnCall = 5000;
  double tolerance = 2.5;
  tolerance = 0.5;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  const ROOT::Minuit2::MnUserParameters &fitParams = fitter.Parameters();
  //  minimum.IsValid()
  //      ROOT::Minuit2::MinimumState minState = minimum.State();
  //      ROOT::Minuit2::MinimumParameters minParams = minState.Parameters();
  
  //    cerr << endl << endl << "EDM=" << minimum.Edm() << endl;
  //    cerr << "MinValue=" <<  minimum.Fval() << endl << endl;
  
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  
  if( !minimum.IsValid() )
  {
    //XXX - should we try to re-fit here? Or do something to handle the
    //      faliure in some reasonable way?
    cerr << endl << endl << "status is not valid"
    << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
    << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
    << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
    << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
    << "\n\tHasValidParameters: " << minimum.HasValidParameters()
    << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
    << endl;
    if( minimum.IsAboveMaxEdm() )
      cout << "\t\tEDM=" << minimum.Edm() << endl;
  }//if( !minimum.IsValid() )
  
  
  vector<double> fitpars = fitParams.Params();
  vector<double> fiterrors = fitParams.Errors();
  chi2Fcn.parametersToPeaks( fitPeaks, &fitpars[0], &fiterrors[0] );
  
  double initialChi2 = chi2Fcn.chi2( &fitpars[0] );
  
  //Lets try to keep whether or not to fit parameters should be the same for
  //  the output peaks as the input peaks.
  //Note that this doesnt account for peaks swapping with eachother in the fit
  assert( fitPeaks.size() == inputPeaks.size() );
  
  //for( size_t i = 0; i < near_peaks.size(); ++i )
  //  fitpeaks[i].inheritUserSelectedOptions( near_peaks[i], true );
  //for( size_t i = 0; i < fixedpeaks.size(); ++i )
  //  fitpeaks[i+near_peaks.size()].inheritUserSelectedOptions( fixedpeaks[i], true );
  
  const double totalNDF = set_chi2_dof( spec, fitPeaks, 0, fitPeaks.size() );
  
  chi2 = initialChi2;
  numDOF = static_cast<int>( std::round(totalNDF) );
  
  for( auto &peak : fitPeaks )
  {
    peak.setFitFor( PeakDef::CoefficientType::Mean, false );
    peak.setFitFor( PeakDef::CoefficientType::Sigma, false );
  }
  
  peaks.swap( fitPeaks );
   */
}//DeconComputeResults decon_compute_peaks( const DeconComputeInput &input )


}//namespace DetectionLimitCalc

