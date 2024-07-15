//  Created by wcjohns on 20110322.
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

/** This file is for developing peak-fit improving code. */

#include "InterSpec_config.h"

#include <map>
#include <deque>
#include <mutex>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/InterSpecServer.h"

using namespace std;

struct FindCandidateSettings
{
  int num_smooth_side_channels = 4; // low res more
  int smooth_polynomial_order = 3;  // highres 3, lowres 2
  double threshold_FOM = 1.3;
  float pos_sum_threshold_sf = -0.01f;
};//struct FindCandidateSettings

bool debug_printout = false;

void find_candidate_peaks( const std::shared_ptr<const SpecUtils::Measurement> data,
                          size_t start_channel,
                          size_t end_channel,
                          std::vector< std::tuple<float,float,float> > &results,
                          const FindCandidateSettings &settings )
{
  results.clear();
  
  if( !data->num_gamma_channels() )
    return;
  
  const size_t nchannel = data->num_gamma_channels();
  const bool highres = PeakFitUtils::is_high_res( data );
  
  if( start_channel >= nchannel )
    start_channel = 0;
  if( end_channel <= start_channel || end_channel >= (nchannel-1) )
    end_channel = nchannel - 2;
  
  const double threshold_FOM = settings.threshold_FOM;
  const float pos_sum_threshold_sf = settings.pos_sum_threshold_sf;
  
  //We will let one bin fluctate negative to avoid fluctations near threshold_FOM
  //Untested for HPGe data.
  //  Currently (20141209) for low res data, this should be kept the same as
  //  in find_roi_for_2nd_deriv_candidate(...) {although all this is under
  //  development}.
  const size_t nFluxuate = highres ? 2 : 2;
  assert( nFluxuate >= 1 );
  
  
  //XXX: using middle energy range so peak finding will always be consistent
  //     despite intput range
  //     Should keep consistent with find_roi_for_2nd_deriv_candidate()
  const size_t midbin = data->num_gamma_channels() / 2;// (start_channel + end_channel) / 2;
  const float midenergy = data->gamma_channel_center( midbin );
  const float midbinwidth = data->gamma_channel_width( midbin );
  
  const int order = settings.smooth_polynomial_order; //highres ? 3 : 2;
  const size_t side_bins = settings.num_smooth_side_channels; //highres ? 4 : std::max( size_t(5), static_cast<size_t>( 0.022f*midenergy/midbinwidth + 0.5f ) );
  
  
  vector<float> second_deriv;
  smoothSpectrum( data, static_cast<int>(side_bins), order, 2, second_deriv );
  
  //Since the peak resolution changes a lot for low-resolution spectra, if
  //  side_bins is to large, it will wipe out low energy features, meaning we
  //  wont detect low energy peaks
  //Note this code should be kept the same as in
  //  find_roi_for_2nd_deriv_candidate(...) while all of this is in development
  if( !highres && side_bins > 5 && nchannel >= 512
     && start_channel < (nchannel/15) )
  {
    //We should also have a minbimum statistics requirment here.
    const size_t index = std::max( (nchannel/15), side_bins );
    vector<float> second_deriv_lower;
    smoothSpectrum( data, 4, order, 2, second_deriv_lower );
    
    for( size_t i = 0; i < (index-side_bins); ++i )
      second_deriv[i] = second_deriv_lower[i];
    
    //transition over 'side_bins' between the smoothings.
    for( size_t i = 0; i < side_bins; ++i )
    {
      const float factor = float(i+1) / float(side_bins+1);
      const size_t current = index - side_bins + i;
      second_deriv[current] = factor*second_deriv[current]
      + (1.0f-factor)*second_deriv_lower[current];
    }
    
  }//if( !highres )
  
  //XXX: the below 1.5 is empiracally found, I'm not entirely sure where
  //     comes from...and infact might be higher
  const double amp_fake_factor = 1.5;
  
  
  const vector<float> &energies = *data->gamma_channel_energies();
  
  size_t minbin = 0, firstzero = 0, secondzero = 0;
  float secondsum = 0.0f, minval = 9999999999.9f;
  
  for( size_t channel = start_channel; channel <= end_channel; ++channel )
  {
    const float secondDeriv = second_deriv[channel]; //Not dividing by binwidth^2 here,
    
    bool secondSumPositive = true;
    float positivesum = 0.0f;
    for( size_t i = 0; i < nFluxuate; ++i )
    {
      if( (channel+i) <= end_channel )
      {
        const bool above = (second_deriv[channel+i] > 0.0f);
        if( above )
          positivesum += second_deriv[channel+i];
        secondSumPositive &= above;
      }
    }
    
    //Rather than using 'pos_sum_threshold_sf*secondsum', it might be better to
    //  use something invlolving sqrt(secondsum) since this makes a bit more
    //  sense.
    //Also, positivesum can probably also thresholded off of some sort of more
    //  absolute quantity
    secondSumPositive &= (positivesum > pos_sum_threshold_sf*secondsum);
    
    if( secondSumPositive && (minval < 99999999999.9)
       && (secondsum!=0.0) && (firstzero>0)
       && ((channel-firstzero)>2) )
    {
      secondzero = channel;
      const double mean = data->gamma_channel_center(minbin);
      const double sigma = 0.5*(data->gamma_channel_center(secondzero)
                                - data->gamma_channel_center(firstzero));
      
      if( debug_printout )
        cout << "secondzero=" << channel << ", energy=" << data->gamma_channel_center(channel) << ", sigma=" << sigma << endl;
      
      
      const double deriv_sigma = 0.5*( secondzero - firstzero );
      const double part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                       boost::math::constants::e<double>() ) )
      / ( deriv_sigma * deriv_sigma );
      const double amplitude = -amp_fake_factor * secondsum / part;
      
      
      const float lowerEnengy = static_cast<float>( mean - 3.0*sigma );
      const float upperEnergy = static_cast<float>( mean + 3.0*sigma );
      
      double data_area = data->gamma_integral( lowerEnengy, upperEnergy );
      double est_sigma = sqrt( std::max(data_area,1.0) );
      
      //In principle we would want to use the (true) continuums area to derive
      //  the est_sigma from, but for practical purposes I think this can give
      //  us false positives fairly often
      
      const double figure_of_merit = 0.68*amplitude / est_sigma;
      
      if( figure_of_merit > threshold_FOM )
      {
        bool passescuts = true;
        if( !highres && energies[minbin] < 130.0f )
        {
          //look 2 sigma forward and make sure the data has dropped enough
          const size_t p2sigmbin = data->find_gamma_channel( mean + 1.5*sigma );
          const float p2binlower = data->gamma_channel_lower( p2sigmbin );
          const float p2binupper = data->gamma_channel_upper( p2sigmbin );
          
          const double p2gausheight = amplitude*PeakDists::gaussian_integral( mean, sigma, p2binlower, p2binupper );
          const float p2contents = data->gamma_channel_content( p2sigmbin );
          
          const float meanbinlower = data->gamma_channel_lower( minbin );
          const float meanbinupper = data->gamma_channel_upper( minbin );
          const double meangausheight = amplitude*PeakDists::gaussian_integral( mean, sigma, meanbinlower, meanbinupper );
          const float meancontents = data->gamma_channel_content( minbin );
          const double expecteddiff = meangausheight - p2gausheight;
          const float actualdiff = meancontents - p2contents;
          
          //We dont expect the continuum to be changing rate of any more than
          //  1/2 the rate of peaks (this rate made up from my head, so take
          //  it with a grain of salt).  Should incorporate stat uncertainites
          //  here as well!
          passescuts = (actualdiff > 0.45*expecteddiff);
        }//if( !highres && energies[minbin] < 130.0f )
        
        if( !highres && passescuts )
        {
          //For really wide peaks, we'll let there be a bit more fluxuation,
          //  since sometimes
          size_t nflux = max( nFluxuate, ((secondzero-firstzero)/side_bins) + 1 );
          
          //Look forward and backwards to see if the sum of the second
          //  derivative, while positive, is roughly comparable to the negative
          //  sum.  We expect the positive-summs on either side of the negative
          //  region to add up to about the same area as the negative region.
          float nextpositivesum = 0.0;
          for( size_t i = channel; i <= end_channel; ++i )
          {
            bool secondSumNegative = true;
            for( size_t j = 0; j < nflux && ((i+j) < end_channel); ++j )
              secondSumNegative &= (second_deriv[i+j] < 0.0f);
            if( secondSumNegative )
              break;
            nextpositivesum += second_deriv[i];
          }//for( size_t i = 0; i <= end_channel; ++i )
          
          
          float prevpositivesum = 0.0;
          for( size_t i = firstzero; i > 0; --i )
          {
            bool secondSumNegative = true;
            for( size_t j = 0; j < nflux && ((i-j) > 0); ++j )
              secondSumNegative &= (second_deriv[i-j] < 0.0f);
            if( secondSumNegative )
              break;
            prevpositivesum += second_deriv[i];
          }//for( size_t i = 0; i <= end_channel; ++i )
          
          //If the current candidate peak is a result of compton backscatter,
          //  then it is likely the next positive region will be large, while
          //  the previous positive region will be small (we expect them to be
          //  about equal, and sum to be equal in magnitide to the negative
          //  region).
          const float nextratio = (-nextpositivesum/secondsum);
          const float prevratio = (-prevpositivesum/secondsum);
          passescuts = (nextratio < 4.0 || prevratio > 0.2)
          && ((nextratio+prevratio)>0.3);
          
        }//if( !highres )
        
        if( passescuts )
          results.push_back( std::tuple<float,float,float>{mean, sigma, amplitude} );
      }//if( region we were just in passed_threshold )
      
      secondsum = 0.0;
      minval = 9999999999.9f;
      minbin = secondzero = firstzero = 0;
    }else
    {
      bool belowzero = true, goingnegative = true, abovezero = true;
      for( size_t i = 0; i < nFluxuate; ++i )
      {
        if( (channel+i+1) < nchannel )
          goingnegative &= (second_deriv[channel+i+1] < 0.0f);
        if( channel >= i )
        {
          belowzero &= (second_deriv[channel-i] <= 0.0f);
          abovezero &= (second_deriv[channel-i] > 0.0f);
        }
      }//for( size_t i = 0; i < nFluxuate; ++i )
      
      if( channel && !firstzero && goingnegative )
      {
        if( debug_printout )
          cout << "firstzero=" << channel << ", energy=" << data->gamma_channel_center(channel) << endl;
        firstzero = channel;
        minbin = channel;
        minval = secondDeriv;
        
        for( size_t i = 1; i < nFluxuate; ++i )
          if( channel >= i )
            secondsum += second_deriv[channel-i];
      }else if( secondSumPositive )
      {
        secondsum = 0.0;
        minval = 9999999999.9f;
        minbin = secondzero = firstzero = 0;
      }
      
      if( firstzero > 0 )
      {
        secondsum += secondDeriv;
        
        if( secondDeriv < minval )
        {
          minbin = channel;
          minval = secondDeriv;
        }
      }//if( firstzero > 0 )
    }//if( we are out of region of interest) / else( in region of )
  }//for( loop over bins )
}//find_candidate_peaks(...)






struct PeakTruthInfo
{
  float energy;
  float cps;
  float area;
  float fwhm;
  
  /**
   Where:
   const double SC = 2.0 + 0.02*(lowSkew + highSkew)*std::pow(energy/661.0, resolutionPower );
   FullWidth = SC * GetFWHM(energy, resolutionOffset, resolution661, resolutionPower );
   */
  float full_width;
  
  /**
   "S.E.", "D.E.", "EscapeXRay", "Peak"
   */
  std::string label;
  
  PeakTruthInfo()
  : energy( 0.0 ), cps( 0.0 ), area( 0.0 ), fwhm( 0.0 ), full_width( 0.0 ), label()
  {
  }
  
  explicit PeakTruthInfo( const std::string &line )
  {
    vector<string> fields;
    SpecUtils::split( fields, line, "," );
    assert( fields.size() == 6 );
    if( fields.size() != 6 )
      throw runtime_error( "Unexpected number of fileds in line: '" + line + "'" );
    
    //"Energy (keV),CPS,Area,FWHM,FullWidth,Label"
    if( !SpecUtils::parse_float( fields[0].c_str(), fields[0].size(), energy ) )
      throw runtime_error( "Failed to parse energy of line: '" + line + "'" );
    
    if( !SpecUtils::parse_float( fields[1].c_str(), fields[1].size(), cps ) )
      throw runtime_error( "Failed to parse cps of line: '" + line + "'" );
    
    if( !SpecUtils::parse_float( fields[2].c_str(), fields[2].size(), area ) )
      throw runtime_error( "Failed to parse area of line: '" + line + "'" );
    
    if( !SpecUtils::parse_float( fields[3].c_str(), fields[3].size(), fwhm ) )
      throw runtime_error( "Failed to parse fwhm of line: '" + line + "'" );
    
    if( !SpecUtils::parse_float( fields[4].c_str(), fields[4].size(), full_width ) )
      throw runtime_error( "Failed to parse full_width of line: '" + line + "'" );
    
    label = fields[5];
    assert( (label == "S.E.") || (label == "D.E.") 
           || (label == "EscapeXRay") || (label == "Peak")
           || (label == "EscapeXRay+Peak") || (label == "D.E.+EscapeXRay")
           || (label == "D.E.+Peak") || (label == "D.E.+S.E.")
           || (label == "Peak+S.E.") );
  }//PeakTruthInfo(...)
};//struct PeakTruthInfo


struct InjectSourceInfo
{
  string file_base_path;
  string src_name;
  
  vector<PeakTruthInfo> source_lines;
  vector<PeakTruthInfo> background_lines;
  
  /** The PCF file for this source. */
  shared_ptr<SpecUtils::SpecFile> spec_file;
  
  /** The source + background spectra, that are Poisson varied. */
  vector<shared_ptr<const SpecUtils::Measurement>> src_spectra;
  
  /** Background that is Poisson varied, and same duration as src spectra. */
  shared_ptr<const SpecUtils::Measurement> short_background;
  
  /** A longer background, with Poisson variation. */
  shared_ptr<const SpecUtils::Measurement> long_background;
  
  /** Source + background without Poisson variation. */
  shared_ptr<const SpecUtils::Measurement> src_no_poisson;
  
  /** Background, with no Poisson variation. */
  shared_ptr<const SpecUtils::Measurement> background_no_poisson;
};//struct InjectSourceInfo


struct DetectorInjectSet
{
  string detector_name;
  string location_name;
  string live_time_name;
  
  vector<InjectSourceInfo> source_infos;
};//struct DetectorInjectSet



struct ExpectedPhotopeakInfo
{
  /** An approximate detection sigma, i.e., 2.33 times this value is kinds 95% DL */
  double nsigma_over_background;
  
  double roi_lower;
  double roi_upper;
  double effective_energy;
  double peak_area;
  double continuum_area;
  
  vector<PeakTruthInfo> gamma_lines;
};//struct ExpectedPhotopeakInfo


struct DataSrcInfo
{
  string detector_name;
  string location_name;
  string live_time_name;
  
  InjectSourceInfo src_info;
  
  vector<ExpectedPhotopeakInfo> expected_photopeaks;
};//struct DataSrcInfo

ExpectedPhotopeakInfo create_expected_photopeak( const InjectSourceInfo &info, const vector<PeakTruthInfo> &lines )
{
  assert( !lines.empty() );
  
  ExpectedPhotopeakInfo peak;
  peak.gamma_lines = lines;
  
  peak.peak_area = 0.0;
  
  // We'll weight the mean and ROI lower/upper by each lines relative area
  //  - not really sure if this is the right thing to do...
  peak.roi_lower = 0.0;
  peak.roi_upper = 0.0;
  peak.effective_energy = 0.0;
  
  for( const auto &l : lines )
  {
    peak.peak_area += l.area;
    peak.roi_lower += l.area * (l.energy - 0.66*l.full_width); // Using 0.66 instead of 0.5 because skew, or whatever, I assume, so none of the dataset ever has more predicted than actually in the data
    peak.roi_upper += l.area * (l.energy + 0.66*l.full_width);
    peak.effective_energy += l.area * l.energy;
  }//for( const auto &l : lines )
  
  peak.roi_lower /= peak.peak_area;
  peak.roi_upper /= peak.peak_area;
  peak.effective_energy /= peak.peak_area;
  
  // Use the non-Poisson varied spectrum to calculate total area
  const double total_area = info.src_no_poisson->gamma_integral( static_cast<float>(peak.roi_lower),
                                                              static_cast<float>(peak.roi_upper) );
    
  assert( ((total_area - peak.peak_area) > -0.000*peak.peak_area)
         || ((total_area - peak.peak_area) > -5.0)
         || (peak.roi_upper > info.src_no_poisson->gamma_energy_max()) );
  peak.continuum_area = std::max( 0.0, total_area - peak.peak_area );
  peak.nsigma_over_background = peak.peak_area / sqrt( std::max(1.0, peak.continuum_area) );
  
  
  //if( peak.effective_energy > 1407.5 && peak.effective_energy < 1408.5 )
  //{
  //  cout << "This one" << endl;
  //}
  
  assert( peak.roi_lower < peak.roi_upper );
  assert( peak.effective_energy < peak.roi_upper );
  assert( peak.effective_energy > peak.roi_lower );
  
  return peak;
}//ExpectedPhotopeakInfo create_expected_photopeak(...)




InjectSourceInfo parse_inject_source_files( const string &base_name )
{
  /*
  CSV file is ordered as:
   DRF  Detector/CZT/1.5cm-2cm-2cm
   Source  Br76_Sh
   Distance  50
   MeasurementTime  1800
   BackgroundLocation  Livermore
   InjectDesc  Source + Background
              
             
   Source Lines:
   Energy (keV)  CPS  Area  FWHM  FullWidth  Label
   100.3  8.39348e-09  1.51083e-05  4.51382  12.7386  D.E.
   178.7  1.89767e-07  0.00034158  5.49774  15.7744  EscapeXRay
   182.23  4.80212e-07  0.000864382  5.53458  15.8892  EscapeXRay
   209.7  0.000884833  1.5927  5.80637  16.7385  Peak
   
   ...
   
   Background Lines:
   Energy (keV)  CPS  Area  FWHM  FullWidth  Label
   0.46  1.65619e-08  2.98114e-05  2.98932  7.44028  EscapeXRay
   0.6  8.78397e-11  1.58111e-07  2.98932  7.47825  EscapeXRay
   3  9.11796e-11  1.64123e-07  2.98932  7.73021  EscapeXRay
   
   */
  
  const string pcf_filename = base_name + ".pcf";
  const string csv_filename = base_name + "_gamma_lines.csv";
  
  auto spec_file = make_shared<SpecUtils::SpecFile>();
  const bool loaded_pcf = spec_file->load_pcf_file( pcf_filename );
  assert( loaded_pcf );
  if( !loaded_pcf )
    throw runtime_error( "failed to load " + pcf_filename );
  
  const vector<shared_ptr<const SpecUtils::Measurement>> meass = spec_file->measurements();
  assert( meass.size() == 16 );
  if( meass.size() != 16 )
    throw runtime_error( "Unexpected number of measurements in " + pcf_filename );
  
  InjectSourceInfo info;
  info.file_base_path = base_name;
  info.src_name = SpecUtils::filename( base_name );
  info.spec_file = spec_file;
  
  if( meass[0]->title() != info.src_name )
    cerr << "Mismatch of src: '" << meass[0]->title() << "', vs '" << info.src_name << "'" << endl;
  
  assert( meass[0]->title() == info.src_name );
  assert( meass[1]->title() == "Background" );
  assert( meass[2]->title() == "Background" );
  assert( meass[3]->title() == info.src_name );
  assert( meass[4]->title() == "Background (non-Poisson)" );
  
  info.src_spectra.push_back( meass[0] );
  info.short_background = meass[1];
  info.long_background = meass[2];
  info.src_no_poisson = meass[3];
  info.background_no_poisson = meass[4];
  
  for( size_t i = 5; i < meass.size(); ++i )
  {
    assert( meass[i]->title() == info.src_name );
    info.src_spectra.push_back( meass[i] );
  }
  
  
  ifstream csv_strm( csv_filename.c_str(), (ios::binary | ios::in) );
  assert( csv_strm.good() );
  if( !csv_strm.good() )
    throw runtime_error( "Failed to open CSV: " + csv_filename );
  
  bool found_src_lines = false, found_background_lines = false;
  string line;
  while( !found_src_lines && SpecUtils::safe_get_line(csv_strm, line) )
    found_src_lines = SpecUtils::istarts_with(line, "Source Lines:");
  
  assert( found_src_lines );
  if( !found_src_lines )
    throw runtime_error( "Failed to find Source Lines in " + csv_filename );
  
  SpecUtils::safe_get_line(csv_strm, line);
  assert( line == "Energy (keV),CPS,Area,FWHM,FullWidth,Label" );
  if( line != "Energy (keV),CPS,Area,FWHM,FullWidth,Label" )
    throw runtime_error( "Unexpected foreground CSV header: " + line );
  
  while( SpecUtils::safe_get_line(csv_strm, line) )
  {
    found_background_lines = SpecUtils::istarts_with(line, "Background Lines:");
    if( found_background_lines )
      break;
    
    SpecUtils::trim( line );
    if( line.empty() )
      continue;
    
    info.source_lines.emplace_back( line );
  }//while( loop over foreground lines )
  
  assert( found_background_lines );
  SpecUtils::safe_get_line(csv_strm, line);
  assert( line == "Energy (keV),CPS,Area,FWHM,FullWidth,Label" );
  if( line != "Energy (keV),CPS,Area,FWHM,FullWidth,Label" )
    throw runtime_error( "Unexpected background CSV header: " + line );
  
  
  while( SpecUtils::safe_get_line(csv_strm, line) )
  {
    SpecUtils::trim( line );
    if( line.empty() )
      continue;
    
    info.background_lines.emplace_back( line );
  }//while( loop over background lines )
  
  return info;
}//InjectSourceInfo parse_inject_source_files( const string &base_name )



int main( int argc, char **argv )
{
  const double start_wall = SpecUtils::get_wall_time();
  const double start_cpu = SpecUtils::get_cpu_time();
  
  string datadir;
  if( datadir.empty() )
  {
    string targetfile = "data/CharacteristicGammas.txt";
    if( !AppUtils::locate_file( targetfile, false, 5, false ) )
    {
      cerr << "Unable to find '" << targetfile << "' directory." << endl;
      return -24;
    }
    
    datadir = SpecUtils::parent_path( targetfile );
  }//if( docroot.empty() )
  
  assert( SpecUtils::is_directory(datadir) );
  
  InterSpec::setStaticDataDirectory( datadir );
  
  const string base_dir = "/Users/wcjohns/rad_ana/peak_area_optimization/peak_fit_accuracy_inject/";
  
  vector<DetectorInjectSet> inject_sets;
  
  for( boost::filesystem::directory_iterator detector_itr(base_dir);
      detector_itr != boost::filesystem::directory_iterator(); ++detector_itr )
  {
    //detector_itr->path().filename()
    
    if( !boost::filesystem::is_directory(detector_itr->status()) )
      continue;
    
    const boost::filesystem::path detector_path = detector_itr->path();
    
    if( detector_path.filename() != "Detective-EX" )
    {
      cerr << "Skipping detector " << detector_path.filename() << endl;
      continue;
    }
    
    cout << "In detector directory: " << detector_path.filename() << endl;
    
    for( boost::filesystem::directory_iterator city_itr(detector_path);
        city_itr != boost::filesystem::directory_iterator(); ++city_itr )
    {
      if( !boost::filesystem::is_directory(city_itr->status()) )
        continue;
      
      const boost::filesystem::path city_path = city_itr->path();
      cout << "In city directory: " << city_path.filename() << endl;
      
      if( city_path.filename() != "Livermore" )
      {
        cout << "Skipping city " << city_itr->path().filename() << endl;
        continue;
      }
      
      // Now loop over {"30_seconds", "300_seconds", "1800_seconds"}
      for( boost::filesystem::directory_iterator livetime_itr(city_path);
          livetime_itr != boost::filesystem::directory_iterator(); ++livetime_itr )
      {
        if( !boost::filesystem::is_directory(livetime_itr->status()) )
          continue;
        
        const boost::filesystem::path livetime_path = livetime_itr->path();
        
        if( livetime_path.filename() != "300_seconds" )
        {
          cout << "Skipping live-time " << livetime_path.filename() << endl;
          continue;
        }
        
        //cout << "In livetime: " << livetime_path << endl;
        
        // Now loop over PCF files
        vector<string> files_to_load_basenames;
        for( boost::filesystem::directory_iterator file_itr(livetime_path);
            file_itr != boost::filesystem::directory_iterator(); ++file_itr )
        {
          if( !boost::filesystem::is_regular_file(file_itr->status()) )
            continue;
          
          const string pcf_name = file_itr->path().string();
          if( !SpecUtils::iends_with(pcf_name, ".pcf") )
            continue;
          
          const string base_name = pcf_name.substr(0, pcf_name.size() - 4 );
          
          //base_name + "_raw_gamma_lines.txt"
          const string csv_name = base_name + "_gamma_lines.csv";
          
          if( !boost::filesystem::is_regular_file(csv_name) )
          {
            cerr << "No gamma lines CSV for " << file_itr->path() << " - skipping" << endl;
            continue;
          }
          
          //debug_printout = true;
          //if( base_name.find( "Eu152_Sh") == string::npos )
          //  continue;
          
          files_to_load_basenames.push_back( base_name );
        }//for( loop over files in directory )
        
        const size_t num_srcs = files_to_load_basenames.size();
        
        inject_sets.push_back( DetectorInjectSet{} );
        
        DetectorInjectSet &injects = inject_sets.back();
        injects.detector_name = detector_path.filename().string();
        injects.location_name = city_path.filename().string();
        injects.live_time_name = livetime_path.filename().string();
        injects.source_infos.resize( num_srcs );
        
        SpecUtilsAsync::ThreadPool pool;
        
        for( size_t file_index = 0; file_index < num_srcs; ++file_index )
        {
          InjectSourceInfo *info = &(injects.source_infos[file_index]);
          const string base_name = files_to_load_basenames[file_index];
          
          pool.post( [info,base_name](){
            *info = parse_inject_source_files( base_name );
          } );
        }//for( loop over files to parse )
        
        pool.join();
        
        // A smoke check to make sure nothing got messed up
        assert( injects.source_infos.size() == num_srcs );
        for( size_t index = 0; index < num_srcs; ++index )
        {
          assert( !injects.source_infos[index].source_lines.empty() );
          assert( !injects.source_infos[index].background_lines.empty() );
          assert( injects.source_infos[index].file_base_path == files_to_load_basenames[index] );
          assert( !injects.source_infos[index].src_spectra.empty() );
          assert( injects.source_infos[index].src_spectra[0]->title() == SpecUtils::filename(files_to_load_basenames[index]) );
        }//for( size_t index = 0; index < injects.source_infos.size(); ++index )
        
        cout << "Parsed " << injects.source_infos.size() << " sources" << endl;
      }//for( loop over live-time directories, livetime_itr )
    }//for( loop over cities, city_itr )
  }//for( loop over detector types )
  
  
  
  // Lets not worry about super-small peaks, even where there is little to no background
  const double min_peak_area = 5;
  
  //Photopeak clusters below this next number of sigma will be discarded, since we really
  //  shouldnt find these peaks
  const double not_expected_thresh_nsigma = 1.0;
  
  vector<DataSrcInfo> input_srcs;
  size_t num_inputs = 0, num_accepted_inputs = 0;
  for( const DetectorInjectSet &inject_set : inject_sets )
  {
    for( const InjectSourceInfo &info : inject_set.source_infos )
    {
      // For the moment, we'll look for visible lines on the first Poisson varied source
      //cout << "For " << inject_set.detector_name << "/"
      //<< inject_set.live_time_name << "/" << inject_set.location_name
      //<< " source " << info.src_name << ": ";
      //cout << "For '" << info.file_base_path << "':" << endl;
      
      num_inputs += 1;
      
      // The "truth" lines are all before random-summing, so the observed peaks will be smaller
      //  in area than the truth CSV says, but if we keep dead-time low, this wont be noticeable.
      //  Right now have chosen 2%, fairly arbitrarily
      //  For EX-100
      //   1%: Lose 14 out of 223 files
      //   2%: Lose  8 out of 223 files
      //   5%: Lose  4 out of 223 files
      //  10%: Lose  4 out of 223 files
      if( info.src_no_poisson->live_time() < 0.98*info.src_no_poisson->real_time() )
      {
        continue;
      }
      
      const shared_ptr<const SpecUtils::Measurement> meas = info.src_spectra[0];
      vector<PeakTruthInfo> lines = info.source_lines;
      lines.insert( end(lines), begin(info.background_lines), end(info.background_lines) );
      
      //for( const PeakTruthInfo &line : lines )
      //  cout << "Source line: " << line.energy << " keV, area=" << line.area << endl;
      
      /**
       - Sort lines from largest area to smallest.
       - cluster, using `~0.75*FWHM`
       - See what is hit, and whats not
       */
      
      std::sort( begin(lines), end(lines), []( const PeakTruthInfo &lhs, const PeakTruthInfo &rhs ){
        return (lhs.area < rhs.area);
      } );
      
      
      const double cluster_fwhm_multiple = 0.5;
      vector<vector<PeakTruthInfo>> clustered_lines;
      while( !lines.empty() )
      {
        const PeakTruthInfo main_line = std::move( lines.back() );
        lines.resize( lines.size() - 1 ); //remove the line we just grabbed
        
        vector<PeakTruthInfo> cluster;
        cluster.push_back( main_line );
        
        // Look all through `lines` for lines within `cluster_fwhm_multiple` of main_line
        deque<size_t> index_to_remove;
        for( size_t i = 0; i < lines.size(); ++i )
        {
          if( fabs(lines[i].energy - main_line.energy) < cluster_fwhm_multiple*main_line.fwhm )
          {
            index_to_remove.push_back( i );
            cluster.push_back( lines[i] );
          }
        }//for( loop over lines to cluster )
        
        for( auto iter = std::rbegin(index_to_remove); iter != std::rend(index_to_remove); ++iter )
        {
//          cout << "Removing " << *iter << ", which has energy " << lines[*iter].energy
//          << ", main line=" << main_line.energy << endl;
          lines.erase( begin(lines) + (*iter) );
        }
        
        clustered_lines.push_back( cluster );
      }//while( !lines.empty() )
      
      
      
      vector<ExpectedPhotopeakInfo> detectable_clusters;
      for( const vector<PeakTruthInfo> &cluster : clustered_lines )
      {
        //cout << "cluster: {";
        //for( const auto &c : cluster )
        //  cout << "{" << c.energy << "," << c.area << "}, ";
        //cout << "}" << endl;
        
        assert( !cluster.empty() );
        if( cluster.empty() )
          continue;
        
        const PeakTruthInfo &main_line = cluster.front();
        if( (main_line.area < 5.0)  // Let be realistic, and require something
           || (fabs(main_line.energy - 478.0) < 1.2)  // Avoid Alpha-Li reaction
           || (fabs(main_line.energy - 511.0) < 1.0 ) // Avoid D.B. 511
              // Shouldn't there be some other broadened reaction lines here???
           || ((main_line.energy + 0.5*main_line.full_width) > info.src_no_poisson->gamma_energy_max()) // Make sure not off upper-end
           || ((main_line.energy - 0.5*main_line.full_width) < info.src_no_poisson->gamma_energy_min()) // Make sure not below lower-end
           || ((main_line.energy - 0.5*main_line.full_width) < 50.0) //eh, kinda arbitrary, maybe shouldnt limit?
           )
        {
          continue;
        }
        
        ExpectedPhotopeakInfo roi_info = create_expected_photopeak( info, cluster );
        
        if( (roi_info.peak_area > min_peak_area)
           && (roi_info.nsigma_over_background > not_expected_thresh_nsigma) )
          detectable_clusters.push_back( std::move(roi_info) );
      }//for( const vector<PeakTruthInfo> &cluster : clustered_lines )
      
      if( !detectable_clusters.empty() )
      {
        num_accepted_inputs += 1;
        DataSrcInfo src_info;
        src_info.detector_name = inject_set.detector_name;
        src_info.location_name = inject_set.location_name;
        src_info.live_time_name = inject_set.live_time_name;
        
        src_info.src_info = info;
        src_info.expected_photopeaks = detectable_clusters;
        
        input_srcs.push_back( src_info );
      }//if( !detectable_clusters.empty() )
    }//for( const InjectSourceInfo &src : source_infos )
  }//for( const DetectorInjectSet &inject_set : inject_sets )
  
  cout << "Used " << num_accepted_inputs << " of total " << num_inputs << " input files." << endl;
  
  const double def_want_nsigma = 8;   // i.e., above 4 sigma, lets weight all peaks the same
  const double lower_want_nsigma = 2; // The number of sigma above which we will positively reward finding a peak
  // Between `def_want_nsigma` and `lower_want_nsigma` we will linearly weight for not finding a peak
  //not_expected_thresh_nsigma - the threshold at which we will start punishing if a peak is not
  //                             expected at this threshold; these source photopeaks have already
  //                             been removed.
  const double found_extra_punishment = 0.25; // 1/this-value gives the trade-off of finding extra peaks, verses not finding peaks
  
  auto eval_settings = [&input_srcs, def_want_nsigma, lower_want_nsigma, found_extra_punishment]( 
                                                            const FindCandidateSettings settings,
                                                            const bool write_n42 )
      -> tuple<double,size_t,size_t,size_t> //<score, num_peaks_found, num_peaks_not_found, num_extra_peaks>
  {
    double sum_score = 0.0;
    size_t num_peaks_not_found = 0, num_extra_peaks = 0, num_peaks_found = 0;
    
    for( const DataSrcInfo &info : input_srcs )
    {
      vector<tuple<float,float,float>> detected_peaks; //{mean, sigma, amplitude}
      shared_ptr<const SpecUtils::Measurement> src_spectrum = info.src_info.src_spectra.front();
      find_candidate_peaks( src_spectrum, 0, 0, detected_peaks, settings );
      
      const vector<tuple<float,float,float>> orig_peak_candidates = detected_peaks;
      
      double score = 0.0;
      size_t num_detected_expected = 0, num_detected_not_expected = 0, num_not_detected = 0;
      
      // First, go through expected peaks, and match up to candidates
      vector<ExpectedPhotopeakInfo> expected_photopeaks = info.expected_photopeaks;
      set<size_t> det_photopeak_remove;
      
      //for( const auto &p : info.expected_photopeaks )
      //  cout << "Expected " << p.effective_energy << " keV nsgima=" << p.nsigma_over_background << endl;
      
      vector<ExpectedPhotopeakInfo> not_detected;
      vector<tuple<float,float,float>> detected_expected, detected_not_expected;
      not_detected.reserve( expected_photopeaks.size() );
      detected_expected.reserve( expected_photopeaks.size() );
      
      
      for( size_t expected_index = 0; expected_index < expected_photopeaks.size(); ++expected_index )
      {
        const ExpectedPhotopeakInfo &expected = expected_photopeaks[expected_index];
        
        vector<pair<tuple<float,float,float>,size_t>> matching_det; //{mean, sigma, amplitude}
        for( size_t det_index = 0; det_index < detected_peaks.size(); ++det_index )
        {
          const tuple<float,float,float> &det_peak = detected_peaks[det_index];
          const float mean = get<0>(det_peak);
          //const float sigma = get<1>(det_peak);
          //const float amp = get<1>(det_peak);
          
          if( (mean > expected.roi_lower) && (mean < expected.roi_upper) )
          {
            matching_det.push_back( make_pair(det_peak,det_index) );
            detected_expected.push_back( det_peak );
          }
        }//for( size_t det_index = 0; det_index < detected_peaks.size(); ++i )
        
        if( matching_det.empty() )
        {
          num_not_detected += 1;
          
          not_detected.push_back( expected );
          /*
           if( expected.nsigma_over_background < lower_want_nsigma )
           score -= 0; //No punishment
           else if( expected.nsigma_over_background < def_want_nsigma )
           score -= ( (def_want_nsigma - expected.nsigma_over_background) / (expected.nsigma_over_background - lower_want_nsigma) );
           else
           score -= 1;
           */
        }else
        {
          num_detected_expected += 1;
          
          // we got a match
          if( expected.nsigma_over_background > def_want_nsigma )
            score += 1.0;
          else if( expected.nsigma_over_background > lower_want_nsigma )
            score += ((def_want_nsigma - expected.nsigma_over_background) / (expected.nsigma_over_background - lower_want_nsigma));
          else
            score += 0.0; //No punishment, but no reward.
        }//if( matching_det.empty() ) / else
        
        for( const auto &i : matching_det )
          det_photopeak_remove.insert( i.second );
      }//for( loop over expected_photopeaks )
      
      for( auto iter = std::rbegin(det_photopeak_remove); iter != std::rend(det_photopeak_remove); ++iter )
      {
        const size_t index = *iter;
        detected_peaks.erase( detected_peaks.begin() + index );
      }
      
      detected_not_expected.swap( detected_peaks );
      num_detected_not_expected = detected_not_expected.size();
      
      score -= found_extra_punishment * num_detected_not_expected;
      
      /*
      cout << "For " << info.src_info.file_base_path
      //<< info.src_info.src_name
      //<< "/" << info.detector_name
      //<< "/" << info.location_name
      //<< "/" << info.live_time_name
      << ":\n"
      << "\tnum_detected_expected=" << num_detected_expected
      << ", num_detected_not_expected=" << num_detected_not_expected
      << ", num_not_detected=" << num_not_detected
      << ", score=" << score << endl;
      */
      
      sum_score += score;
      num_peaks_found += num_detected_expected;
      num_peaks_not_found += num_not_detected;
      num_extra_peaks += num_detected_not_expected;
      
      if( write_n42 )
      {
        //const string src_dir = SpecUtils::parent_path( info.src_info.file_base_path );
        
        string outdir = "/Users/wcjohns/rad_ana/InterSpec/target/peak_fit_improve/build_xcode/output_n42";
        if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;
        
        outdir = SpecUtils::append_path( outdir, info.detector_name );
        if( !SpecUtils::is_directory(outdir) && SpecUtils::create_directory(outdir) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;
        
        outdir = SpecUtils::append_path( outdir, info.location_name );
        if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;
        
        outdir = SpecUtils::append_path( outdir, info.live_time_name );
        if( !SpecUtils::is_directory(outdir) && !SpecUtils::create_directory(outdir) )
          cerr << "Failed to create directory '" << outdir << "'" << endl;
        
        const string out_n42 = SpecUtils::append_path( outdir, info.src_info.src_name ) + ".n42";
        
        SpecMeas output;
        
        output.add_remark( "Failed to find " + std::to_string(num_not_detected) + " peaks that are expected." );
        output.add_remark( "Found " + std::to_string(num_detected_expected) + " peaks that were expected." );
        output.add_remark( "Found " + std::to_string(num_detected_not_expected) + " peaks that were NOT expected." );
        
        shared_ptr<SpecUtils::Measurement> out_with_cand = make_shared<SpecUtils::Measurement>( *src_spectrum );
        out_with_cand->set_sample_number( 1 );
        
        const string title = info.detector_name + "/" + info.src_info.src_name;
        out_with_cand->set_title( title + " + candidate peaks" );
        
        output.add_measurement( out_with_cand, false );
        
        deque<shared_ptr<const PeakDef>> peaks;
        
        vector<tuple<float,float,float>> all_candidates = detected_expected;
        all_candidates.insert( end(all_candidates), begin(detected_not_expected), end(detected_not_expected) );
        for( size_t peak_index = 0; peak_index < all_candidates.size(); ++peak_index )
        {
          const tuple<float,float,float> &p = all_candidates[peak_index]; //{mean, sigma, amplitude}
          const float mean = get<0>(p);
          const float sigma = get<1>(p);
          const float amp = get<2>(p);
          
          const bool is_expected = (peak_index < detected_expected.size());
          
          auto peak = make_shared<PeakDef>( mean, sigma, amp );
          peak->setFitFor( PeakDef::CoefficientType::Mean, false );
          peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
          peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
          peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
          peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
          peak->continuum()->calc_linear_continuum_eqn( src_spectrum, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
          peak->setLineColor( is_expected ? Wt::GlobalColor::darkGreen : Wt::GlobalColor::red );
          
          if( is_expected )
          {
            //Put truth-level info to user-label field
            string info_str;
            for( const ExpectedPhotopeakInfo &exp_info : info.expected_photopeaks )
            {
              if( mean > exp_info.roi_lower && mean < exp_info.roi_upper )
              {
                if( !info_str.empty() )
                  info_str += ".\n";
                
                info_str += "E: A=" + SpecUtils::printCompact(exp_info.peak_area, 3)
                + ", W=" + SpecUtils::printCompact(exp_info.gamma_lines.front().fwhm, 3)
                + ", #s=" + SpecUtils::printCompact(exp_info.nsigma_over_background, 3);
              }
            }
            
            peak->setUserLabel( info_str );
          }else
          {
            peak->setUserLabel( "Not Expected" );
          }
          
          peaks.push_back( peak );
        }//for( const tuple<float,float,float> &p : all_candidates )
        
        
        for( const ExpectedPhotopeakInfo &missing : not_detected )
        {
          const double mean = missing.effective_energy;
          const double sigma = missing.gamma_lines.front().fwhm/2.35482f;
          auto peak = make_shared<PeakDef>( mean, sigma, missing.peak_area );
          peak->setFitFor( PeakDef::CoefficientType::Mean, false );
          peak->setFitFor( PeakDef::CoefficientType::Sigma, false );
          peak->setFitFor( PeakDef::CoefficientType::GaussAmplitude, false );
          peak->continuum()->setType( PeakContinuum::OffsetType::Linear );
          peak->continuum()->setRange( mean - 3*sigma, mean + 3*sigma );
          peak->continuum()->calc_linear_continuum_eqn( src_spectrum, mean, mean - 4*sigma, mean + 4*sigma, 5, 5 );
          peak->setLineColor( Wt::GlobalColor::darkGray );
          peak->setUserLabel( "N_sigma=" + SpecUtils::printCompact(missing.nsigma_over_background, 3) );
          
          peaks.push_back( peak );
        }//for( const ExpectedPhotopeakInfo &missing : not_detected )
        
        
        output.setPeaks( peaks, {1} );
        
        {// begin add second derivative to N42
          auto second_deriv = make_shared<vector<float>>();
          const int side_bins = settings.num_smooth_side_channels;
          const int poly_order = settings.smooth_polynomial_order;
          smoothSpectrum( src_spectrum, side_bins, poly_order, 2, *second_deriv );
          
          shared_ptr<SpecUtils::Measurement> second_deriv_meas = make_shared<SpecUtils::Measurement>();
          second_deriv_meas->set_gamma_counts( second_deriv, 1.0, 1.0 );
          second_deriv_meas->set_energy_calibration( src_spectrum->energy_calibration() );
          second_deriv_meas->set_sample_number( 2 );
          second_deriv_meas->set_title( title + " + smoothed second derivative." );
          output.add_measurement( second_deriv_meas, true );
          
          deque<shared_ptr<const PeakDef>> candidates;
          for( const tuple<float,float,float> &p : orig_peak_candidates )
          {
            auto peak = make_shared<PeakDef>( get<0>(p), get<1>(p), get<2>(p) );
            peak->continuum()->setType( PeakContinuum::OffsetType::Constant );
            peak->continuum()->setPolynomialCoef(0, 0.0);
            candidates.push_back( peak );
          }
          
          output.setPeaks( candidates, {2} );
        }// end add second derivative to N42
        
        
        {
          // I think this will be good for estimating sigma and area
          const int side_bins = 3;
          const int poly_order = 3;
          vector<float> rougher_second_deriv;
          smoothSpectrum( src_spectrum, side_bins, poly_order, 2, rougher_second_deriv );
          shared_ptr<SpecUtils::Measurement> rough_second_deriv_meas = make_shared<SpecUtils::Measurement>();
          auto rough_second_deriv_counts = make_shared<vector<float>>( rougher_second_deriv );
          rough_second_deriv_meas->set_gamma_counts( rough_second_deriv_counts, 1.0, 1.0 );
          rough_second_deriv_meas->set_energy_calibration( src_spectrum->energy_calibration() );
          rough_second_deriv_meas->set_sample_number( 3 );
          rough_second_deriv_meas->set_title( title + " + rougher smoothed second derivative." );
          output.add_measurement( rough_second_deriv_meas, true );
        }
        
        //ofstream outstrm( out_n42.c_str(), ios::out | ios::binary );
        //if( !outstrm )
        //  cerr << "Failed to open '" << out_n42 << "'" << endl;
        //output.write( outstrm, output.sample_numbers(), output.detector_names(), SpecUtils::SaveSpectrumAsType::N42_2012 );
        
        output.save2012N42File( out_n42, [=](){
          cerr << "Failed to write '" << out_n42 << "'" << endl;
        });
        
        //cout << "Wrote '" << out_n42 << endl;
      }//if( write_n42 )
    }//for( const DataSrcInfo &info : input_srcs )
    
    //cout << "Avrg score: " << sum_score/input_srcs.size() << endl;
    
    return {sum_score / input_srcs.size(), num_peaks_found, num_peaks_not_found, num_extra_peaks};
  };//eval_settings lambda
  
  
  /**
   Default settings: Avrg score: 24.4401
   */
  FindCandidateSettings best_settings;
  best_settings.num_smooth_side_channels = 4; // low res more
  best_settings.smooth_polynomial_order = 3;  // highres 3, lowres 2
  best_settings.threshold_FOM = 1.3;
  best_settings.pos_sum_threshold_sf = -0.01f;
  
  double best_sum_score = -1.0E-6;
  size_t best_score_num_peaks_found = 0, best_score_num_peaks_not_found = 0, best_score_num_extra_peaks = 0;
  
  std::mutex score_mutex;
  auto eval_settings_fcn = [&]( const FindCandidateSettings settings ){
    
    const tuple<double,size_t,size_t,size_t> result = eval_settings( settings, false );
  
    const double score = std::get<0>(result);
    const size_t num_peaks_found = std::get<1>(result);
    const size_t num_peaks_not_found = std::get<2>(result);
    const size_t num_extra_peaks = std::get<3>(result);
    
    std::lock_guard<std::mutex> lock( score_mutex );
    if( score > best_sum_score )
    {
      best_sum_score = score;
      best_settings = settings;
      best_score_num_peaks_found = num_peaks_found;
      best_score_num_peaks_not_found = num_peaks_not_found;
      best_score_num_extra_peaks = num_extra_peaks;
    }
  };
  
  
  best_settings.num_smooth_side_channels = 13; // low res more
  best_settings.smooth_polynomial_order = 2;  // highres 3, lowres 2
  best_settings.threshold_FOM = 1.3;
  best_settings.pos_sum_threshold_sf = 0.16f;
  eval_settings_fcn( best_settings );
  
  /*
  SpecUtilsAsync::ThreadPool pool;
  size_t num_posted = 0;
  for( int num_side = 2; num_side < 14; ++num_side )
  {
    for( int poly_order = 2; poly_order < 4; ++poly_order )
    {
      for( double threshold_FOM = 0.6; threshold_FOM < 6; threshold_FOM += 0.1 )
      {
        for( float pos_sum_threshold_sf = -0.2f; pos_sum_threshold_sf < 0.01; pos_sum_threshold_sf += 0.01 )
        {
          ++num_posted;
          FindCandidateSettings settings;
          settings.num_smooth_side_channels = num_side; // low res more
          settings.smooth_polynomial_order = poly_order;  // highres 3, lowres 2
          settings.threshold_FOM = threshold_FOM;
          settings.pos_sum_threshold_sf = pos_sum_threshold_sf;
          
          pool.post( [eval_settings_fcn,settings,&score_mutex](){
            try
            {
              eval_settings_fcn( settings );
            }catch( std::exception &e )
            {
              //std::lock_guard<std::mutex> lock( score_mutex );
              //cerr << "Caught exception: " << e.what() << ", for:" << endl
              //<< "\tsettings.num_smooth_side_channels = " << settings.num_smooth_side_channels << ";" << endl
              //<< "\tsettings.smooth_polynomial_order = " << settings.smooth_polynomial_order << ";" << endl
              //<< "\tsettings.threshold_FOM = " << settings.threshold_FOM << ";" << endl
              //<< "\tsettings.pos_sum_threshold_sf = " << settings.pos_sum_threshold_sf << ";" << endl
              //<< endl;
            }
          } );
        }//for( float pos_sum_threshold_sf = -0.1f; pos_sum_threshold_sf < 0.1; pos_sum_threshold_sf += 0.01 )
      }//for( double threshold_FOM = 0.8; threshold_FOM < 2.5; threshold_FOM += 0.1 )
    }//for( int poly_order = 1; poly_order < 4; ++poly_order )
  }//for( int num_side = 1; num_side < 8; ++num_side )
  
  cout << "Posted " << num_posted << " evaluations." << endl;
  
  pool.join();
  */

  cout << "found_extra_punishment = " << found_extra_punishment << endl;
  cout << "Best settings had score " << best_sum_score << ", with values:" << endl
  << "\tsettings.num_smooth_side_channels = " << best_settings.num_smooth_side_channels << ";" << endl
  << "\tsettings.smooth_polynomial_order = " << best_settings.smooth_polynomial_order << ";" << endl
  << "\tsettings.threshold_FOM = " << best_settings.threshold_FOM << ";" << endl
  << "\tsettings.pos_sum_threshold_sf = " << best_settings.pos_sum_threshold_sf << ";" << endl
  << "And:\n"
  << "\tnum_peaks_found:" << best_score_num_peaks_found << endl
  << "\tnum_peaks_not_found:" << best_score_num_peaks_not_found << endl
  << "\tnum_extra_peaks:" << best_score_num_extra_peaks << endl
  << endl << endl;
   
  
  //cerr << "Setting best_settings instead of actually finding them!" << endl;
  //best_settings.num_smooth_side_channels = 12;
  //best_settings.smooth_polynomial_order = 2;
  //best_settings.threshold_FOM = 1.8;
  //best_settings.pos_sum_threshold_sf = -0.04;
  
  eval_settings( best_settings, true );
  cout << "Wrote N42 with best settings." << endl;
  
  const double end_wall = SpecUtils::get_wall_time();
  const double end_cpu = SpecUtils::get_cpu_time();
  
  cout << "Ran in {wall=" << (end_wall - start_wall)
        << ", cpu=" << (end_cpu - start_cpu) << "} seconds" << endl;
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


