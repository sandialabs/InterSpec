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
#include <mutex>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecServer.h"

using namespace std;


enum class CoarseResolutionType : int
{
  /** Csi, NaI*/
  Low,
  
  /** LaBr, CZT */
  Medium,
  
  /** HPGe */
  High,
  
  /** Unknown */
  Unknown
};//enum class CoarseResolutionType

const char *title_str( const CoarseResolutionType type )
{
  switch( type )
  {
    case CoarseResolutionType::Low:     return "Low Resolution";
    case CoarseResolutionType::Medium:  return "Medium Resolution";
    case CoarseResolutionType::High:    return "High Resolution";
    case CoarseResolutionType::Unknown: return "Unknown Resolution";
  }//switch( type )
  
  assert( 0 );
  return "";
}//const char *title_str( const CoarseResolutionType type )


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
  
  const string base_dir = "/Users/wcjohns/rad_ana/peak_area_optimization/resolution_study_inject/";
  
  const vector<string> num_channels{ "512_channels", "1024_channels", "2048_channels", "4096_channels", "8192_channels" };
  const vector<string> det_types{ "CZT", "CsI", "HPGe", "LaBr3", "NaI" };
  // Detector names
  const vector<string> time_strs{ "30_seconds", "300_seconds", "1800_seconds" };
  
  SpecUtilsAsync::ThreadPool pool;
  
  std::mutex channel_correlations_mutex;
  std::map<CoarseResolutionType,vector<double>> channel_correlations;
  
  for( const string &nchan_str : num_channels )
  {
    const string nchan_path = SpecUtils::append_path(base_dir, nchan_str);
    for( const string &det_type_str : det_types )
    {
      CoarseResolutionType expected_type = CoarseResolutionType::Unknown;
      if( (det_type_str == "CsI") || (det_type_str == "NaI") )
        expected_type = CoarseResolutionType::Low;
      else if( (det_type_str == "CZT") || (det_type_str == "LaBr3") )
        expected_type = CoarseResolutionType::Medium;
      else if( det_type_str == "HPGe" )
        expected_type = CoarseResolutionType::High;
      else
        throw runtime_error( "Unknown detector type: " + det_type_str );
      
      assert( expected_type != CoarseResolutionType::Unknown );
      
      const string det_type_path = SpecUtils::append_path(nchan_path, det_type_str);
      const vector<string> detector_paths =  SpecUtils::ls_directories_in_directory( det_type_path );
      
      for( const string &det_path : detector_paths )
      {
        const string det_name = SpecUtils::filename(det_path);
        
        for( const string &time_str : time_strs )
        {
          const string final_parent_path = SpecUtils::append_path(det_path, time_str);
          const vector<string> pcf_files = SpecUtils::recursive_ls( final_parent_path, ".pcf" );
          
          cout << "For " << det_name << " (" << det_type_str << " - " << time_str << ") there are " << pcf_files.size() << " files." << endl;
          
          const auto analyze_file = [&channel_correlations_mutex, &channel_correlations]( const string filename, const CoarseResolutionType expected  ){
            SpecUtils::SpecFile file;
            if( !file.load_file(filename, SpecUtils::ParserType::Pcf ) )
            {
              cerr << ("Failed to open '" + filename + "'") << endl;
              return;
            }
            
            assert( file.num_measurements() == 2 );
            if( file.num_measurements() != 2 )
            {
              cerr << ("File '" + filename + "' had " + std::to_string(file.num_measurements()) 
                       + " measurements (expected 2)") << endl;
              return;
            }
            
            shared_ptr<const SpecUtils::Measurement> foreground = file.measurements()[0];
            shared_ptr<const SpecUtils::Measurement> background = file.measurements()[1];
            
            assert( foreground && background && (background->title() == "Background") );
            
            /*
            secondDerivativePeakCanidates( const std::shared_ptr<const SpecUtils::Measurement> data,
                                          size_t start_channel,
                                          size_t end_channel,
                                          std::vector< std::tuple<float,float,float> > &results );
             */
            
            assert( !!foreground->gamma_channel_contents() );
            const vector<float> &channel_contents = *foreground->gamma_channel_contents();
            double corr = 0.0;
            for( size_t i = 1; i < channel_contents.size(); ++i )
            {
              const float &prev_val = channel_contents[i-1];
              const float &this_val = channel_contents[i];
              const float diff = fabs(this_val - prev_val);
              const float avrg = 0.5f*(prev_val + this_val);
              const double this_corr = diff / sqrt( ((avrg > 1.0) ? avrg : 1.0f) );
              corr += this_corr;
            }
            
            corr /= channel_contents.size();
            
            
            std::lock_guard<std::mutex> lock( channel_correlations_mutex );
            channel_correlations[expected].push_back( corr );
          };//analyze_file lambda
          
          
          for( const string input : pcf_files )
          {
            pool.post( [=](){ analyze_file( input, expected_type ); } );
          }//for( const string input : pcf_files )
        }//for( const string &time_str : time_strs )
      }//for( const string &det_path : detector_paths )
    }//for( const string &det_type_str : det_types )
  }//for( const string &nchan_str : num_channels )
  
  pool.join();
  
  double min_corr = std::numeric_limits<double>::max(), max_corr = 0.0;
  
  for( const auto &t_v : channel_correlations )
  {
    for( const auto &v : t_v.second )
    {
      min_corr = std::min( min_corr, v );
      max_corr = std::max( max_corr, v );
    }
  }
  
  const size_t num_channel = 512;
  std::map<CoarseResolutionType,vector<float>> channel_correlation_hists;
  for( const auto &t_v : channel_correlations )
    channel_correlation_hists[t_v.first] = vector<float>( num_channel, 0.0f );
  
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_full_range_fraction( num_channel, {0.0f, static_cast<float>(max_corr)}, {} );
  
  for( const auto &t_v : channel_correlations )
  {
    for( const auto &v : t_v.second )
    {
      const double channel = cal->channel_for_energy(v);
      const int index = std::min( std::max(0,static_cast<int>(floor(channel))), static_cast<int>(num_channel) - 1 );
      channel_correlation_hists[t_v.first][index] += 1;
    }
  }
  
  SpecUtils::SpecFile output_hists;
  for( const auto &t_v : channel_correlations )
  {
    auto meas = make_shared<SpecUtils::Measurement>();
    auto counts = make_shared<vector<float>>( channel_correlation_hists[t_v.first] );
    meas->set_gamma_counts( counts, 1.0f, 1.0f );
    meas->set_energy_calibration( cal );
    meas->set_title( title_str(t_v.first) );
    output_hists.add_measurement( meas, false );
  }
  output_hists.cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );
  
  {// Begin write output
    const char *outname = "correlation.n42";
    ofstream coor_output( outname );
    if( output_hists.write_2012_N42( coor_output ) )
      cout << "Wrote '" << outname << "'" << endl;
    else
      cout << "Failed to write '" << outname << "'" << endl;
  }// End write output
  
  /*
  void smoothSpectrum( const std::vector<float> &spectrum, const int side_bins,
                      const int order, const int derivative,
                      std::vector<float> &results );
  void secondDerivativePeakCanidates( const std::shared_ptr<const SpecUtils::Measurement> data,
                                      size_t start_channel,
                                      size_t end_channel,
                                     std::vector< std::tuple<float,float,float> > &results );
  */
  
  const double end_wall = SpecUtils::get_wall_time();
  const double end_cpu = SpecUtils::get_cpu_time();
  
  cout << "Ran in {wall=" << (end_wall - start_wall)
        << ", cpu=" << (end_cpu - start_cpu) << "} seconds" << endl;
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


