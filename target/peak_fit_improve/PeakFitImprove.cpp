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
#include <numeric>
#include <iostream>

#include <boost/program_options.hpp>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MakeDrfFit.h"
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


float nai_fwhm_fcn( const float energy )
{
  static const vector<float> nai_fwhm_coefs{ -4.0f, 6.3f, 0.6f };   //"NaI 3x3"
  return DetectorPeakResponse::peakResolutionFWHM( energy,
                DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, nai_fwhm_coefs );
}//float nai_fwhm_fcn( const float energy )


float labr_fwhm_fcn( const float energy )
{
  static const vector<float> labr_fwhm_coefs{ 5.0f, 3.0f, 0.55f };  //"LaBr 10%"
  return DetectorPeakResponse::peakResolutionFWHM( energy,
                DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, labr_fwhm_coefs );
}//float labr_fwhm_fcn( const float energy )

float hpge_fwhm_fcn( const float energy )
{
  static const vector<float> hpge_fwhm_coefs{ 1.55f, 0.25f, 0.35f };//"HPGe 40%"
  return DetectorPeakResponse::peakResolutionFWHM( energy,
    DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, hpge_fwhm_coefs );
}//float hpge_fwhm_fcn( const float energy )


CoarseResolutionType coarse_resolution_from_peaks( const vector<shared_ptr<const PeakDef>> &peaks )
{
  size_t num_peaks = 0;
  double max_sig = 0.0;
  double low_w = 0, med_w = 0, high_w = 0, all_w = 0;
  
  for( const auto &p : peaks )
  {
    if( !p || (p->mean() > 3000) || (p->mean() < 50) )
      continue;
    
    num_peaks += 1;
    
    const double drf_low_fwhm = nai_fwhm_fcn( p->mean() );
    const double drf_med_fwhm = labr_fwhm_fcn( p->mean() );
    const double drf_high_fwhm = hpge_fwhm_fcn( p->mean() );
    
    const double chi_dof = p->chi2dof();
    const double stat_sig = p->peakArea() / p->peakAreaUncert();
    max_sig = std::max( max_sig, stat_sig );
    
    const double w = std::min( stat_sig, 10.0 );// / std::max( chi_dof, 0.5 );
    
    all_w += w;
    
    const double low_diff = fabs(p->fwhm() - drf_low_fwhm);
    const double med_diff = fabs(p->fwhm() - drf_med_fwhm);
    const double high_diff = fabs(p->fwhm() - drf_high_fwhm);
    if( (low_diff < med_diff) && (low_diff < high_diff) )
      low_w += w;
    else if( med_diff < high_diff )
      med_w += w;
    else
      high_w += w;
  }//for( const auto &p : peak_candidates )
  
  if( (num_peaks == 1) && (max_sig < 5) )
    return CoarseResolutionType::Unknown;
  
  if( all_w <= 0.0 )
    return CoarseResolutionType::Unknown;
  
  //if( low_w > high_w )
  //  return CoarseResolutionType::Low;
  
  
  if( (low_w > med_w) && (low_w > high_w) )
    return CoarseResolutionType::Low;
  
  if( (med_w > high_w) )
    return CoarseResolutionType::Medium;
  
  return CoarseResolutionType::High;
}//CoarseResolutionType coarse_resolution_from_peaks( const vector<shared_ptr<const PeakDef>> &peaks )


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
  //const vector<string> num_channels{ "2048_channels" };
  const vector<string> det_types{ "CZT", "CsI", "HPGe", "LaBr3", "NaI" };
  //const vector<string> det_types{ "NaI" };
  //const vector<string> time_strs{ "30_seconds", "300_seconds", "1800_seconds" };
  const vector<string> time_strs{ "300_seconds" };

  
  //Also see:
  // void expected_peak_width_limits( const float energy, const bool highres, float &min_sigma_width_kev, float &max_sigma_width_kev )
  
  SpecUtilsAsync::ThreadPool pool;
  
  std::mutex intemed_res_mutex;
  std::map<CoarseResolutionType,vector<double>> intemed_res;
  std::map<CoarseResolutionType,tuple<int,int,int,int>> peak_to_drf_comp;
  
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
     
      if( expected_type == CoarseResolutionType::High )
      {
        if( (nchan_str != "4096_channels") && (nchan_str != "8192_channels") )
          continue;
      }
      
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
          
          //cout << "For " << det_name << " (" << det_type_str << " - " << time_str << ") there are "
          //     << pcf_files.size() << " files." << endl;
          
          const auto analyze_file = [&intemed_res_mutex, &intemed_res, &peak_to_drf_comp]( const string filename, const CoarseResolutionType expected  ){
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
            
            const size_t nchannel = foreground->num_gamma_channels();
            assert( nchannel > 64 );
            
            //vector<tuple<float,float,float>> peak_candidates; //{mean,sigma,area}
            //secondDerivativePeakCanidates( foreground, 0, nchannel - 1, peak_candidates );
            
            //const double x0 = foreground->gamma_channel_lower(0);
            //const double x1 = foreground->gamma_channel_lower(nchannel-1);
            //const double ncausalitysigma = 0.0;
            //const double stat_threshold = 0.0;
            //const double hypothesis_threshold = 0.0;
            //const bool isRefit = false;
            //const vector<PeakDef> peak_candidates = fitPeaksInRange( x0, x1,
            //                              ncausalitysigma, stat_threshold, hypothesis_threshold,
            //                              {}, foreground, {}, isRefit );
            
            vector<shared_ptr<const PeakDef>> peak_candidates
              = ExperimentalAutomatedPeakSearch::search_for_peaks( foreground, nullptr, {}, true );
            
            //cout << "'" << filename << "'" << endl;
            //cout << "peak_candidates.size()=" << peak_candidates.size() << endl;
            
            double cs137_fwhm = 0.0;
            if( peak_candidates.size() > 0 )
            {
              try
              {
                auto peaks = make_shared< deque<shared_ptr<const PeakDef>> >();
                //for( const auto &p : peak_candidates )
                //  peaks->push_back( make_shared<PeakDef>(get<0>(p), get<1>(p), get<2>(p)) );
                
                //for( const auto &p : peak_candidates )
                //  peaks->push_back( make_shared<PeakDef>(p) );
                
                size_t num_closer_NaI = 0, num_closer_HPGe = 0;
                for( const auto &p : peak_candidates )
                {
                  peaks->push_back( p );
                  const double nai_fwhm = nai_fwhm_fcn( p->mean() );
                  const double hpge_fwhm = hpge_fwhm_fcn( p->mean() );
                  
                  if( fabs(nai_fwhm - p->fwhm()) < fabs(hpge_fwhm - p->fwhm()) )
                    num_closer_NaI += 1;
                  else
                    num_closer_HPGe += 1;
                }
                
                DetectorPeakResponse::ResolutionFnctForm form = DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial;
                if( peak_candidates.size() == 1 )
                  form = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
                
                const bool highResolution = (num_closer_HPGe > num_closer_NaI); //This only matters if LLS fit to FWHM fails, or uses GADRAS FWHM
                
                // Lets go through 
                
                const int sqrtEqnOrder = 2; //
                vector<float> fwhm_coeffs, fwhm_coeff_uncerts;
                MakeDrfFit::performResolutionFit( peaks, form, highResolution, sqrtEqnOrder,
                                                 fwhm_coeffs, fwhm_coeff_uncerts );
                
                double initial_cs137_fwhm = DetectorPeakResponse::peakResolutionFWHM( 661, form, fwhm_coeffs );
                //cout << "Initial FWHM=" << 100*initial_cs137_fwhm/661 << "%" << endl;
                if( IsNan(initial_cs137_fwhm) && (peak_candidates.size() > 1) )
                {
                  // This can happen if we have 2 or 3 peaks, below 661 keV, so that by the time we
                  //  get up to 661, the equation is invalid
                  peaks->erase( begin(*peaks) );
                  
                  //if( peak_candidates.size() == 1 )
                    form = DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn;
                  
                  MakeDrfFit::performResolutionFit( peaks, form, highResolution, sqrtEqnOrder,
                                                   fwhm_coeffs, fwhm_coeff_uncerts );
                  
                  initial_cs137_fwhm = DetectorPeakResponse::peakResolutionFWHM( 661, form, fwhm_coeffs );
                  //cout << "After removing lowest energy peak FWHM=" << 100*initial_cs137_fwhm/661 << "%" << endl;
                }//if( IsNan(initial_cs137_fwhm) && (peak_candidates.size() > 1) )
                
                if( IsNan(initial_cs137_fwhm) || ((100*initial_cs137_fwhm/661) < 1) )
                {
                  cout << "Peaks: ";
                  for( const auto &p : peak_candidates )
                    cout << "{m=" << p->mean() << ", fwhm=" << p->fwhm() << "}, ";
                  cout << endl;
                  cout << "Coefs: {";
                  for( const auto l : fwhm_coeffs )
                    cout << l << ", ";
                  cout << "}" << endl;
                  
                  if( IsNan(initial_cs137_fwhm) )
                  {
                    cout << "Not sure..." << endl;
                  }
                }
                
                // Go through and remove outliers....
                if( peaks->size() > 5 )
                {
                  auto new_peaks = make_shared< deque<shared_ptr<const PeakDef>> >();
                  vector<pair<double,shared_ptr<const PeakDef>>> distances;
                  for( const auto &p : *peaks )
                  {
                    double pred_fwhm = 0.0;
                    const double mean_MeV = 0.001 * p->mean();  //kSqrtPolynomial
                    for( size_t i = 0; i < fwhm_coeffs.size(); ++i )
                      pred_fwhm += fwhm_coeffs[i] * std::pow( mean_MeV, 1.0*i );
                    pred_fwhm = sqrt( pred_fwhm );
                    //pred_fwhm = fwhm_coeffs[0] + fwhm_coeffs[1]*sqrt( p->mean() ); //
                    
                    const double frac_diff = fabs( p->fwhm() - pred_fwhm ) / p->fwhm();
                    if( !IsNan(frac_diff) && !IsInf(frac_diff) )
                      distances.emplace_back( frac_diff, p );
                  }//for( const auto &p : initial_fit_peaks )
                  
                  std::sort( begin(distances), end(distances),
                    []( const pair<double,shared_ptr<const PeakDef>> &lhs,
                        const pair<double,shared_ptr<const PeakDef>> &rhs) -> bool {
                      return lhs.first > rhs.first;
                    } );

                  // Limit to un-selecting max of 20% of peaks (arbitrarily chosen), if the deviate
                  // more than 25% from the fit (again, arbitrarily chosen).
                  const size_t max_remove = static_cast<size_t>( std::ceil( 0.2*distances.size() ) );
                  for( size_t index = 0; index < distances.size(); ++index )
                  {
                    if( (distances[index].first < 0.25) || (index > max_remove) )
                      new_peaks->push_back( distances[index].second );
                  }
                  
                  MakeDrfFit::performResolutionFit( new_peaks, form, highResolution, sqrtEqnOrder,
                                                   fwhm_coeffs, fwhm_coeff_uncerts );
                }//if( peaks->size() > 5 )
                
                cs137_fwhm = DetectorPeakResponse::peakResolutionFWHM( 661, form, fwhm_coeffs );
                
                //cout << "Final FWHM=" << 100*cs137_fwhm/661 << "%" << endl;
              }catch( std::exception &e )
              {
                cerr << ("Caught exception: " + string(e.what()) + "\n");
              }
            }//if( peak_candidates.size() > 0 )
            
            
            const CoarseResolutionType type_from_peaks = coarse_resolution_from_peaks( peak_candidates );
            
            if( !peak_candidates.empty() && (type_from_peaks == CoarseResolutionType::Low && expected == CoarseResolutionType::High) )
            {
              cerr << ("File '" + filename + "' failed.\n");
              cerr << endl;
            }
            
            std::lock_guard<std::mutex> lock( intemed_res_mutex );
            intemed_res[expected].push_back( 100 * cs137_fwhm / 661.0 );
            
            switch( type_from_peaks )
            {
              case CoarseResolutionType::Low:
                get<1>(peak_to_drf_comp[expected]) += 1;
                break;
              case CoarseResolutionType::Medium:
                get<2>(peak_to_drf_comp[expected]) += 1;
                break;
              case CoarseResolutionType::High:
                get<3>(peak_to_drf_comp[expected]) += 1;
                break;
              case CoarseResolutionType::Unknown:
                get<0>(peak_to_drf_comp[expected]) += 1;
                break;
            }//switch( type_from_peaks )
          };//analyze_file lambda
          
          for( const string input : pcf_files )
          {
            pool.post( [=](){ analyze_file( input, expected_type ); } );
            //analyze_file( input, expected_type );
          }//for( const string input : pcf_files )
        }//for( const string &time_str : time_strs )
      }//for( const string &det_path : detector_paths )
    }//for( const string &det_type_str : det_types )
  }//for( const string &nchan_str : num_channels )
  
  pool.join();
  
  double min_corr = std::numeric_limits<double>::max(), max_corr = 0.0;
  
  for( const auto &t_v : intemed_res )
  {
    for( const auto &v : t_v.second )
    {
      min_corr = std::min( min_corr, v );
      max_corr = std::max( max_corr, v );
    }
  }
  
  max_corr = std::min( max_corr, 25.0 );
  
  const size_t num_channel = 1024;
  std::map<CoarseResolutionType,vector<float>> channel_correlation_hists;
  for( const auto &t_v : intemed_res )
    channel_correlation_hists[t_v.first] = vector<float>( num_channel, 0.0f );
  
  auto cal = make_shared<SpecUtils::EnergyCalibration>();
  cal->set_full_range_fraction( num_channel, {0.0f, static_cast<float>(max_corr)}, {} );
  
  for( const auto &t_v : intemed_res )
  {
    for( const auto &v : t_v.second )
    {
      const double channel = cal->channel_for_energy(v);
      const int index = std::min( std::max(0,static_cast<int>(floor(channel))), static_cast<int>(num_channel) - 1 );
      channel_correlation_hists[t_v.first][index] += 1;
    }
  }
  
  SpecUtils::SpecFile output_hists;
  for( const auto &t_v : intemed_res )
  {
    auto meas = make_shared<SpecUtils::Measurement>();
    auto counts = make_shared<vector<float>>( channel_correlation_hists[t_v.first] );
    const double ncounts = std::accumulate( begin(*counts), end(*counts), 0.0f );
    meas->set_gamma_counts( counts, ncounts, ncounts );
    meas->set_energy_calibration( cal );
    meas->set_title( title_str(t_v.first) );
    output_hists.add_measurement( meas, false );
  }
  
  
  for( const auto &t_v : peak_to_drf_comp )
  {
    auto meas = make_shared<SpecUtils::Measurement>();
    auto counts = make_shared<vector<float>>( 16, 0.0 );
    (*counts)[1] = get<0>( t_v.second );
    (*counts)[3] = get<1>( t_v.second );
    (*counts)[5] = get<2>( t_v.second );
    (*counts)[7] = get<3>( t_v.second );
    
    const double ncounts = std::accumulate( begin(*counts), end(*counts), 0.0f );
    meas->set_gamma_counts( counts, ncounts, ncounts );
    
    auto cal = make_shared<SpecUtils::EnergyCalibration>();
    cal->set_full_range_fraction( 16, {0.0f, 16.0f}, {} );
    
    meas->set_energy_calibration( cal );
    meas->set_title( string(title_str(t_v.first)) + " Nearest Type" );
    output_hists.add_measurement( meas, false );
  }
  const SpecUtils::time_point_t time = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
  output_hists.set_uuid( SpecUtils::to_iso_string(time) );
  output_hists.cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );
  
  {// Begin write output
    const char *outname = "results.n42";
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


