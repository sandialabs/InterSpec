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
#include <atomic>
#include <chrono>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>

//#include "/external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"
#include <Wt/Json/Value>
#include <Wt/Json/Array>
#include <Wt/Json/Parser>
#include <Wt/Json/Object>
#include <Wt/Json/Serializer>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/D3SpectrumExport.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/MakeDrfFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/ReferenceLineInfo.h"

#include "InitialFit_GA.h"
#include "PeakFitImproveData.h"

using namespace std;

namespace PeakFitImproveData
{

namespace G2k
{
  vector<G2kPeak> g2k_peaks_from_file( istream &strm )
  {
    bool found_peaks_start = false;
    string line;
    while( !found_peaks_start && SpecUtils::safe_get_line(strm, line) )
    {
      SpecUtils::trim( line );
      if( !found_peaks_start
         && SpecUtils::icontains( line, "*                              Peak Search Report                               *") )
      {
        while( !found_peaks_start && SpecUtils::safe_get_line(strm, line) )
        {
          SpecUtils::trim( line );
          if( SpecUtils::istarts_with( line, "-----" ) )
          {
            found_peaks_start = true;
            break;
          }
        }
      }
    }//while( !found_peaks_start && SpecUtils::safe_get_line(strm, line) )

    if( !found_peaks_start )
      throw runtime_error( "Failed to find peaks part of report" );

    vector<G2kPeak> result;
    while( SpecUtils::safe_get_line(strm, line) )
    {
      SpecUtils::trim( line );
      if( line.empty() )
        continue;

      if( SpecUtils::istarts_with(line, "Errors" ) )
        break;

      vector<string> fields;
      SpecUtils::split( fields, line, " \t");

      assert( (fields.size() == 8) || (fields.size() == 9) );
      if( (fields.size() != 8) && (fields.size() != 9) )
        throw runtime_error( "Got line w/o 8 or 9 fields: '" + line + "'" );

      assert( (fields.size() == 8) || (fields[0] == "M") || (fields[0] == "m") );
      if( (fields.size() == 9) && ((fields[0] != "M") && (fields[0] != "m")) )
        throw runtime_error( "First field is not M or m: '" + line + "'" );

      const char mutltiplet = (fields.size() == 9) ? fields[0].at(0) : ' ';
      if( fields.size() == 9 )
      {
        line = line.substr( 1 );
        SpecUtils::trim( line );
      }

      vector<float> values;
      SpecUtils::split_to_floats( line.c_str(), values, " \t", false );
      assert( values.size() == 8 );

      G2kPeak peak;
      peak.multiplet = mutltiplet;
      peak.PeakNum = static_cast<int>( values[0] );
      peak.Energy = values[1];
      peak.PeakSig = values[2];
      peak.FWHM = values[3];
      peak.ContinuumCounts = values[4];
      peak.NetPeakArea = values[5];
      peak.NetAreaError = 0.5*values[6]; //the RPT files give uncert at 2-sigma
      peak.NetCountRate = values[7];

      assert( result.empty() || ((result.back().PeakNum + 1) == peak.PeakNum) );
      result.push_back( peak );
    }//while( SpecUtils::safe_get_line(strm, line) )

    return result;
  };//auto g2k_peaks_from_file
}//namespace G2k


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
  peak.effective_fwhm = 0.0;

  double weighted_variance_numer = 0.0, weighted_variance_denom = 0.0;
  double largest_line_area = 0.0, fwhm_of_largest_line = 0.0;
  for( const auto &l : lines )
  {
    peak.peak_area += l.area;
    peak.roi_lower += l.area * (l.energy - 0.66*l.full_width); // Using 0.66 instead of 0.5 because skew, or whatever, I assume, so none of the dataset ever has more predicted than actually in the data
    peak.roi_upper += l.area * (l.energy + 0.66*l.full_width);
    peak.effective_energy += l.area * l.energy;

    if( l.area > largest_line_area )
    {
      largest_line_area = l.area;
      fwhm_of_largest_line = l.fwhm;
    }

    weighted_variance_numer += l.area * pow(l.fwhm/2.35482, 3.0);
    weighted_variance_denom += l.area * l.fwhm/2.35482;
  }//for( const auto &l : lines )


  peak.roi_lower /= peak.peak_area;
  peak.roi_upper /= peak.peak_area;
  peak.effective_energy /= peak.peak_area;

  // We have to calc `variance_from_seperation_num` after we have already calculated peak.effective_energy
  double variance_from_seperation_num = 0.0;
  for( const auto &l : lines )
    variance_from_seperation_num += l.area * (l.full_width/2.35482) * std::pow(l.energy - peak.effective_energy, 2.0);

  peak.effective_fwhm = 2.35482 * sqrt( (weighted_variance_numer + variance_from_seperation_num)
                             / weighted_variance_denom );

  //cout << "For truth-peak " << peak.effective_energy << " keV, Effective_FWHM="
  //     << peak.effective_fwhm << " vs largest line for peak had FWHM=" << fwhm_of_largest_line << endl;

  // Use the non-Poisson varied spectrum to calculate total area
  const double total_area = info.src_no_poisson->gamma_integral( static_cast<float>(peak.roi_lower),
                                                              static_cast<float>(peak.roi_upper) );

  if( !( ((total_area - peak.peak_area) > -8.0*sqrt(peak.peak_area))
         || ((total_area - peak.peak_area) > -5.0)
         || (peak.roi_upper > info.src_no_poisson->gamma_energy_max()) ) )
  {
    // This looks to trigger for very large are peak (like 500k counts), where there is a good
    //   amount of peak-skew.
    cerr << "Warning: create_expected_photopeak: For peak at " << peak.effective_energy << " keV,"
    << " got peak.peak_area=" << peak.peak_area << ", while histogram total_area=" << total_area
    << " between " << peak.roi_lower << " and " << peak.roi_upper << " keV, for \n\t"
    << info.file_base_path << "." << endl;
  }

  assert( (total_area > 1.0E5) //large skew with such areas
         || ((total_area - peak.peak_area) > -8.0*sqrt(peak.peak_area))
         || ((total_area - peak.peak_area) > -5.0)
         || (peak.roi_upper > info.src_no_poisson->gamma_energy_max())
         || (info.src_no_poisson->live_time() < 0.85*info.src_no_poisson->real_time()) );

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


InjectSourceInfo parse_inject_source_files( const string &base_name, const std::vector<std::pair<float,float>> &deviation_pairs )
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
  //assert( loaded_pcf );
  if( !loaded_pcf )
    throw runtime_error( "failed to load " + pcf_filename );

  const vector<shared_ptr<const SpecUtils::Measurement>> meass = spec_file->measurements();
  assert( meass.size() == 16 );
  if( meass.size() != 16 )
    throw runtime_error( "Unexpected number of measurements in " + pcf_filename );


  // The PCF files dont contain the nonlinear deviation pairs that GADRAS created the spectrum files with, so we will
  //  add them in here - if we have them.
  if( !deviation_pairs.empty() )
  {
    map<shared_ptr<const SpecUtils::EnergyCalibration>,shared_ptr<const SpecUtils::EnergyCalibration>> updated_cals;
    for( const shared_ptr<const SpecUtils::Measurement> &m : meass )
    {
      shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = m ? m->energy_calibration() : nullptr;
      assert( orig_cal && orig_cal->valid() );
      if( !orig_cal || !orig_cal->valid() )
        continue;

      const auto prev_pos = updated_cals.find(orig_cal);
      if( prev_pos != end(updated_cals) )
      {
        spec_file->set_energy_calibration( prev_pos->second, m );
        continue;
      }

      assert( orig_cal->deviation_pairs().empty() );

      assert( orig_cal->type() == SpecUtils::EnergyCalType::FullRangeFraction );
      if( orig_cal->type() != SpecUtils::EnergyCalType::FullRangeFraction )
        throw runtime_error( "unexpected energy cal type!" );

      auto new_cal = make_shared<SpecUtils::EnergyCalibration>();
      new_cal->set_full_range_fraction( orig_cal->num_channels(), orig_cal->coefficients(), deviation_pairs );

      updated_cals[orig_cal] = new_cal;
      spec_file->set_energy_calibration( new_cal, m );
    }//for( const shared_ptr<const SpecUtils::Measurement> &m : meass )
  }//if( !deviation_pairs.empty() )


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


std::tuple<std::vector<DetectorInjectSet>,std::vector<DataSrcInfo>> load_inject_data_with_truth_info(
                                                                        const std::string &base_dir,
                                                                        const std::vector<std::string> &wanted_detectors,
                                                                        const std::vector<std::string> &live_times,
                                                                        const std::vector<std::string> &wanted_cities )
{
  vector<DetectorInjectSet> inject_sets;

  std::mutex input_srcs_mutex;
  vector<DataSrcInfo> input_srcs;
  vector<bool> used_detector( wanted_detectors.size(), false ), used_live_times( live_times.size(), false );

  for( boost::filesystem::directory_iterator detector_itr(base_dir);
      detector_itr != boost::filesystem::directory_iterator(); ++detector_itr )
  {
    //detector_itr->path().filename()

    if( !boost::filesystem::is_directory(detector_itr->status()) )
      continue;

    const boost::filesystem::path detector_path = detector_itr->path();

    bool is_wanted_det = false;
    for( size_t i = 0; i < wanted_detectors.size(); ++i )
    {
      const string &det = wanted_detectors[i];
      const bool want = (detector_path.filename() == det);
      is_wanted_det |= want;
      used_detector[i] = (used_detector[i] || want);
    }

    if( !is_wanted_det )
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

      const auto wanted_pos = std::find(begin(wanted_cities), end(wanted_cities), city_path.filename() );
      if( !wanted_cities.empty() && (wanted_pos == end(wanted_cities)) )
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


        //if( PeakFitImprove::debug_printout )
        {

          bool is_wanted_lt = false;
          for( size_t i = 0; i < live_times.size(); ++i )
          {
            const string &lt = live_times[i];
            const bool wanted = (livetime_path.filename() == lt);
            is_wanted_lt |= wanted;
            used_live_times[i] = (used_live_times[i] || wanted);
          }

          if( !is_wanted_lt )
          {
            cerr << "Skipping live-time " << livetime_path.filename() << endl;
            continue;
          }
        }//if( PeakFitImprove::debug_printout )

        //cout << "In livetime: " << livetime_path << endl;
        const boost::filesystem::path deviation_path = livetime_path / "Deviation.gadras";

        vector<pair<float,float>> deviation_pairs;
        if( boost::filesystem::is_regular_file(deviation_path) )
        {
          const string dev_path_str = deviation_path.string();
          std::ifstream file( dev_path_str.c_str(), ios::in | ios::binary );
          assert( file.is_open() );
          if( !file.is_open() )
            std::runtime_error( "Error: Could not open file " + dev_path_str );
          string line;
          while( SpecUtils::safe_get_line(file, line, 1024) )
          {
            SpecUtils::trim(line);
            if( line.empty() || !isnumber(line[0]) )
              continue;
            std::istringstream iss(line);
            float first, second;
            if (iss >> first >> second)
              deviation_pairs.emplace_back(first, second);
            else
              throw runtime_error( "Failed to parse Deviation.gadras" );
          }
        }//if( boost::filesystem::is_regular_file(deviation_path) )


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
            assert( 0 );
            cerr << "No gamma lines CSV for " << file_itr->path() << " - skipping" << endl;
            continue;
          }

          if( PeakFitImprove::debug_printout )
          {
            const vector<string> wanted_sources{
              //"Ac225_Unsh",
              "Eu152_Sh",
              //"Fe59_Phant",
              //"Am241_Unsh",
              //"Ho166m_Unsh"
            };
            bool wanted = wanted_sources.empty();
            for( const string &src : wanted_sources )
            {
              const bool want = (base_name.find(src) != string::npos);
              wanted |= want;
              if( want )
                cout << "Setting info for '" << src << "'" << endl;
            }
            if( !wanted )
              continue;
          }//if( PeakFitImprove::debug_printout )

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

          pool.post( [info,base_name,deviation_pairs](){
            try
            {
              *info = PeakFitImproveData::parse_inject_source_files( base_name, deviation_pairs );
            }catch( std::exception &e )
            {
              //Baltimore/1800_seconds/Re188_Sh.pcf doesnt load
              cerr << "Failed to load file: " << e.what() << endl;
            }
          } );
        }//for( loop over files to parse )

        pool.join();

        //injects.source_infos.erase( std::remove_if( begin(injects.source_infos), end(injects.source_infos),
        // []( const InjectSourceInfo &info ){
        //  return !info.spec_file;
        //}), end(injects.source_infos) );

        // A smoke check to make sure nothing got messed up - let 5 files not parse though
        assert( injects.source_infos.size() == num_srcs );
        size_t num_sucess = 0;
        for( size_t index = 0; index < num_srcs; ++index )
        {
          if( !injects.source_infos[index].spec_file )
            continue;

          assert( !injects.source_infos[index].source_lines.empty() );
          assert( !injects.source_infos[index].background_lines.empty() );
          assert( injects.source_infos[index].file_base_path == files_to_load_basenames[index] );
          assert( !injects.source_infos[index].src_spectra.empty() );
          assert( injects.source_infos[index].src_spectra[0]->title() == SpecUtils::filename(files_to_load_basenames[index]) );

          num_sucess += 1;
        }//for( size_t index = 0; index < injects.source_infos.size(); ++index )
        assert( (num_sucess + 5) >= num_srcs );

        cout << "Parsed " << num_sucess << " sources" << endl;
      }//for( loop over live-time directories, livetime_itr )
    }//for( loop over cities, city_itr )
  }//for( loop over detector types )

  std::atomic<size_t> num_inputs = 0, num_accepted_inputs = 0;

  size_t num_queued = 0;
  SpecUtilsAsync::ThreadPool pool;
  for( const DetectorInjectSet &inject_set : inject_sets )
  {
    for( const InjectSourceInfo &info : inject_set.source_infos )
    {
      const auto create_DataSrcInfo = [&num_inputs, &num_accepted_inputs, &inject_set, &input_srcs_mutex, &input_srcs]( const InjectSourceInfo &info ){
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

        if( !info.spec_file || !info.src_no_poisson || info.src_spectra.empty() || !info.short_background )
        {
          cerr << "load_inject_data_with_truth_info: Skipping '" << info.src_name << "' - as it looks like it didnt read in" << endl;
          return;
        }

        //#if( !RETURN_PeakDef_Candidates )
#if( !WRITE_ALL_SPEC_TO_HTML )
#warning "Only using <2% dead-time files"
        if( info.src_no_poisson->live_time() < 0.98*info.src_no_poisson->real_time() )
        {
          return;
        }
#endif //#if( !WRITE_ALL_SPEC_TO_HTML )
        //#else
        //#warning "Using all Live-Times for peak search"
        //#endif

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

          assert( !cluster.empty() );
          if( cluster.empty() )
            continue;

          const PeakTruthInfo &main_line = cluster.front();
          if( (main_line.area < 4.0)  // Let be realistic, and require something
             || (fabs(main_line.energy - 478.0) < 1.2)  // Avoid Alpha-Li reaction
             || (fabs(main_line.energy - 511.0) < 1.0 ) // Avoid D.B. 511
             // Shouldn't there be some other broadened reaction lines here???
             || ((main_line.energy + 0.5*main_line.full_width) > info.src_no_poisson->gamma_energy_max()) // Make sure not off upper-end
             || ((main_line.energy - 0.5*main_line.full_width) < info.src_no_poisson->gamma_energy_min()) // Make sure not below lower-end
             || ((main_line.energy - 0.5*main_line.full_width) < 5.0) //eh, kinda arbitrary, maybe shouldnt limit?
             )
          {
            //cout << "Skipping cluster at " << main_line.energy << endl;
            continue;
          }

          ExpectedPhotopeakInfo roi_info = PeakFitImproveData::create_expected_photopeak( info, cluster );

          if( (roi_info.peak_area > JudgmentFactors::min_truth_peak_area)
             && (roi_info.nsigma_over_background > JudgmentFactors::min_truth_nsigma) )
            detectable_clusters.push_back( std::move(roi_info) );
          //else
          //  cout << "Discarding undetectable cluster at " << roi_info.effective_energy << endl;
        }//for( const vector<PeakTruthInfo> &cluster : clustered_lines )

        std::sort( begin(detectable_clusters), end(detectable_clusters),
                  []( const ExpectedPhotopeakInfo &lhs, const ExpectedPhotopeakInfo &rhs ) -> bool {
          return lhs.effective_energy < rhs.effective_energy;
        } );

        if( PeakFitImprove::debug_printout )
        {
          /*
           for( const ExpectedPhotopeakInfo &roi_info : detectable_clusters )
           {
           cout << "Expected ROI: {energy: " << roi_info.effective_energy
           << ", fwhm: " << roi_info.gamma_lines.front().fwhm
           << ", PeakArea: " << roi_info.peak_area
           << ", ContinuumArea: " << roi_info.continuum_area
           << ", NSigma: " << roi_info.nsigma_over_background
           << ", ROI: [" << roi_info.roi_lower << ", " << roi_info.roi_upper << "]"
           << "}"
           << endl;

           //cout << "cluster: {";
           //for( const auto &c : cluster )
           //  cout << "{" << c.energy << "," << c.area << "}, ";
           //cout << "}" << endl;
           }//for( const ExpectedPhotopeakInfo &roi_info : detectable_clusters )
           */
        }//if( PeakFitImprove::debug_printout )


        if( !detectable_clusters.empty() )
        {
          num_accepted_inputs += 1;
          DataSrcInfo src_info;
          src_info.detector_name = inject_set.detector_name;
          src_info.location_name = inject_set.location_name;
          src_info.live_time_name = inject_set.live_time_name;

          src_info.src_info = info;
          src_info.expected_photopeaks = detectable_clusters;

          std::lock_guard<std::mutex> lock( input_srcs_mutex );

          input_srcs.push_back( src_info );
        }//if( !detectable_clusters.empty() )
      };//create_DataSrcInfo lambda

      num_queued += 1;
      if( num_queued > 100 )
      {
        pool.join();
        num_queued = 0;
      }

      pool.post( [&](){
        create_DataSrcInfo( info );
      } );
    }//for( const InjectSourceInfo &info : source_infos )
  }//for( const DetectorInjectSet &inject_set : inject_sets )

  pool.join();


  for( size_t i = 0; i < wanted_detectors.size(); ++i )
  {
    if( !used_detector[i] )
    {
      cerr << "Failed to use detector '" << wanted_detectors[i] << "'!" << endl;
      assert( 0 );
    }
  }

  for( size_t i = 0; i < used_live_times.size(); ++i )
  {
    if( !used_live_times[i] )
    {
      cerr << "Failed to use live-time '" << used_live_times[i] << "'!" << endl;
      assert( 0 );
    }
  }

  cout << "Used " << num_accepted_inputs.load() << " of total " << num_inputs.load() << " input files." << endl;

  std::lock_guard<std::mutex> lock( input_srcs_mutex ); //dont really need, I dont think

  return tuple<vector<DetectorInjectSet>,vector<DataSrcInfo>>{ std::move(inject_sets), std::move(input_srcs) };
}//std::tuple<std::vector<DetectorInjectSet>,std::vector<DataSrcInfo>> load_inject_data_with_truth_info( )


#if( WRITE_ALL_SPEC_TO_HTML )
void write_html_summary( const vector<DataSrcInfo> &src_infos )
{
  ofstream output( "mikes_inject.html" );

  D3SpectrumExport::write_html_page_header( output, "GADRAS Simulations of likely field nuclides" );

  output << "<body>" << endl;

  //Hack putting <style></style> block in HTML body, but wahtever, seems ot work
  output <<"<style>"
  << ".TopLinesTable{ margin-top: 25px; margin-left: auto; margin-right: auto; border-collapse: collapse; border: 1px solid black; }" << endl
  << "table, th, td{ border: 1px solid black; }" << endl
  << "fieldset{width: 90vw; margin-left: auto; margin-right: auto; margin-top: 20px;}" << endl
  << "</style>" << endl;


  for( size_t i = 0; i < src_infos.size(); ++i )
  {
    const DataSrcInfo &info = src_infos[i];
    const InjectSourceInfo &src_info = info.src_info;

    const string &detector_name = info.detector_name;
    const string &location_name = info.location_name;
    const string &live_time_name = info.live_time_name;


    const vector<PeakTruthInfo> &source_lines = src_info.source_lines;
    const vector<PeakTruthInfo> &background_lines = src_info.background_lines;
    //src_info.spec_file;
    const vector<shared_ptr<const SpecUtils::Measurement>> &src_spectra = src_info.src_spectra;
    assert( !src_spectra.empty() );
    const shared_ptr<const SpecUtils::Measurement> &spectrum = src_spectra.front(); //We'll just plot the first spectrum only
    const shared_ptr<const SpecUtils::Measurement> &short_background = src_info.short_background;
    const shared_ptr<const SpecUtils::Measurement> &long_background = src_info.long_background;
    //src_info.src_no_poisson;
    //src_info.background_no_poisson;

    const vector<ExpectedPhotopeakInfo> &expected_photopeaks = info.expected_photopeaks;


    ReferenceLineInfo ref_line_info;
    ref_line_info.m_validity = ReferenceLineInfo::InputValidity::Valid;
    ref_line_info.m_source_type = ReferenceLineInfo::SourceType::OneOffSrcLines;
    ref_line_info.m_has_coincidences = false;

    double max_peak_area = 0.0;
    for( const ExpectedPhotopeakInfo &line : expected_photopeaks )
    {
      if( line.nsigma_over_background < 1.0 ) //Dont show 1-sigma peaks
        continue;

      ReferenceLineInfo::RefLine ref_line;
      ref_line.m_energy = line.effective_energy;
      ref_line.m_normalized_intensity = line.peak_area;
      max_peak_area = std::max( max_peak_area, line.peak_area );
      ref_line.m_drf_factor = 1.0;
      ref_line.m_shield_atten = 1.0;
      ref_line.m_particle_sf_applied = 1.0;
      ref_line.m_color = Wt::WColor( Wt::GlobalColor::darkBlue );
      ref_line.m_decay_intensity = line.peak_area;
      ref_line.m_particle_type = ReferenceLineInfo::RefLine::Particle::Gamma;
      ref_line.m_source_type = ReferenceLineInfo::RefLine::RefGammaType::Normal;
      ref_line.m_attenuation_applies = false;

      ref_line_info.m_ref_lines.push_back( ref_line );
    }//for( const ExpectedPhotopeakInfo &line : expected_photopeaks )

    if( max_peak_area > 1.0 )
    {
      for( auto &p : ref_line_info.m_ref_lines )
        p.m_normalized_intensity /= max_peak_area;
    }

    string ref_line_json;
    ref_line_info.toJson( ref_line_json );

    std::string title = "";
    std::string dataTitle = "";
    bool useLogYAxis = true, showVerticalGridLines = false, showHorizontalGridLines = false;
    bool legendEnabled = true, compactXAxis = true;
    bool showPeakUserLabels = false, showPeakEnergyLabels = false, showPeakNuclideLabels = false, showPeakNuclideEnergyLabels = false;
    bool showEscapePeakMarker = false, showComptonPeakMarker = false, showComptonEdgeMarker = false, showSumPeakMarker = false;
    bool backgroundSubtract = false;
    float xMin = 0, xMax = 3000;
    std::map<std::string,std::string> refernce_lines_json;
    refernce_lines_json["TruthPeaks"] = ref_line_json;

    D3SpectrumExport::D3SpectrumChartOptions options( title, "Energy (keV)", "Counts/Channel",
                                                     dataTitle, useLogYAxis,
                                                     showVerticalGridLines, showHorizontalGridLines,
                                                     legendEnabled, compactXAxis,
                                                     showPeakUserLabels, showPeakEnergyLabels, showPeakNuclideLabels,
                                                     showPeakNuclideEnergyLabels, showEscapePeakMarker, showComptonPeakMarker,
                                                     showComptonEdgeMarker, showSumPeakMarker, backgroundSubtract,
                                                     xMin, xMax, refernce_lines_json );

    D3SpectrumExport::D3SpectrumOptions foreground_opts, background_opts;
    foreground_opts.line_color = "black";
    background_opts.line_color = "steelblue";
    foreground_opts.title = src_info.src_name;
    background_opts.title = "Background";
    foreground_opts.display_scale_factor = 1.0;
    background_opts.display_scale_factor = spectrum->live_time() / long_background->live_time();
    foreground_opts.spectrum_type = SpecUtils::SpectrumType::Foreground;
    background_opts.spectrum_type = SpecUtils::SpectrumType::Background;

    const string div_id = "chart_" + std::to_string(i);


    output << "<fieldset style=\"\">" << endl
    << "<legend>" << src_info.src_name << "</legend>" << endl;

    output << "<div id=\"" << div_id << "\" class=\"chart\" oncontextmenu=\"return false;\"></div>" << endl;  // Adding the main chart div


    output << "<script>" << endl;

    D3SpectrumExport::write_js_for_chart( output, div_id, options.m_dataTitle, options.m_xAxisTitle, options.m_yAxisTitle );

    std::vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
    measurements.emplace_back( spectrum.get(), foreground_opts );
    measurements.emplace_back( long_background.get(), background_opts );

    write_and_set_data_for_chart( output, div_id, measurements );


    output << R"delim(
    const resizeChart)delim" << i << R"delim( = function(){
      let height = window.innerHeight;
      let width = document.documentElement.clientWidth;
      let el = spec_chart_)delim" << div_id << R"delim(.chart;
      el.style.width = 0.8*width + "px";
      el.style.height = Math.min(500,Math.max(250, Math.min(0.4*width,height-175))) + "px";
      el.style.marginLeft = 0.05*width + "px";
      el.style.marginRight = 0.05*width + "px";
      
    )delim"
    << "  spec_chart_" << div_id << R"delim(.handleResize();
    };
    
    window.addEventListener('resize', resizeChart)delim" << i << R"delim();
    )delim" << endl;

    write_set_options_for_chart( output, div_id, options );

    output << "spec_chart_" << div_id << ".setReferenceLines( reference_lines_" << div_id << " );" << endl;

    output << "resizeChart" << i << "();" << endl;
    output << "</script>" << endl;

    vector<ExpectedPhotopeakInfo> photopeaks = expected_photopeaks;
    std::sort( begin(photopeaks), end(photopeaks), [](const ExpectedPhotopeakInfo &lhs, const ExpectedPhotopeakInfo &rhs ){
      return lhs.peak_area > rhs.peak_area;
    });

    output << "<table class=\"TopLinesTable\" style=\"\">" << endl;
    output << "<tr><th>Energy (keV)</th><th>Truth Area</th><th>Truth CPS</th><th>ROI Lower</th><th>ROI Upper</th><th>Continuum Area</th></tr>" << endl;
    output << "<caption>Top gamma lines in spectrum</caption>" << endl;
    for( size_t i = 0; (i < photopeaks.size()) && (i < 20); ++i )
    {
      const ExpectedPhotopeakInfo &p = photopeaks[i];
      output << "<tr>"
      << "<td>" << p.effective_energy << "</td>"
      << "<td>" << p.peak_area << "</td>"
      << "<td>" << p.peak_area / spectrum->live_time() << "</td>"
      << "<td>" << p.roi_lower << "</td>"
      << "<td>" << p.roi_upper << "</td>"
      << "<td>" << p.continuum_area << "</td>"
      << "</tr>"
      << endl;
    }

    if( photopeaks.size() > 20 )
      output << "<tr><td colspan=\"6\">Plus " << (photopeaks.size() - 20) << " more peaks</td></tr>" << endl;

    output << "</table>" << endl;

    output << "</fieldset>" << endl;
  }//for( size_t i = 0; i < src_info.size(); ++i )

  output << "</body>" << endl;
  output << "</html>" << endl;
}//void write_html_summary( const vector<DataSrcInfo> &src_info )
#endif //WRITE_ALL_SPEC_TO_HTML

}//namespace PeakFitImproveData
