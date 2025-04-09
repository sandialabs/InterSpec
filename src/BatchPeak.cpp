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

#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <exception>

#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/BatchInfoLog.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;


namespace BatchPeak
{

void propagate_energy_cal( const shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                           shared_ptr<SpecUtils::Measurement> &to_spectrum,
                           shared_ptr<SpecMeas> &to_specfile,
                           const set<int> &used_sample_nums )
{
  assert( to_spectrum );
  assert( energy_cal );
  assert( to_specfile );
  
  shared_ptr<const SpecUtils::EnergyCalibration> original_cal = to_spectrum->energy_calibration();
  assert( original_cal );
  if( !energy_cal )
    throw runtime_error( "Missing energy in from spectrum." );
  
  const size_t num_spec_chan = to_spectrum->num_gamma_channels();
  shared_ptr<const SpecUtils::EnergyCalibration> updated_cal;
  if( energy_cal == original_cal )
  {
    // We already have this energy cal, nothing to do here
    updated_cal = energy_cal;
  }else
  {
    if( energy_cal->num_channels() == num_spec_chan )
    {
      to_spectrum->set_energy_calibration( energy_cal );
      updated_cal = energy_cal;
    }else
    {
      auto new_cal = make_shared<SpecUtils::EnergyCalibration>();
      
      switch( energy_cal->type() )
      {
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          new_cal->set_polynomial( num_spec_chan, energy_cal->coefficients(), energy_cal->deviation_pairs() );
          break;
          
        case SpecUtils::EnergyCalType::FullRangeFraction:
          new_cal->set_full_range_fraction( num_spec_chan, energy_cal->coefficients(), energy_cal->deviation_pairs() );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
          if( num_spec_chan > energy_cal->num_channels() )
          {
            throw std::runtime_error( " its lower energy channel calibration, and exemplar has"
                                     " fewer channels." );
            new_cal.reset();
          }else
          {
            vector<float> channel_energies = *new_cal->channel_energies();
            channel_energies.resize( num_spec_chan + 1 );
            new_cal->set_lower_channel_energy( num_spec_chan, channel_energies );
          }
          break;
          
        case SpecUtils::EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
      }//switch( energy_cal->type() )
      
      updated_cal = new_cal;
      if( new_cal )
        to_spectrum->set_energy_calibration( new_cal );
    }//if( num channels match exemplar ) / else
  }//if( energy_cal == original_cal )
  
  // Update `to_specfile` as well, as we may write it back out as a N42
  if( updated_cal && to_specfile )
  {
    const vector<string> &det_names = to_specfile->detector_names();
    
    if( used_sample_nums.size() > 1 )
    {
      // Translate peaks from old energy, to new energy
      auto peaks = to_specfile->peaks( used_sample_nums );
      if( peaks && peaks->size() && (original_cal != updated_cal) )
      {
        const deque<shared_ptr<const PeakDef>> new_peaks
               = EnergyCal::translatePeaksForCalibrationChange( *peaks, original_cal, updated_cal );
        to_specfile->setPeaks( new_peaks, used_sample_nums );
      }
    }//if( used_sample_nums.size() > 1 )
    
    for( const int sample : used_sample_nums )
    {
      shared_ptr<const SpecUtils::EnergyCalibration> prev_cal;
      for( const string &det : det_names )
      {
        auto m = to_specfile->measurement( sample, det );
        if( m && (m->num_gamma_channels() == updated_cal->num_channels()) )
        {
          prev_cal = m->energy_calibration();
          if( prev_cal != updated_cal )
            to_specfile->set_energy_calibration( updated_cal, m );
        }
      }//for( const string &det : det_names )
      
      if( prev_cal && (used_sample_nums.size() > 1) && (prev_cal != updated_cal) )
      {
        // Translate peaks from old energy, to new energy, only if we havent already done it
        auto peaks = to_specfile->peaks( {sample} );
        if( peaks && peaks->size() )
        {
          const deque<shared_ptr<const PeakDef>> new_peaks
          = EnergyCal::translatePeaksForCalibrationChange( *peaks, prev_cal, updated_cal );
          to_specfile->setPeaks( new_peaks, {sample} );
        }//if( peaks && peaks->size() )
      }//if( prev_cal && (used_sample_nums.size() > 1) )
    }//for( const int sample : used_sample_nums )
  }//if( updated_cal )
}//propagate_energy_cal(...)
  
void fit_energy_cal_from_fit_peaks( shared_ptr<SpecUtils::Measurement> &raw, vector<PeakDef> peaks )
{
  if( !raw || raw->num_gamma_channels() < 16 )
    throw runtime_error( "update_gain_from_peak: invalid input spectrum" );
  
  vector<EnergyCal::RecalPeakInfo> peakinfos;
  
  shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = raw->energy_calibration();
  assert( orig_cal && orig_cal->valid() && (orig_cal->coefficients().size() > 1) );
  
  const SpecUtils::EnergyCalType energy_cal_type = (orig_cal && orig_cal->valid()) 
                                                    ? orig_cal->type()
                                                    : SpecUtils::EnergyCalType::InvalidEquationType;
  switch( energy_cal_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
      if( orig_cal->coefficients().size() < 2 )
        throw std::logic_error( "Somehow the energy calibration has less than two coefficients." );
      break;
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      throw std::logic_error( "The original calibration must be either polynomial or full-range-fraction." );
      break;
  }//switch( exemplar_cal->type() )
  
  
  for( const PeakDef &peak : peaks )
  {
    if( !peak.hasSourceGammaAssigned() )
    {
      cerr << "Warning: peak at " << peak.mean() << " keV doesnt have a source assigned, so will"
      << " not be used for energy calibration" << endl;
      continue;
    }
    
    
    if( !peak.useForEnergyCalibration() )  //This defaults to true, so if this is false, then user selected it
      continue;
    
    EnergyCal::RecalPeakInfo recalpeak;
    recalpeak.peakMean = peak.mean();
    recalpeak.peakMeanUncert = peak.meanUncert();
    recalpeak.peakMeanBinNumber = orig_cal->channel_for_energy( peak.mean() );
    recalpeak.photopeakEnergy = peak.gammaParticleEnergy();
    
    peakinfos.push_back( recalpeak );
  }//for( const double peak_energy : peak_energies )
  
  std::sort( begin(peakinfos), end(peakinfos), 
    []( const EnergyCal::RecalPeakInfo &lhs, const EnergyCal::RecalPeakInfo &rhs ) -> bool {
    return lhs.peakMean < rhs.peakMean;
  } );
  
  if( peakinfos.empty() )
    throw runtime_error( "No peaks selected for use in energy calibration." );
  
  const size_t nchannels = raw->num_gamma_channels();
  
  const vector<float> &orig_coefs = orig_cal->coefficients();
  const vector<pair<float,float>> &dev_pairs = orig_cal->deviation_pairs();
  const size_t num_coefs = orig_coefs.size();
  assert( num_coefs >= 2 );
  
  vector<float> fit_coefs_uncert;
  vector<float> fit_coefs = orig_coefs;
  vector<bool> fitfor( num_coefs, false );
  
  // If we only have a couple peaks, then we cant fit for like 4 coefficients; we'll
  //  just hardcode a rough heuristic for what coefficients to fit for, because I think
  //  bothering the user with this level of detail is probably a bit too much.
  //
  //  But also note that we have already fit peaks, so we must be reasonably close
  //  right now anyway, so we can be a bit more liberal in terms of fitting more
  //  coefficients than a user might normally do during an interactive session.
  //
  //  We dont want to update both gain and offset from like the two Co60 peaks, so
  //  we'll count the effective number of peaks, requiring them to be separated a bit.
  //  We will require the separation to be max( 100, min(0.1*total-energy-range, 200)) keV,
  //  with the 0.1, 100 and 200, all being entirely arbitrary, but hopefully reasonable.
  size_t num_effective_peaks = 1;
  const double energy_range = (orig_cal->upper_energy() - orig_cal->lower_energy());
  const double sep_dist = std::max( 100.0, std::min(0.1*energy_range, 200.0) );
  size_t last_peak = 0;
  for( size_t peak_index = 1; peak_index < peakinfos.size(); ++peak_index )
  {
    const double delta_energy = peakinfos[peak_index].peakMean - peakinfos[last_peak].peakMean;
    assert( delta_energy >= 0.0 );
    if( delta_energy >= sep_dist )
    {
      num_effective_peaks += 1;
      last_peak = peak_index;
    }
  }//for( loop over peaks to count how many are separated by at least `sep_dist` )
  
  if( num_effective_peaks == 1 )
  {
    fitfor[1] = true; // Only fit gain
  }else
  {
    for( size_t index = 0; (index < num_effective_peaks) && (index < fitfor.size()); ++index )
      fitfor[index] = true;
    // TODO: we could consider not fitting for offset if `peakinfos[0].peakMean` is less than 
    //       ~200 keV or something, but for the moment we wont, because we know we are close
    //       in energy calibration, and we know the peaks the user is interested in, so overfitting
    //       isnt as large of a concern as not lining up the ROI edges as much
  }
  
  shared_ptr<SpecUtils::EnergyCalibration> updated_cal = make_shared<SpecUtils::EnergyCalibration>();
  
  switch( energy_cal_type )
  {
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      EnergyCal::fit_energy_cal_poly( peakinfos, fitfor, nchannels, dev_pairs, fit_coefs, fit_coefs_uncert );
      updated_cal->set_polynomial( nchannels, fit_coefs, dev_pairs );
      break;
      
    case SpecUtils::EnergyCalType::FullRangeFraction:
      EnergyCal::fit_energy_cal_frf( peakinfos, fitfor, nchannels, dev_pairs, fit_coefs, fit_coefs_uncert );
      updated_cal->set_full_range_fraction( nchannels, fit_coefs, dev_pairs );
      break;
      
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      assert( 0 );
      break;
  }//switch( exemplar_cal->type() )
  
  raw->set_energy_calibration( updated_cal );
  
#if( PERFORM_DEVELOPER_CHECKS )
  cout << "Updated energy calibration using ROIs from exemplar.\n\tCoefficients:\n";
  assert( fit_coefs.size() == orig_cal->coefficients().size() );
  for( size_t i = 0; i < fit_coefs.size(); ++i )
    cout << "\t\t" << std::setprecision(6) << std::setw(12) << orig_cal->coefficients().at(i)
    << "\t-->\t" << std::setprecision(6) << std::setw(12) << fit_coefs[i] << endl;
  
  cout << "This moved peak energies:" << endl;
  for( const auto &peak : peaks )
  {
    const double energy = peak.mean();
    const double orig = orig_cal->channel_for_energy( energy );
    const double now = updated_cal->channel_for_energy( energy );
    cout << "\t" << std::setprecision(6) << std::setw(12) << energy << " keV from channel "
    << std::setprecision(6) << std::setw(12) << orig << " to "
    << std::setprecision(6) << std::setw(12) << now << endl;
  }
  
  cout << endl << endl;
#endif
}//void fit_energy_cal_from_fit_peaks(...)

  
void fit_peaks_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                        const ::BatchPeak::BatchPeakFitOptions &options )
{
  vector<string> warnings;
  
  if( files.empty() )
    throw runtime_error( "No input files specified." );
  
  if( !options.output_dir.empty() && !SpecUtils::is_directory(options.output_dir) )
    throw runtime_error( "Output directory ('" + options.output_dir + "'), is not a directory." );
    
  if( options.write_n42_with_results && options.output_dir.empty() )
    throw runtime_error( "If you specify to write N42 files with the fit peaks, you must specify an output directory." );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    const char *message = "Unable to load the nuclide decay database."
    " Run the executable from a working directory where the 'data' folder is "
    " located - the program expects to load data/sandia.decay.xml.";
    
    throw runtime_error( message );
  }//if( !db )
  
  const vector<pair<string,string>> spec_chart_js_and_css = BatchInfoLog::load_spectrum_chart_js_and_css();
  
  // Load report templates, and setup inja::environment
  inja::Environment env = BatchInfoLog::get_default_inja_env( options );
  
  
  nlohmann::json summary_json;
  
  BatchInfoLog::add_exe_info_to_json( summary_json );
  BatchInfoLog::add_peak_fit_options_to_json( summary_json, options );
  
  summary_json["ExemplarFile"] = exemplar_filename;
  if( !exemplar_sample_nums.empty() )
    summary_json["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
  summary_json["InputFiles"] = files;
  for( const pair<string,string> &key_val : spec_chart_js_and_css )
    summary_json[key_val.first] = key_val.second;
  
  std::shared_ptr<const SpecMeas> cached_exemplar_n42;
  for( const string filename : files )
  {
    const BatchPeak::BatchPeakFitResult fit_results
                 = fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, nullptr, {}, options );
    
    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.exemplar;
    warnings.insert(end(warnings), begin(fit_results.warnings), end(fit_results.warnings) );
    
    nlohmann::json data;
    BatchInfoLog::add_exe_info_to_json( data );
    BatchInfoLog::add_peak_fit_options_to_json( data, options );
    data["ExemplarFile"] = exemplar_filename;
    if( !exemplar_sample_nums.empty() )
      data["ExemplarSampleNumbers"] = vector<int>{begin(exemplar_sample_nums), end(exemplar_sample_nums)};
    data["Filepath"] = filename;
    data["Filename"] = SpecUtils::filename( filename );
    data["ParentDir"] = SpecUtils::parent_path( filename );
    
    BatchInfoLog::add_peak_fit_results_to_json( data, fit_results );
    
    summary_json["Files"].push_back( data );
    
    for( const pair<string,string> &key_val : spec_chart_js_and_css )
      data[key_val.first] = key_val.second;
    
    for( string tmplt : options.report_templates )
    {
      try
      {
        const string rpt = BatchInfoLog::render_template( tmplt, env,
                            BatchInfoLog::TemplateRenderType::PeakFitIndividual, options, data );
        
        if( options.to_stdout && !SpecUtils::iequals_ascii(tmplt, "html" ) )
          cout << "\n\n" << rpt << endl << endl;
        
        if( !options.output_dir.empty() )
        {
          const string out_file
                    = BatchInfoLog::suggested_output_report_filename( filename, tmplt,
                                  BatchInfoLog::TemplateRenderType::PeakFitIndividual, options );
          
          if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
          {
            warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                               " See the '--overwrite-output-files' option to force writing." );
          }else
          {
#ifdef _WIN32
            const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
            std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
            std::ofstream output( out_file.c_str(), ios::binary | ios::out);
#endif
            if( !output )
              warnings.push_back( "Failed to open report output '" + out_file + "', for writing.");
            else
              output.write( rpt.c_str(), rpt.size() );
          }//if( is file ) / else write file
        }//if( !options.output_dir.empty() )
      }catch( inja::InjaError &e )
      {
        const string msg = "Error templating results (" + e.type + ": line "
        + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
        + " of '" + tmplt + "'): " + e.message + ". While processing '" + filename + "'.";
        
        cerr << msg << endl;
        warnings.push_back( msg );
      }catch( std::exception &e )
      {
        cerr << "Error templating results: " << e.what() << endl;
        warnings.push_back( "Error templating results: " + string(e.what()) );
      }
    }//for( const string &tmplt : options.report_templates )
    
    if( !options.output_dir.empty() && options.create_json_output )
      BatchInfoLog::write_json( options, warnings, filename, data );
    
    if( !fit_results.success )
      continue;
    
    assert( fit_results.measurement );
    
    if( options.write_n42_with_results && fit_results.measurement )
    {
      string outn42 = SpecUtils::append_path(options.output_dir, SpecUtils::filename(filename) );
      if( !SpecUtils::iequals_ascii(SpecUtils::file_extension(filename), ".n42") )
        outn42 += ".n42";
      
      if( SpecUtils::is_file( outn42 ) && !options.overwrite_output_files )
      {
        warnings.push_back( "Not writing '" + outn42 + "', as it would overwrite a file."
                           " See the '--overwrite-output-files' option to force writing." );
      }else
      {
        if( !fit_results.measurement->save2012N42File( outn42 ) )
          warnings.push_back( "Failed to write '" + outn42 + "'.");
        else
          cout << "Have written '" << outn42 << "' with peaks" << endl;
      }
    }//if( options.write_n42_with_peaks )
    
    deque<shared_ptr<const PeakDef>> fit_peaks = fit_results.fit_peaks;
    if( options.show_nonfit_peaks )
    {
      for( const auto p : fit_results.unfit_exemplar_peaks )
      {
        auto np = make_shared<PeakDef>(*p);
        np->setAmplitude( 0.0 );
        np->setAmplitudeUncert( 0.0 );
        np->setSigmaUncert( 0.0 );
        auto cont = make_shared<PeakContinuum>( *np->continuum() );
        const size_t num_cont_pars = PeakContinuum::num_parameters(cont->type());
        for( size_t i = 0; i < num_cont_pars; ++i )
        {
          cont->setPolynomialCoef( i, 0.0 );
          cont->setPolynomialUncert( i, -1.0 );
        }
        np->setContinuum( cont );
        np->set_coefficient( -1.0, PeakDef::Chi2DOF );
        fit_peaks.push_back(np);
      }
      std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMeanShrdPtr );
    }//if( make_nonfit_peaks_zero )
    
    
    if( !options.output_dir.empty() && options.create_csv_output )
    {
      const string leaf_name = SpecUtils::filename(filename);
      string outcsv = SpecUtils::append_path(options.output_dir, leaf_name) + ".CSV";
      
      if( SpecUtils::is_file(outcsv) && !options.overwrite_output_files )
      {
        warnings.push_back( "Not writing '" + outcsv + "', as it would overwrite a file."
                           " See the '--overwrite-output-files' option to force writing." );
      }else
      {
#ifdef _WIN32
        const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(outcsv);
        std::ofstream output_csv( woutcsv.c_str() );
#else
        std::ofstream output_csv( outcsv.c_str() );
#endif
        
        if( !output_csv )
        {
          warnings.push_back( "Failed to open '" + outcsv + "', for writing.");
        }else
        {
          PeakModel::write_peak_csv( output_csv, leaf_name, 
                                    PeakModel::PeakCsvType::Full, fit_peaks, fit_results.spectrum );
          cout << "Have written '" << outcsv << "'" << endl;
        }
      }//if( SpecUtils::is_file( outcsv ) ) / else
    }//if( !options.output_dir.empty() )
    
    if( options.to_stdout )
    {
      const string leaf_name = SpecUtils::filename(filename);
      cout << "peaks for '" << leaf_name << "':" << endl;
      PeakModel::write_peak_csv( cout, leaf_name, PeakModel::PeakCsvType::Full,
                                fit_peaks, fit_results.spectrum );
      cout << endl;
    }
  }//for( const string filename : files )
    
  // Add any encountered errors to output summary JSON
  for( const string &warn : warnings )
    summary_json["Warnings"].push_back( warn );
  
  // Now write summary report(s)
  for( const string &summary_tmplt : options.summary_report_templates )
  {
    try
    {
      const string rpt = BatchInfoLog::render_template( summary_tmplt, env,
                       BatchInfoLog::TemplateRenderType::PeakFitSummary, options, summary_json );
      
      if( options.to_stdout && !SpecUtils::iequals_ascii(summary_tmplt, "html" ) )
        cout << "\n\n" << rpt << endl << endl;
      
      if( !options.output_dir.empty() )
      {
        const string out_file
                    = BatchInfoLog::suggested_output_report_filename( "", summary_tmplt,
                                    BatchInfoLog::TemplateRenderType::PeakFitSummary, options );
        
        if( SpecUtils::is_file(out_file) && !options.overwrite_output_files )
        {
          warnings.push_back( "Not writing '" + out_file + "', as it would overwrite a file."
                             " See the '--overwrite-output-files' option to force writing." );
        }else
        {
#ifdef _WIN32
          const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(out_file);
          std::ofstream output( woutcsv.c_str(), ios::binary | ios::out );
#else
          std::ofstream output( out_file.c_str(), ios::binary | ios::out );
#endif
          if( !output )
            throw runtime_error( "Failed to open summary peak fit report output, '" + out_file + "'" );
          
          output.write( rpt.c_str(), rpt.size() );
        }//if( file exists and dont overwrite ) / else
      }
    }catch( inja::InjaError &e )
    {
      const string msg = "Error templating summary peak fit output (" + e.type + ": line "
      + std::to_string(e.location.line) + ", column " + std::to_string(e.location.column)
      + " of '" + summary_tmplt + "'): " + e.message + ".";
      
      cerr << msg << endl;
      warnings.push_back( msg );
    }catch( std::exception &e )
    {
      warnings.push_back( "Error making summary peak fit output: " + string(e.what()) );
    }
  }//if( !options.summary_report_template.empty() )
  
  
  if( !options.output_dir.empty() && options.create_json_output )
    BatchInfoLog::write_json( options, warnings, "", summary_json );
  
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;
}//fit_peaks_in_files(...)
  
  
void get_exemplar_spectrum_and_peaks(
            std::shared_ptr<const SpecUtils::Measurement> &exemplar_spectrum,
            std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &exemplar_peaks,
            std::set<int> &exemplar_sample_nums,
            const std::shared_ptr<const SpecMeas> &exemplar_n42 )
{
  const vector<string> det_names = exemplar_n42->detector_names();
  const set<set<int>> withPeakSampleNums = exemplar_n42->sampleNumsWithPeaks();
    
  if( exemplar_n42->measurements().empty() )
  {
    throw logic_error( "Exemplar spectrum file did not have any measurements." );
  }else if( !exemplar_sample_nums.empty() )
  {
    const set<set<int>> withPeakSampleNums = exemplar_n42->sampleNumsWithPeaks();
    if( withPeakSampleNums.count(exemplar_sample_nums) )
      throw runtime_error( "The specified exemplar sample numbers did not have peaks associated with them." );
      
    exemplar_peaks = exemplar_n42->peaks( exemplar_sample_nums );
    exemplar_spectrum = exemplar_n42->sum_measurements( exemplar_sample_nums, det_names, nullptr );
      
    if( !exemplar_peaks || exemplar_peaks->empty() || !exemplar_spectrum )
      throw runtime_error( "The specified exemplar sample numbers did not have peaks, or spectra couldnt be summed." );
  }else if( exemplar_n42->measurements().size() == 1 )
  {
    exemplar_spectrum = exemplar_n42->measurements().front();
    if( !exemplar_spectrum )
      throw logic_error( "Unexpected invalid exemplar spectrum." );
      
    exemplar_sample_nums.insert( exemplar_spectrum->sample_number() );
    exemplar_peaks = exemplar_n42->peaks( exemplar_sample_nums );
    if( !exemplar_peaks || exemplar_peaks->empty() )
      throw logic_error( "Exemplar spectrum did not contain any peaks." );
  }else
  {
    set<set<int>> foregroundPeaks, backgroundPeaks, otherPeaks;
    for( const set<int> &samples : withPeakSampleNums )
    {
      auto peaks = exemplar_n42->peaks( samples );
      if( !peaks || peaks->empty() )
        continue;
        
      auto m = exemplar_n42->sum_measurements( samples, det_names, nullptr );
      if( !m )
        continue;
      
      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity:
        case SpecUtils::SourceType::Calibration:
          otherPeaks.insert( samples );
          break;
          
        case SpecUtils::SourceType::Background:
          backgroundPeaks.insert( samples );
          break;
          
        case SpecUtils::SourceType::Foreground:
        case SpecUtils::SourceType::Unknown:
          foregroundPeaks.insert( samples );
          break;
      }//switch( m->source_type() )
    }//for( const set<int> &samples : withPeakSampleNums )
    
    if( foregroundPeaks.size() > 1 )
      throw runtime_error( "Ambiguous which peaks to use from exemplar file" );
    
    if( foregroundPeaks.empty() && backgroundPeaks.empty() )
      throw runtime_error( "No valid peaks exemplar file"
                          + string(otherPeaks.empty() ? "." : " (intrinsic and calibration spectra in files are ignored).") );
    
    if( foregroundPeaks.empty() && (backgroundPeaks.size() != 1) )
      throw runtime_error( "Ambiguous which peaks to use from exemplar file; multiple background spectra with peaks." );
    
    if( foregroundPeaks.size() == 1 )
    {
      exemplar_sample_nums = *begin(foregroundPeaks);
    }else if( backgroundPeaks.size() == 1 )
    {
      exemplar_sample_nums = *begin(backgroundPeaks);
    }else
    {
      throw logic_error( "Error getting peaks from exemplar." );
    }
    
    exemplar_peaks = exemplar_n42->peaks( exemplar_sample_nums );
    exemplar_spectrum = exemplar_n42->sum_measurements( exemplar_sample_nums, det_names, nullptr );
  }//if( sample nums specified ) / else ( single meas ) / else ( search for peaks )
}//void get_exemplar_spectrum_and_peaks(...)
  

BatchPeak::BatchPeakFitResult fit_peaks_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          std::shared_ptr<SpecMeas> cached_spectrum,
                          std::set<int> foreground_sample_numbers,
                          const BatchPeakFitOptions &options )
{
  shared_ptr<const SpecMeas> exemplar_n42 = cached_exemplar_n42;
  
  if( !exemplar_n42 )
  {
    // Read in the N42 file exported from InterSpec.
    //  This file should have good energy calibration applied, have the peaks fit that you care
    //  about, and have exactly one spectrum
    auto exemplar = make_shared<SpecMeas>();
    const bool exemplar_is_n42 = exemplar->load_N42_file( exemplar_filename );
        
    if( !exemplar_is_n42 && (options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background) )
      throw runtime_error( "Exemplar file wasnt an N42 file, but using its energy cal was specified - not allowed." );
    
    if( exemplar_is_n42 )
      exemplar_n42 = exemplar;
  }//if( !cached_exemplar_n42 )
  
  std::shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> exemplar_peaks;
  if( exemplar_n42 )
  {
    get_exemplar_spectrum_and_peaks( exemplar_spectrum, exemplar_peaks,
                                    exemplar_sample_nums, exemplar_n42 );
    
    assert( exemplar_peaks );
    assert( exemplar_spectrum );
    if( !exemplar_peaks || !exemplar_spectrum )
      throw logic_error( "Logic error retrieving spectrum from exemplar." );
    
    if( !exemplar_peaks || !exemplar_spectrum )
      throw logic_error( "Logic error retrieving peaks from exemplar." );
    
    if( (options.use_exemplar_energy_cal || options.use_exemplar_energy_cal_for_background)
       && (exemplar_spectrum->energy_calibration_model() == SpecUtils::EnergyCalType::InvalidEquationType) )
      throw runtime_error( "Exemplar spectrum doesnt have a valid energy calibration." );
    
    if( exemplar_spectrum->num_gamma_channels() < 64  )
      throw runtime_error( "Exemplar spectrum doesnt have enough gamma channels." );
  }//if( exemplar_is_n42 )
  
  assert( !exemplar_n42 || exemplar_spectrum );
  assert( !exemplar_n42 || (exemplar_peaks && !exemplar_peaks->empty()) );
  
  
  BatchPeakFitResult results;
  results.file_path = exemplar_filename;
  results.options = options;
  results.exemplar = exemplar_n42;
  results.exemplar_sample_nums = exemplar_sample_nums;
  //if( exemplar_peaks )
  //  results.exemplar_peaks = *exemplar_peaks;
  results.exemplar_spectrum = exemplar_spectrum;
  results.success = false;
  
  {
    shared_ptr<SpecMeas> specfile;
    
    if( cached_spectrum )
    {
      specfile = cached_spectrum;
    }else
    {
      specfile = make_shared<SpecMeas>();
      const bool loaded = specfile->load_file( filename, SpecUtils::ParserType::Auto, filename );
      if( !loaded || !specfile->num_measurements() )
      {
        results.warnings.push_back( "Couldnt read in '" + filename + "' as a spectrum file -- skipping." );
        
        return results;
      }//if( !loaded )
    }//if( cached_spectrum ) / else
    
    results.measurement = specfile;
    
    shared_ptr<SpecUtils::Measurement> spec;
    const vector<string> det_names = specfile->detector_names();
    
    set<int> used_sample_nums;
    if( !foreground_sample_numbers.empty() )
    {
      try
      {
        spec = specfile->sum_measurements( foreground_sample_numbers, det_names, nullptr );
        used_sample_nums = foreground_sample_numbers;
      }catch( std::exception &e )
      {
        results.warnings.push_back( "Invalid sample numbers specified to sum: " + string(e.what()) );
        results.success = false;
        return results;
      }
    }else
    {
      if( specfile->sample_numbers().size() == 1 )
      {
        used_sample_nums = specfile->sample_numbers();
        spec = specfile->sum_measurements( used_sample_nums, det_names, nullptr );
      }else
      {
        set<int> foregroundSamples, backgroundSamples, unknownSamples, otherSamples;
        for( const int sample_num : specfile->sample_numbers() )
        {
          for( const string det_name : det_names )
          {
            auto m = specfile->measurement( sample_num, det_name );
            if( !m )
              continue;
            switch( m->source_type() )
            {
              case SpecUtils::SourceType::IntrinsicActivity:
              case SpecUtils::SourceType::Calibration:
                otherSamples.insert( sample_num );
                break;
              case SpecUtils::SourceType::Background:
                backgroundSamples.insert( sample_num );
                break;
              case SpecUtils::SourceType::Foreground:
                foregroundSamples.insert( sample_num );
                break;
              case SpecUtils::SourceType::Unknown:
                unknownSamples.insert( sample_num );
                break;
            }//switch( m->source_type() )
          }//for( const string det_name : det_names )
        }//for( const int sample_num : specfile.sample_numbers() )
        
        if( foregroundSamples.size() == 1 )
        {
          used_sample_nums = foregroundSamples;
        }else if( unknownSamples.size() == 1 )
        {
          used_sample_nums = unknownSamples;
        }else if( backgroundSamples.size() == 1 )
        {
          used_sample_nums = backgroundSamples;
        }else if( otherSamples.size() == 1 )
        {
          used_sample_nums = otherSamples;
        }else
        {
          results.warnings.push_back( "Spectrum file '" + filename + "' was ambiguous of which spectrum to use for peak fitting." );
          return results;
        }
        
        spec = specfile->sum_measurements( used_sample_nums, det_names, nullptr );
      }//if( specfile.sample_numbers().size() == 1 ) / else
    }//if( !foreground_sample_numbers.empty() ) / else
    
    assert( spec );
    if( !spec )
    {
      results.warnings.push_back( "Spectrum file '" + filename + "' failed to extract wanted spectrum." );
      return results;
    }
    
    
    if( options.use_exemplar_energy_cal )
    {
      try
      {
        propagate_energy_cal( exemplar_spectrum->energy_calibration(), spec, specfile, used_sample_nums );
      }catch( std::exception &e )
      {
        results.warnings.push_back( "Not using exemplar energy calibration for '" + filename + "': "
                                   + std::string(e.what()) );
      }
    }//if( options.use_exemplar_energy_cal )
    
    
    if( !spec || (spec->num_gamma_channels() < 64)
       || (spec->energy_calibration_model() == SpecUtils::EnergyCalType::InvalidEquationType) )
    {
      results.warnings.push_back( "Failed to get spectrum from file '" + filename + "' -- skipping." );
      return results;
    }
    
    set<int> back_sample_nums;
    shared_ptr<SpecMeas> background_n42;
    if( !options.background_subtract_file.empty() )
    {
      background_n42 = make_shared<SpecMeas>();
      if( !background_n42->load_file( options.background_subtract_file, SpecUtils::ParserType::Auto ) )
        throw runtime_error( "Couldnt open background file '" + options.background_subtract_file + "'" );
      
      if( options.background_subtract_samples.empty() )
      {
        back_sample_nums = background_n42->sample_numbers();
        if( back_sample_nums.size() != 1 )
          throw runtime_error( "There should only be a single sample in background subtract file." );
      }else
      {
        back_sample_nums = options.background_subtract_samples;
      }
      
      if( background_n42->num_measurements() == 1 )
      {
        assert( !back_sample_nums.empty() );
        if( !back_sample_nums.empty()
           && ((*begin(back_sample_nums)) != background_n42->measurements()[0]->sample_number() ) )
        {
          results.warnings.push_back( "Specified background sample number invalid." );
          return results;
        }
        
        results.background = make_shared<SpecUtils::Measurement>( *(background_n42->measurements()[0]) );
      }else
      {
        try
        {
          const vector<string> &dets = background_n42->detector_names();
          results.background = background_n42->sum_measurements( back_sample_nums, dets, nullptr );
        }catch( std::exception &e )
        {
          results.warnings.push_back( "Failed to sum spectrum from background '"
                                     + options.background_subtract_file + "' -- skipping '"
                                     + filename + "'." );
          return results;
        }//try / catch to sum background
        
        if( !results.background->energy_calibration()
           || !results.background->energy_calibration()->valid()
           || (results.background->num_gamma_channels() < 16)
           || (results.background->live_time() <= 1.0E-3) )
        {
          results.warnings.push_back( "Background '"
                                     + options.background_subtract_file
                                     + "' didnt have energy calibration, too few channels, or no live-time"
                                     " -- skipping '" + filename + "'." );
          return results;
        }//
        
        if( options.use_exemplar_energy_cal_for_background )
        {
          try
          {
            propagate_energy_cal( exemplar_spectrum->energy_calibration(), results.background, background_n42, back_sample_nums );
          }catch( std::exception &e )
          {
            results.warnings.push_back( "Not using exemplar energy calibration for background of '" + filename + "': "
                                       + std::string(e.what()) );
          }
        }//if( options.use_exemplar_energy_cal_for_background )
      }//if( background->num_measurements() == 1 ) / else
      
      
      try
      {
        const bool no_neg = true;
        const bool do_round = false;
        
        const bool sf = spec->live_time() / results.background->live_time();
        
        shared_ptr<const vector<float>> fore_counts = spec->gamma_counts();
        shared_ptr<const vector<float>> back_counts = results.background->gamma_counts();
        
        // Make sure back_counts has the same energy calibration and fore_counts, so we can subtract
        //  on a bin-by-bin basis
        if( results.background->energy_calibration() != spec->energy_calibration()
           && (*results.background->energy_calibration()) != (*spec->energy_calibration()) )
        {
          auto new_backchan = make_shared<vector<float>>( fore_counts->size(), 0.0f );
          SpecUtils::rebin_by_lower_edge( *results.background->channel_energies(), *back_counts,
                                         *spec->channel_energies(), *new_backchan );
          back_counts = new_backchan;
        }
        
        // Create what will be the background subtracted foreground
        auto back_sub_counts = make_shared<vector<float>>( *fore_counts );
        
        //back_counts and fore_counts should always be the same size, but we'll be a bit verbose anyway
        assert( back_counts->size() == fore_counts->size() );
        const size_t nchann = std::min( back_counts->size(), fore_counts->size() );
        
        // Do the actual background subtraction
        for( size_t i = 0; i < nchann; ++i )
        {
          float &val = (*back_sub_counts)[i];
          val -= sf*(*back_counts)[i];
          
          if( no_neg )
            val = std::max( 0.0f, val );
          
          if( do_round )
            val = std::round( val );
        }//for( size_t i = 0; i < nchann; ++i )
        
        // Create a new Measurement object, based on the old foreground
        auto newspec = make_shared<SpecUtils::Measurement>( *spec );
        newspec->set_gamma_counts( back_sub_counts, spec->live_time(), spec->real_time() );
        vector<string> remarks = spec->remarks();
        remarks.push_back( "This spectrum has been background subtracted in InterSpec" );
        newspec->set_remarks( remarks );
        newspec->set_sample_number( 1 );
        
        spec = newspec;
        
        results.measurement->removeAllPeaks();
        results.measurement->remove_measurements( results.measurement->measurements() );
        results.measurement->add_measurement( spec, true );
      }catch( std::exception &e )
      {
        results.warnings.push_back( "Error background subtracting: '" + string(e.what())
                                   + "' -- skipping '" + filename + "'." );
        return results;
      }//try / catch
    }//if( !options.background_subtract_file.empty() )
    
    
    results.sample_numbers = used_sample_nums;
    results.spectrum = spec;
    
    vector<PeakDef> starting_peaks;
    
    if( exemplar_peaks )
    {
      for( const auto p : *exemplar_peaks )
        starting_peaks.push_back( *p );
    }else
    {
      // Open the input CSV file - and get those peaks.
      assert( !exemplar_n42 );
      
#ifdef _WIN32
      const std::wstring wfilename = SpecUtils::convert_from_utf8_to_utf16(exemplar_filename);
      std::ifstream input( wfilename.c_str() );
#else
      std::ifstream input( exemplar_filename.c_str() );
#endif
      
      if( !input.is_open() )
        throw runtime_error( "Exemplar peak file, '" + exemplar_filename + "', could not be opened." );
      
      try
      {
        starting_peaks = PeakModel::csv_to_candidate_fit_peaks( spec, input );
      }catch( std::exception &e )
      {
        throw runtime_error( "Invalid input exemplar CSV peak file: " + string(e.what()) );
      }
      
      if( starting_peaks.empty() )
      {
        results.warnings.push_back( "No candidate peaks for file '" + filename + "', perhaps peaks from CSV were not good candidate peaks -- skipping file." );
        
        return results;
      }
      
    }//if( !exemplar_peaks )
    
    if( starting_peaks.empty() )
    {
      results.warnings.push_back( "No candidate peaks for '" + filename + "' - maybe a programming logic error?" );
      
      return results;
    }
    
    // Sort the peaks by mean even though its probably already sorted
    std::sort( begin(starting_peaks), end(starting_peaks), &PeakDef::lessThanByMean );
    
    results.exemplar_peaks.clear();
    vector<shared_ptr<const PeakDef>> exemplar_peaks;
    set<shared_ptr<const PeakDef>> unused_exemplar_peaks;
    for( const auto &p : starting_peaks )
    {
      auto ep = make_shared<PeakDef>(p);
      results.exemplar_peaks.push_back( ep );
      exemplar_peaks.push_back( ep );
      unused_exemplar_peaks.insert( ep );
      results.unfit_exemplar_peaks.push_back( ep ); //We will clear this later on if unsuccessful
    }
    
    //  We are re-using functions called by the GUI InterSpec, so there are some extra arguments
    //  that arent totally applicable to us right now.
    const double lower_energy = starting_peaks.front().mean() - 0.1;
    const double uppper_energy = starting_peaks.back().mean() + 0.1;
    
    // We are not refitting peaks, because the areas may be wildly different.
    const bool isRefit = false;
    
    const double ncausalitysigma = 0.0;
    const double stat_threshold = options.peak_stat_threshold;
    const double hypothesis_threshold = options.peak_hypothesis_threshold;
    
    vector<PeakDef> candidate_peaks;
    for( const auto &p : starting_peaks )
    {
      PeakDef peak = p;
      
      // If you had selected for certain peak quantities to not be fit for in exemplar spectrum
      //  in InterSpec, then by default those quantities also wouldnt be fit for here.
      //  You can alter this similar to below.
      //  For the moment we'll just make sure the peak amplitude will be fit for.
      
      //peak.setFitFor( PeakDef::Mean, true );
      //peak.setFitFor( PeakDef::Sigma, false );
      peak.setFitFor( PeakDef::GaussAmplitude, true );
      
      
      // Lets also make sure starting amplitude is something halfway reasonable,
      //  and continuum coefficients are reasonable starting values
      if( peak.gausPeak() )
      {
        std::shared_ptr<PeakContinuum> continuum = peak.continuum();
        assert( continuum );
        if( continuum && continuum->isPolynomial() )
        {
          const PeakContinuum::OffsetType origType = continuum->type();
          continuum->calc_linear_continuum_eqn( spec, peak.mean(), peak.lowerX(), peak.upperX(), 2, 2 );
          continuum->setType( origType );
        }//if( continuum )
        
        const double mean = peak.mean(), fwhm = peak.fwhm();
        const double data_area = spec->gamma_integral( mean - fwhm, mean + fwhm );
        
        if( (data_area > 1) && (peak.amplitude() > data_area) )
        {
          double cont_area = continuum->offset_integral(  mean - fwhm, mean + fwhm, spec );
          if( (cont_area > 0.0) && (cont_area < data_area) )
          {
            peak.setAmplitude( data_area - cont_area );
          }else
          {
            peak.setAmplitude( 0.25*data_area );
          }
          
        }//if( exemplar peak is clearly much larger than data )
      }//if( peak.gausPeak() )
      
      candidate_peaks.push_back( peak );
    }//for( const auto &p : exemplar_peaks )
    
    const bool isHPGe = PeakFitUtils::is_high_res( spec );
    results.original_energy_cal = spec ? spec->energy_calibration() : nullptr;
    
    if( options.refit_energy_cal )
    {
      // We will refit the energy calibration - maybe a few times - to really hone in on things
      for( size_t i = 0; i < 5; ++i )
      {
        vector<PeakDef> energy_cal_peaks = candidate_peaks;
        for( auto &peak : energy_cal_peaks )
        {
          peak.setFitFor( PeakDef::Mean, true );
          peak.setFitFor( PeakDef::Sigma, true );
          peak.setFitFor( PeakDef::GaussAmplitude, true );
        }
        
        vector<PeakDef> peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                                stat_threshold, hypothesis_threshold,
                                                energy_cal_peaks, spec, {}, isRefit, isHPGe );
        try
        {
          fit_energy_cal_from_fit_peaks( spec, peaks );
        }catch( std::exception &e )
        {
          const string msg = "Energy calibration not performed: " + string(e.what());
          auto pos = std::find( begin(results.warnings), end(results.warnings), msg );
          if( pos == end(results.warnings) )
            results.warnings.push_back( msg );
        }
      }//for( size_t i = 0; i < 1; ++i )
      
      // Propagate the updated energy cal to the result file
      assert( spec && spec->energy_calibration() && spec->energy_calibration()->valid() );
      shared_ptr<const SpecUtils::EnergyCalibration> new_cal = spec ? spec->energy_calibration() : nullptr;
      if( new_cal && new_cal->valid() )
      {
        try
        {
          propagate_energy_cal( new_cal, spec, results.measurement, {} );
        }catch( std::exception &e )
        {
          results.warnings.push_back( "Failed to propagate fit energy calibration in '" + filename + "'." );
        }
      }else
      {
        results.warnings.push_back( "Failed to fit an appropriate energy calibration in '" + filename + "'." );
      }
      
      results.refit_energy_cal = spec ? spec->energy_calibration() : nullptr;
    }//if( options.refit_energy_cal )
    
    vector<PeakDef> fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                                stat_threshold, hypothesis_threshold,
                                                candidate_peaks, spec, {}, isRefit, isHPGe );
    
    // Re-fit the peaks again a few more times
    for( size_t i = 0; i < 3; ++i )
    {
      fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                  stat_threshold, hypothesis_threshold,
                                  fit_peaks, spec, {}, true, isHPGe );
    }
    
    //cout << "Fit for the following " << fit_peaks.size() << " peaks (the exemplar file had "
    //<< starting_peaks.size() <<  ") from the raw spectrum:"
    //<< endl;
    
    for( PeakDef &p: fit_peaks )
    {
      // Find nearest exemplar peak, and we'll use this to set nuclides, colors, and such
      const double fit_mean = p.mean();
      
      shared_ptr<const PeakDef> exemplar_parent;
      for( const auto &exemplar : exemplar_peaks )
      {
        const double exemplar_mean = exemplar->mean();
        
        const double energy_diff = fabs( fit_mean - exemplar_mean );
        // We will require the fit peak to be within 0.5 FWHM (arbitrarily chosen distance)
        //  of the exemplar peak, and we will use the exemplar peak closest in energy to the
        //  fit peak
        if( ((energy_diff < 0.5*p.fwhm()) || (energy_diff < 0.5*exemplar->fwhm()))
           && (!exemplar_parent || (energy_diff < fabs(exemplar_parent->mean() - fit_mean))) )
        {
          exemplar_parent = exemplar;
        }
      }//for( loop over exemplars to find original peak corresponding to the fit peak 'p' )
      
      if( exemplar_parent )
      {
        unused_exemplar_peaks.erase( exemplar_parent );
        p.inheritUserSelectedOptions( *exemplar_parent, false );
      }else
      {
        results.warnings.push_back( "In '" + filename + "', failed to find exemplar peak"
                                   " corresponding to peak fit with mean=" + std::to_string( p.mean() ) + " keV." );
        //cout << "\tmean=" << p.mean() << ", FWHM=" << p.fwhm() << ", area=" << p.amplitude() << endl;
      }//if( exemplar_parent ) / else
    }//for( PeakDef &p: fit_peaks )
    
    //cout << endl << endl;
    
    // Now we will associate the fit peaks to the spectrum and save an N42 file you can open up in
    //  InterSpec and inspect the fits.
    
    deque<shared_ptr<const PeakDef>> fit_peaks_ptrs;
    for( const auto &p: fit_peaks )
      fit_peaks_ptrs.push_back( make_shared<const PeakDef>(p) );
    
    if( options.background_subtract_file.empty() )
    {
      specfile->setPeaks( fit_peaks_ptrs, used_sample_nums );
    }else
    {
      assert( specfile->num_measurements() == 1 );
      specfile->setPeaks( fit_peaks_ptrs, specfile->sample_numbers() );
    }
   
    results.success = true;
    results.measurement = specfile;
    results.spectrum = spec;
    results.sample_numbers = used_sample_nums;
    results.fit_peaks = fit_peaks_ptrs;
    results.unfit_exemplar_peaks.clear();
    results.unfit_exemplar_peaks.insert( end(results.unfit_exemplar_peaks),
                                        begin(unused_exemplar_peaks), end(unused_exemplar_peaks) );
    std::sort( begin(results.unfit_exemplar_peaks), end(results.unfit_exemplar_peaks),
              &PeakDef::lessThanByMeanShrdPtr );
  }
  
  return results;
}//fit_peaks_in_file(...)
  
}//namespace BatchPeak

