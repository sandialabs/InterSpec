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

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;

namespace BatchPeak
{

void fit_energy_cal_from_fit_peaks( shared_ptr<SpecUtils::Measurement> &raw, vector<PeakDef> peaks, const size_t num_coefs )
{
  if( !raw || raw->num_gamma_channels() < 16 )
    throw runtime_error( "update_gain_from_peak: invalid input spectrum" );
  
  vector<EnergyCal::RecalPeakInfo> peakinfos;
  
  shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = raw->energy_calibration();
  assert( orig_cal && orig_cal->valid() && (orig_cal->coefficients().size() > 1) );
  
  for( const PeakDef &peak : peaks )
  {
    if( !peak.hasSourceGammaAssigned() )
    {
      cerr << "Warning: peak at " << peak.mean() << " keV doesnt have a source assigned, so will"
      << " not be used for energy calibration" << endl;
      continue;;
    }
    
    EnergyCal::RecalPeakInfo recalpeak;
    recalpeak.peakMean = peak.mean();
    recalpeak.peakMeanUncert = peak.meanUncert();
    recalpeak.peakMeanBinNumber = orig_cal->channel_for_energy( peak.mean() );
    recalpeak.photopeakEnergy = peak.gammaParticleEnergy();
    
    peakinfos.push_back( recalpeak );
  }//for( const double peak_energy : peak_energies )
  
  const size_t nchannels = raw->num_gamma_channels();
  
  const vector<pair<float,float>> &dev_pairs = orig_cal->deviation_pairs();
  
  vector<float> fit_coefs_uncert;
  vector<float> fit_coefs = orig_cal->coefficients();
  
  vector<bool> fitfor( orig_cal->coefficients().size(), false );
  if( fitfor.size() < num_coefs )
    fitfor.resize( num_coefs );
  for( size_t i = 0; i < num_coefs; ++i )
    fitfor[i] = true;
  
  shared_ptr<SpecUtils::EnergyCalibration> updated_cal = make_shared<SpecUtils::EnergyCalibration>();
  
  switch( orig_cal->type() )
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
}//void fit_energy_cal_from_fit_peaks(...)

  
void fit_peaks_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          const BatchPeakFitOptions &options )
{
  vector<string> warnings;
  
  if( files.empty() )
    throw runtime_error( "No input files specified." );
  
  if( !options.output_dir.empty() && !SpecUtils::is_directory(options.output_dir) )
    throw runtime_error( "Output directory ('" + options.output_dir + "'), is not a directory." );
    
  if( options.write_n42_with_peaks && options.output_dir.empty() )
    throw runtime_error( "If you specify to write N42 files with the fit peaks, you must specify an output directory." );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    const char *message = "Unable to load the nuclide decay database."
    " Run the executable from a working directory where the 'data' folder is "
    " located - the program expects to load data/sandia.decay.xml.";
    
    throw runtime_error( message );
  }//if( !db )
  
  std::shared_ptr<const SpecMeas> cached_exemplar_n42;
  for( const string filename : files )
  {
    const BatchPeakFitResult fit_results
                 = fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, options );
    
    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.exemplar;
    warnings.insert(end(warnings), begin(fit_results.warnings), end(fit_results.warnings) );
    
    if( !fit_results.success )
      continue;
    
    assert( fit_results.measurement );
    
    if( options.write_n42_with_peaks && fit_results.measurement )
    {
      string outn42 = SpecUtils::append_path(options.output_dir, SpecUtils::filename(filename) );
      if( !SpecUtils::iequals_ascii(SpecUtils::file_extension(filename), ".n42") )
        outn42 += ".n42";
      
      if( SpecUtils::is_file( outn42 ) )
      {
        warnings.push_back( "Not writing '" + outn42 + "', as it would overwrite a file.");
      }else
      {
        if( !fit_results.measurement->save2012N42File( outn42 ) )
          warnings.push_back( "Failed to write '" + outn42 + "'.");
        else
          cout << "Have written '" << outn42 << "' with peaks" << endl;
      }
    }//if( options.write_n42_with_peaks )
    
    if( !options.output_dir.empty() )
    {
      const string leaf_name = SpecUtils::filename(filename);
      string outcsv = SpecUtils::append_path(options.output_dir, leaf_name) + ".CSV";
      
      if( SpecUtils::is_file( outcsv ) )
      {
        warnings.push_back( "Not writing '" + outcsv + "', as it would overwrite a file.");
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
          PeakModel::write_peak_csv( output_csv, leaf_name, fit_results.fit_peaks, fit_results.spectrum );
          cout << "Have written '" << outcsv << "'" << endl;
        }
      }//if( SpecUtils::is_file( outcsv ) ) / else
    }//if( !options.output_dir.empty() )
    
    if( options.to_stdout )
    {
      const string leaf_name = SpecUtils::filename(filename);
      cout << "peaks for '" << leaf_name << "':" << endl;
      PeakModel::write_peak_csv( cout, leaf_name, fit_results.fit_peaks, fit_results.spectrum );
      cout << endl;
    }
  }//for( const string filename : files )
    
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;
  
  return EXIT_SUCCESS;
}//fit_peaks_in_files(...)
  
  

BatchPeak::BatchPeakFitResult fit_peaks_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
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
        
    if( !exemplar_is_n42 && options.use_exemplar_energy_cal )
      throw runtime_error( "Exemplar file wasnt an N42 file, but using its energy cal was specified - not allowed." );
    
    if( exemplar_is_n42 )
      exemplar_n42 = exemplar;
  }//if( !cached_exemplar_n42 )
  
  std::shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> exemplar_peaks;
  if( exemplar_n42 )
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
    
    assert( exemplar_peaks );
    assert( exemplar_spectrum );
    if( !exemplar_peaks || !exemplar_spectrum )
      throw logic_error( "Logic error retrieving spectrum or peaks from exemplar." );
    
    
    
    if( options.use_exemplar_energy_cal
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
    shared_ptr<SpecMeas> specfile = make_shared<SpecMeas>();
    const bool loaded = specfile->load_file( filename, SpecUtils::ParserType::Auto, filename );
    if( !loaded || !specfile->num_measurements() )
    {
      results.warnings.push_back( "Couldnt read in '" + filename + "' as a spectrum file -- skipping." );
      
      return results;
    }//if( !loaded )
    
    results.measurement = specfile;
    
    set<int> used_sample_nums;
    const vector<string> det_names = specfile->detector_names();
    
    shared_ptr<SpecUtils::Measurement> spec;
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
    
    assert( spec );
    if( !spec )
    {
      results.warnings.push_back( "Spectrum file '" + filename + "' failed to extract wanted spectrum." );
      return results;
    }
    
    if( options.use_exemplar_energy_cal )
    {
      assert( exemplar_spectrum );
      shared_ptr<const SpecUtils::EnergyCalibration> energy_cal = exemplar_spectrum->energy_calibration();
      
      const size_t num_spec_chan = spec->num_gamma_channels();
      if( energy_cal->num_channels() == num_spec_chan )
      {
        spec->set_energy_calibration( energy_cal );
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
              new_cal.reset();
              results.warnings.push_back( "Not using exemplar energy calibration for '" + filename + "'"
                                         " - its lower energy channel calibration, and exemplar has"
                                         " fewer channels." );
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
        
        if( new_cal )
          spec->set_energy_calibration( new_cal );
      }//if( num channels match exemplar ) / else
    }//if( options.use_exemplar_energy_cal )
    
    
    if( !spec || (spec->num_gamma_channels() < 64)
       || (spec->energy_calibration_model() == SpecUtils::EnergyCalType::InvalidEquationType) )
    {
      results.warnings.push_back( "Failed to get spectrum from file '" + filename + "' -- skipping." );
      return results;
    }
    
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
    for( const auto &p : starting_peaks )
      results.exemplar_peaks.push_back( make_shared<PeakDef>(p) );
    
    //  We are re-using functions called by the GUI InterSpec, so there are some extra arguments
    //  that arent totally applicable to us right now.
    const double lower_energy = starting_peaks.front().mean() - 0.1;
    const double uppper_energy = starting_peaks.back().mean() + 0.1;
    
    // We are not refitting peaks, because the areas may be wildly different.
    const bool isRefit = false;
    
    //Use default for peak fit filters
    const double ncausalitysigma = 0.0, stat_threshold  = 0.0, hypothesis_threshold = 0.0;
    
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
      candidate_peaks.push_back( peak );
    }//for( const auto &p : exemplar_peaks )
    
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
                                                energy_cal_peaks, spec, {}, isRefit );
        fit_energy_cal_from_fit_peaks( spec, peaks, 4 );
      }//for( size_t i = 0; i < 1; ++i )
    }//if( options.refit_energy_cal )
    
    vector<PeakDef> fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                                stat_threshold, hypothesis_threshold,
                                                candidate_peaks, spec, {}, isRefit );
    
    // Could re-fit the peaks again...
    fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                stat_threshold, hypothesis_threshold,
                                fit_peaks, spec, {}, true );
    
    //cout << "Fit for the following " << fit_peaks.size() << " peaks (the exemplar file had "
    //<< starting_peaks.size() <<  ") from the raw spectrum:"
    //<< endl;
    
    for( PeakDef &p: fit_peaks )
    {
      // Find nearest exemplar peak, and we'll use this to set nuclides, colors, and such
      shared_ptr<PeakDef> exemplar_parent;
      for( const auto &exemplar : starting_peaks )
      {
        const double fit_mean = p.mean();
        const double energy_diff = fabs( fit_mean - exemplar.mean() );
        // We will require the fit peak to be within 0.25 FWHM (arbitrarily chosen distance)
        //  of the exemplar peak, and we will use the exemplar peak closest in energy to the
        //  fit peak
        if( (energy_diff < 0.5*p.fwhm())
           && (!exemplar_parent || (energy_diff < fabs(exemplar_parent->mean() - fit_mean))) )
        {
          exemplar_parent = make_shared<PeakDef>( exemplar );
        }
      }//for( loop over exemplars to find original peak corresponding to the fit peak 'p' )
      
      if( exemplar_parent )
        p.inheritUserSelectedOptions( *exemplar_parent, false );
      else
        results.warnings.push_back( "In '" + filename + "', failed to find exemplar peak"
          " corresponding to peak fit with mean=" + std::to_string( p.mean() ) + " keV." );
      //cout << "\tmean=" << p.mean() << ", FWHM=" << p.fwhm() << ", area=" << p.amplitude() << endl;
    }
    
    //cout << endl << endl;
    
    // Now we will associate the fit peaks to the spectrum and save an N42 file you can open up in
    //  InterSpec and inspect the fits.
    
    deque<shared_ptr<const PeakDef>> fit_peaks_ptrs;
    for( const auto &p: fit_peaks )
      fit_peaks_ptrs.push_back( make_shared<const PeakDef>(p) );
    
    specfile->setPeaks( fit_peaks_ptrs, used_sample_nums );
   
    results.success = true;
    results.measurement = specfile;
    results.spectrum = spec;
    results.sample_numbers = used_sample_nums;
    results.fit_peaks = fit_peaks_ptrs;
  }
  
  return results;
}//fit_peaks_in_file(...)
  
}//namespace BatchPeak

