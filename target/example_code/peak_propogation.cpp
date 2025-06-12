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
#include "InterSpec_config.h"

#include <deque>
#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SandiaDecay/SandiaDecay.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

// Forward declarations
#ifdef _WIN32
void getUtf8Args( int &argc, char ** &argv );
#endif

/** Fits a single peak near energy \p peak_energy, and uses that to update the energy gain.
 */
void update_gain_from_peak( SpecMeas &specfile, const vector<float> &peak_energies );


/** Fits the `num_coefs` leading energy calibration coefficients, using peaks you have already fit to data.
 */
void fit_energy_cal_from_fit_peaks( SpecMeas &specfile, vector<PeakDef> peaks, const size_t num_coefs );

/** Example program that reads in an N42 file exported from InterSpec where you have performed energy calibration and fit for
 peaks of interest - this is the "exemplar" file.  This program then reads in a "raw" spectrum, applies the energy calibration from the
 exemplar file, and tries to fit the same peaks as from the other file, and then saves those peaks to a N42 file you can open in InterSpec,
 as well as to a CSV file like you would export from InterSpec.
 
 This program also has a that tries to update the linear energy calibration based on the 609.31 keV Ra226 peak, which is meant as
 a minor correction - you can comment out this step if you dont want or need it.
 
 If instead you wanted to read in candidate peaks from a CSV file exported from InterSpec, you could use
 the PeakMode::csv_to_candidate_fit_peaks(...) function to read in the candidate peaks, instead of getting them from the exemplar
 N42 file.
 */
int main( int argc, char **argv )
{
#ifdef _WIN32
  getUtf8Args( argc, argv );
#endif

  //  You need to supply exactly three arguments to this program
  if( argc != 4 )
  {
    cout << "Program usage: " << argv[0]
         << " /path/to/cal_spectrum_with_fit_peaks.n42"
         << " /path/to/test_spectrum.spe"
         << " /path/to/output_peaks"
         << endl;
    
    return EXIT_FAILURE;
  }//if( argc != 4 )
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
  {
    cerr << "\n\nUnable to load the decay database - nuclides will not be able to be assigned"
            " to peaks.  Run the executable from a working directory where the 'data' folder is "
            " located - the program expects to load data/sandia.decay.xml" << endl;
  }
  
  
  // The "exemplar" spectrum file is the file you exported from InterSpec as a N42-2012 file, that
  //  has exactly one gamma spectrum in it, with some fit peaks.
  const string exemplar_n42_filename = argv[1];
   
  // The "raw" spectrum file is the measurement from your detector you want to fit the same peaks
  //  as the exemplar file has in it.  It can be any format spectrum file (SPE, CSV, PCF, N42, etc),
  //  but it must have exactly one gamma spectrum, with the same number of channels as the exemplar.
  const string raw_spec_filename = argv[2];
   
  // The output will be a CSV of the peaks fit, as well as a N42-2012 spectrum file with the peaks
  //  fit as well, so you can load it into InterSpec easily to check things.
  const string output_filename = argv[3];
  
  
  try
  {
    // Sanity check we wont overwrite any of the input files
    const string output_n42_filename = output_filename + ".n42";
    const string output_csv_filename = output_filename + ".CSV";
    
    if( (output_n42_filename == raw_spec_filename )
        || (output_n42_filename == exemplar_n42_filename)
        || (output_csv_filename == raw_spec_filename)
        || (output_csv_filename == exemplar_n42_filename) )
      throw runtime_error( "The output name would cause overwriting one of the input files." );
    
    
    // Read in the N42 file exported from InterSpec.
    //  This file should have good energy calibration applied, have the peaks fit that you care
    //  about, and have exactly one spectrum
  
    SpecMeas exemplar_n42;
    if( !exemplar_n42.load_N42_file( exemplar_n42_filename ) )
      throw runtime_error( "Unable to parse file '" + exemplar_n42_filename + "'" );
    
    if( exemplar_n42.measurements().size() != 1 )
      throw runtime_error( "Expected exactly one spectrum in '" + exemplar_n42_filename + "'" );
    
    const shared_ptr<const SpecUtils::Measurement> exemplar = exemplar_n42.measurements()[0];
    assert( exemplar );
    
    if( exemplar->num_gamma_channels() < 128 )
      throw runtime_error( "'" + exemplar_n42_filename + "' doesnt contain a gamma spectrum" );
    
    shared_ptr<const SpecUtils::EnergyCalibration> exemplar_cal = exemplar->energy_calibration();
    
    // Only supporting polynomial or FRF energy calibration types right now
    switch( exemplar_cal->type() )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
        break;
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "'" + exemplar_n42_filename
                             + "' must be have either polynomial or full range fraction energy calibration" );
        break;
    }//switch( exemplar_cal->type() )
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> exemplar_peaks_ptr = exemplar_n42.peaks( {exemplar->sample_number()} );
    
    if( !exemplar_peaks_ptr || exemplar_peaks_ptr->empty() )
      throw runtime_error( "'" + exemplar_n42_filename + "' didnt contain any fit peaks" );
    
    deque<shared_ptr<const PeakDef>> exemplar_peaks = *exemplar_peaks_ptr; //for convience
    
    // Sort the deque by peak mean even though its probably already sorted
    std::sort( begin(exemplar_peaks), end(exemplar_peaks), &PeakDef::lessThanByMeanShrdPtr );
    
    // If we're here, the input exemplar spectrum file has a gamma spectrum, valid energy
    //  calibration, and some fit peaks.
    
    // Lets open the "raw" spectrum file
    SpecMeas raw_n42;
    if( !raw_n42.load_file( raw_spec_filename, SpecUtils::ParserType::Auto ) )
      throw runtime_error( "Unable to parse file '" + raw_spec_filename + "'" );
    
    if( raw_n42.measurements().size() != 1 )
      throw runtime_error( "Expected exactly one spectrum in '" + raw_spec_filename + "'" );
    
    const shared_ptr<const SpecUtils::Measurement> raw = raw_n42.measurements()[0];
    assert( raw );
    
    if( exemplar->num_gamma_channels() != raw->num_gamma_channels() )
      throw runtime_error( "The spectrum files have different numbers of channels" );
    
    // If we are here the "raw" file is probably compatible with the exemplar
    
    // Apply the energy calibration from the exemplar to the raw
    raw_n42.set_energy_calibration( exemplar_cal, raw );
    
    
    // Lets update the energy using the 609 keV peak
    update_gain_from_peak( raw_n42, {609.31} );
    
    
    // Define inputs to the peak fit.
    
    //  We are re-using functions called by the GUI InterSpec, so there are some extra arguments
    //  that arent totally applicable to us right now.
    const double lower_energy = exemplar_peaks.front()->mean() - 0.1;
    const double uppper_energy = exemplar_peaks.back()->mean() + 0.1;
    
    // We are not refitting peaks, because the areas may be wildly different.
    const bool isRefit = false;
    
    //Use default for peak fit filters
    const double ncausalitysigma = 0.0, stat_threshold  = 0.0, hypothesis_threshold = 0.0;
    
    vector<PeakDef> candidate_peaks;
    for( const auto &p : exemplar_peaks )
    {
      PeakDef peak = *p;
      
      // If you had selected for certain peak quantities to not be fit for in exemplar spectrum
      //  in InterSpec, then by default those quantities also wouldnt be fit for here.
      //  You can alter this similar to below.
      //  For the moment we'll just make sure the peak amplitude will be fit for.
      
      //peak.setFitFor( PeakDef::Mean, true );
      //peak.setFitFor( PeakDef::Sigma, false );
      peak.setFitFor( PeakDef::GaussAmplitude, true );
      
      const SandiaDecay::Nuclide * const nuc = peak.parentNuclide();
      cout << "Peak at " << peak.mean() << " keV has nuclide: " << (nuc ? nuc->symbol : string("null")) << endl;
      
      candidate_peaks.push_back( peak );
    }//for( const auto &p : exemplar_peaks )
    
    
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
                                              energy_cal_peaks, raw, {}, isRefit );
      fit_energy_cal_from_fit_peaks( raw_n42, peaks, 4 );
    }//for( size_t i = 0; i < 1; ++i )
    
    
    vector<PeakDef> fit_peaks = fitPeaksInRange( lower_energy, uppper_energy, ncausalitysigma,
                                                 stat_threshold, hypothesis_threshold,
                                                 candidate_peaks, raw, {}, isRefit );
    
    cout << "Fit for the following " << fit_peaks.size() << " peaks (the exemplar file had "
         << exemplar_peaks.size() <<  ") from the raw spectrum:"
         << endl;
    
    for( PeakDef &p: fit_peaks )
    {
      // Find nearest exemplar peak, and we'll use this to set nuclides, colors, and such
      shared_ptr<const PeakDef> exemplar_parent;
      for( const auto &exemplar : exemplar_peaks )
      {
        const double fit_mean = p.mean();
        const double energy_diff = fabs( fit_mean - exemplar->mean() );
        // We will require the fit peak to be within 0.25 FWHM (arbitrarily chosen distance) 
        //  of the exemplar peak, and we will use the exemplar peak closest in energy to the
        //  fit peak
        if( (energy_diff < 0.25*p.fwhm()) 
            && (!exemplar_parent 
                 || (energy_diff < fabs(exemplar_parent->mean() - fit_mean))) ) 
        {
          exemplar_parent = exemplar;
        }
      }//for( loop over exemplars to find original peak cooresponding to the fit peak 'p' )

      if( exemplar_parent )
        p.inheritUserSelectedOptions( *exemplar_parent, false );
      else 
        cerr << "Failed to find exemplar peak cooresponding to peak fit with mean=" << p.mean() 
             << " - will not set parent nuclide/color/etc of fit peak." << endl;

      cout << "\tmean=" << p.mean() << ", FWHM=" << p.fwhm() << ", area=" << p.amplitude() << endl;
    }
    
    cout << endl << endl;
    
    // Now we will associate the fit peaks to the spectrum and save an N42 file you can open up in
    //  InterSpec and inspect the fits.
    
    deque<shared_ptr<const PeakDef>> fit_peaks_ptrs;
    for( const auto &p: fit_peaks )
      fit_peaks_ptrs.push_back( make_shared<const PeakDef>(p) );
    
    raw_n42.setPeaks( fit_peaks_ptrs, {raw->sample_number()} );
    
    if( !raw_n42.save2012N42File( output_n42_filename ) )
      throw runtime_error( "Failed to write output N42 file" );
    
    cout << "Have written '" << output_n42_filename << "'" << endl;
    
    ofstream output_csv( output_csv_filename.c_str(), ios::out | ios::binary );
    if( !output_csv )
      throw runtime_error( "Failed to open output CSV file '"
                           + output_csv_filename + "' for writing" );
    
    string rawfilename = SpecUtils::filename(raw_spec_filename);
    
    PeakModel::write_peak_csv( output_csv, rawfilename, fit_peaks_ptrs, raw );
    
    cout << "Have written '" << output_csv_filename << "'" << endl;
  }catch( std::exception &e )
  {
    cerr << "Encountered an error: " << e.what() << endl;
    return EXIT_FAILURE;
  }//try / catch
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


void fit_energy_cal_from_fit_peaks( SpecMeas &specfile, vector<PeakDef> peaks, const size_t num_coefs )
{
  if( specfile.measurements().size() != 1 )
    throw runtime_error( "update_gain_from_peak: expected one spectrum in file" );
  
  vector<EnergyCal::RecalPeakInfo> peakinfos;
  
  auto raw = specfile.measurements()[0];
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
  
  specfile.set_energy_calibration( updated_cal, raw );
  
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




void update_gain_from_peak( SpecMeas &specfile, const vector<float> &peak_energies )
{
  if( specfile.measurements().size() != 1 )
    throw runtime_error( "update_gain_from_peak: expected one spectrum in file" );
  
  vector<EnergyCal::RecalPeakInfo> peakinfos;
  
  auto raw = specfile.measurements()[0];
  shared_ptr<const SpecUtils::EnergyCalibration> orig_cal = raw->energy_calibration();
  assert( orig_cal && orig_cal->valid() && (orig_cal->coefficients().size() > 1) );
  
  
  for( const double peak_energy : peak_energies )
  {
    const double pixelPerKev = 2;
    vector<shared_ptr<const PeakDef>> fit_peaks = searchForPeakFromUser( peak_energy, pixelPerKev, raw, {} ).first;
    
    if( fit_peaks.empty() )
      throw runtime_error( "Could not fit a peaks near " + std::to_string(peak_energy) + " for energy calibration" );
    
    if( fit_peaks.size() != 1 )
      throw runtime_error( "Logic error fitting energy calibration peak." );
    
    auto peak = fit_peaks[0];
    
    cout << "Peak fit for energy calibration has mean=" << peak->mean()
    << ", FWHM=" << peak->fwhm() << ", area=" << peak->peakArea()
    << endl;
  
    // We could easily adjust "gain" by just taking ratios, but we'll invoke the full energy
    //  calibration fitting mechanism for demonstration, but also incase we want to be able to use
    //  multiple peaks later.
    
    EnergyCal::RecalPeakInfo recalpeak;
    recalpeak.peakMean = peak->mean();
    recalpeak.peakMeanUncert = peak->meanUncert();
    recalpeak.peakMeanBinNumber = orig_cal->channel_for_energy( peak->mean() );
    recalpeak.photopeakEnergy = peak_energy;
    
    peakinfos.push_back( recalpeak );
  }//for( const double peak_energy : peak_energies )
    
  const size_t nchannels = raw->num_gamma_channels();
  
  const vector<pair<float,float>> &dev_pairs = orig_cal->deviation_pairs();
  
  vector<float> fit_coefs_uncert;
  vector<float> fit_coefs = orig_cal->coefficients();
  vector<bool> fitfor( orig_cal->coefficients().size(), false );
  fitfor[1] = true; // only "fit" for gain
  
  
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
  
  specfile.set_energy_calibration( updated_cal, raw );

  cout << "Updated energy calibration coefficients:\n";
  assert( fit_coefs.size() == orig_cal->coefficients().size() );
  for( size_t i = 0; i < fit_coefs.size(); ++i )
    cout << "\t" << std::setprecision(10) << orig_cal->coefficients().at(i) << "\t-->\t" << fit_coefs[i] << endl;
  
  cout << endl << endl;
}//void update_gain_from_peak( ... )

/*
void refit_peaks(  )
{
  try
  {
    PeakModel * const model = interspec->peakModel();
    const shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const shared_ptr<const SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground);
    
    
    if( !model || !data || !foreground ) //shouldnt ever happen
    {
      passMessage( "No data loaded to refit", "", WarningWidget::WarningMsgInfo );
      return;
    }
    
    //shared_ptr<const PeakDef> peak = interspec->nearestPeak( rightClickEnergy );
    shared_ptr<const PeakDef> peak = model->nearestPeak( rightClickEnergy );
    if( !peak )
    {
      passMessage( "There was no peak to refit", "", WarningWidget::WarningMsgInfo );
      return;
    }
    
    vector<PeakDef> inputPeak, fixedPeaks, outputPeak;
    const vector<shared_ptr<const PeakDef>> peaksInRoi = model->peaksSharingRoi( peak );
    const vector<shared_ptr<const PeakDef>> peaksNotInRoi  = model->peaksNotSharingRoi( peak );
    
    assert( peaksInRoi.size() >= 1 );
    
    for( const auto &m : peaksInRoi )
      inputPeak.push_back( *m );
    
    for( const auto &m : peaksNotInRoi )
      fixedPeaks.push_back( *m );
    
    std::sort( inputPeak.begin(), inputPeak.end(), &PeakDef::lessThanByMean );
    
    
    
    if( inputPeak.size() > 1 )
    {
      const shared_ptr<const DetectorPeakResponse> &detector = foreground->detector();
      const PeakShrdVec result = refitPeaksThatShareROI( data, detector, peaksInRoi, PeakFitLM::PeakFitLMOptions::SmallRefinementOnly );
      
      if( result.size() == inputPeak.size() )
      {
        for( size_t i = 0; i < result.size(); ++i )
        fixedPeaks.push_back( *result[i] );
        std::sort( fixedPeaks.begin(), fixedPeaks.end(), &PeakDef::lessThanByMean );
        model->setPeaks( fixedPeaks );
        return;
      }else
      {
        cerr << "refit_peaks_from_right_click was not successful" << endl;
      }//if( result.size() == inputPeak.size() ) / else
    }//if( inputPeak.size() > 1 )
    
    
    
    if( outputPeak.size() != inputPeak.size() )
    {
      WStringStream msg;
      msg << "Failed to refit peak (became insignificant), from "
      << int(inputPeak.size()) << " to " << int(outputPeak.size()) << " peaks";
      passMessage( msg.str(), "", WarningWidget::WarningMsgInfo );
      return;
    }//if( outputPeak.size() != 1 )
    
    if( inputPeak.size() > 1 )
    {
      fixedPeaks.insert( fixedPeaks.end(), outputPeak.begin(), outputPeak.end() );
      std::sort( fixedPeaks.begin(), fixedPeaks.end(), &PeakDef::lessThanByMean );
      model->setPeaks( fixedPeaks );
    }else
    {
      assert( !outputPeak.empty() );
      
      model->updatePeak( peak, outputPeak[0] );
    }//if( inputPeak.size() > 1 )
  }catch( std::exception &e )
  {
    passMessage( "Sorry, error encountered refitting ROI.", "", WarningWidget::WarningMsgInfo );
    cerr << "Error encountered refitting ROI: " << e.what() << endl;
  }
}//void refit_peaks_from_right_click(...)
*/


#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include <stdio.h>
#include <shellapi.h>

#include "SpecUtils/StringAlgo.h"

/** Get command line arguments encoded as UTF-8.
    This function just leaks the memory
 
 Note that environment variables are not in UTF-8, we could account for this
 similar to:
 wchar_t *wenvstrings = GetEnvironmentStringsW();
 ...
 */
void getUtf8Args( int &argc, char ** &argv )
{
  LPWSTR *argvw = CommandLineToArgvW( GetCommandLineW(), &argc );
  if( !argvw )
  {
    std::cout << "CommandLineToArgvW failed - good luck" << std::endl;
    return ;
  }
  
  argv = (char **)malloc(sizeof(char *)*argc);
    
  for( int i = 0; i < argc; ++i)
  {
    printf("Argument: %d: %ws\n", i, argvw[i]);
  
    const std::string asutf8 = SpecUtils::convert_from_utf16_to_utf8( argvw[i] );
    argv[i] = (char *)malloc( sizeof(char)*(asutf8.size()+1) );
    strcpy( argv[i], asutf8.c_str() );
  }//for( int i = 0; i < argc; ++i)

  // Free memory allocated for CommandLineToArgvW arguments.
  LocalFree(argvw);
}//void processCustomArgs()
#endif
