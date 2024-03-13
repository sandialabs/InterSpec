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

#include "rapidxml/rapidxml.hpp"

#include "Minuit2/MnUserParameters.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"


using namespace std;

namespace BatchActivity
{

shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf_name )
{
  if( drf_file.empty() && drf_name.empty() )
    return nullptr;
    
  //  Check if we are over-riding the DRF to use
  if( !drf_file.empty() )
  {
    // If user specified a directory, assume its a GADRAS DRF
    if( SpecUtils::is_directory(drf_file) )
    {
      try
      {
        auto drf = make_shared<DetectorPeakResponse>();
        drf->fromGadrasDirectory( drf_file );
        return drf;
      }catch( std::exception &e )
      {
        throw runtime_error( "Directory '" + drf_file + "' wasnt a GADRAS DRF: " + string(e.what()) );
      }
    }//if( SpecUtils::is_directory(drf_file) )
    
    
    // Try a URI file
    if( SpecUtils::is_file(drf_file) && (SpecUtils::file_size(drf_file) < 64*1024) )
    {
      try
      {
        std::vector<char> data;
        SpecUtils::load_file_data( drf_file.c_str(), data );
        string url_query( begin(data), end(data) );
        SpecUtils::trim( url_query );
        auto drf = DetectorPeakResponse::parseFromAppUrl( url_query );
        if( drf )
          return drf;
      }catch( std::exception &e )
      {
        // Was not a file with URI as its contents
      }
    }else
    {
      try
      {
        string url_query = drf_file;
        SpecUtils::trim( url_query );
        auto drf = DetectorPeakResponse::parseFromAppUrl( url_query );
        if( drf )
          return drf;
      }catch( std::exception &e )
      {
        // Was not a URI
      }
    }//if( was a small file ) / else
    
    
    
    // Try a CSV/TSV file that may have multiple DRF in it
    fstream input( drf_file.c_str(), ios::binary | ios::in );
    if( !input )
      throw runtime_error( "Couldnt open DRF file '" + drf_file + "'." );
      
      
    vector<shared_ptr<DetectorPeakResponse>> drfs;
    try
    {
      vector<string> credits;
        
      // This should also cover files that could be parsed by
      //  `DetectorPeakResponse::parseSingleCsvLineRelEffDrf(string &)`
        
      DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );
    }catch(...)
    {
    }//try / catch for parseMultipleRelEffDrfCsv
      
      
    if( drfs.size() == 1 )
    {
      if( drf_name.empty() )
      {
        return drfs[0];
      }else
      {
        if( (drfs[0]->name() != drf_name) || (drfs[0]->description() != drf_name) )
          throw runtime_error( "The DRF in '" + drf_file + "' is named '"
                              + drfs[0]->name() + "', but you specified a name of '"
                                + drf_name + "'" );
        return drfs[0];
      }//if( drf_name.empty() ) / else
    }else if( drfs.size() > 1 )
    {
      string names;
      for( size_t i = 0; i < drfs.size(); ++i )
      {
        names += string(names.empty() ? "" : ", ") + "'" + drfs[i]->name() + "'";
        if( !drf_name.empty() && (drfs[i]->name() == drf_name) )
          return drfs[i];  // TODO: check if multiple DRFs with same name
      }
        
      if( drf_name.empty() )
        throw runtime_error( "No DRF name specified; possible DRFs in '" + drf_file
                              + "', are: " + names );
    }//if( drfs.size() == 1 ) / else
    
    
    // Try a Rel Eff CSV file
    input.seekg( 0, ios::beg );
      
    try
    {
      auto det = DrfSelect::parseRelEffCsvFile( drf_file );
      if( det )
        return det;
    }catch( std::exception &e )
    {
      
    }//try catch
    
    
    // Try a ECC file
    try
    {
      input.seekg( 0, ios::beg );
        
      auto answer = DetectorPeakResponse::parseEccFile( input );
      if( std::get<0>(answer) )
        return std::get<0>(answer);
    }catch( std::exception &e )
    {
      // Not a .ECC file
    }
  }else if( !drf_name.empty() )  //if( !drf_file.empty() )
  {
    try
    {
      string manufacturer = "";
      string model = drf_name;
      auto drf = DrfSelect::initARelEffDetector( SpecUtils::DetectorType::Unknown,
                                                manufacturer, model );
      if( drf )
        return drf;
    }catch( std::exception &e )
    {
      
    }//try /catch to load a standard Detector by name
    
    
    // TODO: Look for previous DRFs in the DB, that the user has used
      
  }//if( !drf_file.empty() ) / else if( !drf_name.empty() )
  
  throw runtime_error( "Could not load specified detector efficiency function."
                        + (drf_file.empty() ? string("") : " Filename='" + drf_file + "'.")
                        + (drf_name.empty() ? string("") : " Name='" + drf_name + "'.") );
  
  return nullptr;
}//shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf-name )
  
  
void fit_activities_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          const BatchActivityFitOptions &options )
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
    const BatchActivityFitResult fit_results
                 = fit_activities_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, options );
    
    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.m_exemplar_file;
    warnings.insert(end(warnings), begin(fit_results.m_warnings), end(fit_results.m_warnings) );
    
    if( fit_results.m_result_code == BatchActivity::BatchActivityFitResult::ResultCode::Success )
    {
      cout << "Success analyzing '" << filename << "'!" << endl;
      assert( fit_results.m_fit_results );
      for( const ShieldingSourceFitCalc::IsoFitStruct &nuc : fit_results.m_fit_results->fit_src_info )
      {
        cout << (nuc.fitActivity ? "\t" : "[act not fit]") << nuc.nuclide->symbol << ": "
             << PhysicalUnits::printToBestActivityUnits(nuc.activity)
             << " +- " << PhysicalUnits::printToBestActivityUnits(nuc.activityUncertainty);
        
        if( nuc.ageIsFittable )
        {
          if( nuc.fitAge )
          {
            cout << " - fit age " << PhysicalUnits::printToBestTimeUnits(nuc.age) << " +- "
                 << PhysicalUnits::printToBestActivityUnits(nuc.ageUncertainty);
          }else
          {
            cout << " [used age=" << PhysicalUnits::printToBestTimeUnits(nuc.age) << "]";
          }
        }//if( nuc.ageIsFittable )
        
        
        
        switch( nuc.sourceType )
        {
          case ShieldingSourceFitCalc::ModelSourceType::Point:
            break;
            
          case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
            break;
            
          case ShieldingSourceFitCalc::ModelSourceType::Trace:
            break;
        }//switch( nuc.sourceType )
        
        cout << endl;
      }//for( const ShieldingSourceFitCalc::IsoFitStruct &nuc : m_fit_results.fit_src_info )
      
      
      for( const ShieldingSourceFitCalc::FitShieldingInfo &shield : fit_results.m_fit_results->final_shieldings )
      {
        //const char *GammaInteractionCalc::to_str( shield.m_geometry );
        
        /** Dimensions of this shielding; the meaning of the entries differs depending on the geometry,
         or if a generic material.
         
         Spherical: ['Thickness', n/a, n/a]
         Cylinder:  ['Radius','Length',n/a]
         Rectangle: ['Width','Height','Depth']
         Generic:   ['AtomicNumber','ArealDensity',n/a]
         */
        /*
        double m_dimensions[3];
        double m_dimensionUncerts[3];
        bool m_fitDimensions[3];
        std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractionUncerts;
        std::map<const SandiaDecay::Nuclide *,double> m_traceSourceActivityUncerts;
        
        bool m_isGenericMaterial;
        bool m_forFitting;
        std::shared_ptr<const Material> m_material;
        bool m_fitMassFrac;
        std::map<const SandiaDecay::Nuclide *,double> m_nuclideFractions;
        std::vector<TraceSourceInfo> m_traceSources;
        */
      }//for( loop fit_results.m_fit_results->fit_src_info )
      
    }else
    {
      cout << "Failure: " << fit_results.m_error_msg << endl;
    }
    
    /*
    blah blah blah
    
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
          PeakModel::write_peak_csv( output_csv, leaf_name, fit_peaks, fit_results.spectrum );
          cout << "Have written '" << outcsv << "'" << endl;
        }
      }//if( SpecUtils::is_file( outcsv ) ) / else
    }//if( !options.output_dir.empty() )
    
    if( options.to_stdout )
    {
      const string leaf_name = SpecUtils::filename(filename);
      cout << "peaks for '" << leaf_name << "':" << endl;
      PeakModel::write_peak_csv( cout, leaf_name, fit_peaks, fit_results.spectrum );
      cout << endl;
    }
     */
  }//for( const string filename : files )
    
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;
}//fit_activities_in_files(...)
  
  

BatchActivityFitResult fit_activities_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          const BatchActivityFitOptions &options )
{
  //  TODO: allow specifying, not just in the exemplar N42.  Also note InterSpec defines a URL encoding for model as well
  BatchActivityFitResult result;
  result.m_options = options;
  result.m_result_code = BatchActivityFitResult::ResultCode::UnknownStatus;
  result.m_filename = filename;
  result.m_exemplar_sample_numbers = exemplar_sample_nums;
  
  
  
  if( !cached_exemplar_n42 )
  {
    auto tmp = make_shared<SpecMeas>();
    if( !tmp->load_file( exemplar_filename, SpecUtils::ParserType::Auto ) )
    {
      result.m_error_msg = "Could not load exemplar '" + exemplar_filename + "'.";
      result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenExemplar;
      return result;
    }//if( !meas->load_file( filename, ParserType::Auto ) )
    
    cached_exemplar_n42 = tmp;
  }//if( !cached_exemplar_n42 )
  
  result.m_exemplar_file = cached_exemplar_n42;
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    result.m_error_msg = "Could not initialize nuclide DecayDataBase.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  MaterialDB matdb;
  const string materialfile = SpecUtils::append_path( InterSpec::staticDataDirectory(), "MaterialDataBase.txt" );
  try
  {
    matdb.parseGadrasMaterialFile( materialfile, db, false );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not initialize shielding material database.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  const Material * const iron = matdb.material("Fe (iron)");
  if( !iron )
  {
    result.m_error_msg = "Couldnt access materials from shielding database.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntInitializeStaticResources;
    return result;
  }
  
  if( !cached_exemplar_n42 )
  {
    result.m_error_msg = "No exemplar spectrum file provided";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoExemplar;
    return result;
  }//if( !cached_exemplar_n42 )
  
  
  auto specfile = make_shared<SpecMeas>();
  if( !specfile->load_file( filename, SpecUtils::ParserType::Auto ) )
  {
    result.m_error_msg = "Could not load foreground '" + filename + "'.";
    result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenInputFile;
    return result;
  }//if( !meas->load_file( filename, ParserType::Auto ) )
  result.m_foreground_file = specfile;
  
  
  shared_ptr<SpecMeas> backfile;
  if( !options.background_subtract_file.empty() )
  {
    if( options.background_subtract_file == filename )
    {
      backfile = specfile;
      result.m_background_file = specfile;
    }else if( options.background_subtract_file == exemplar_filename )
    {
      //backfile = cached_exemplar_n42;
      result.m_background_file = cached_exemplar_n42;
    }else
    {
      auto tmp = make_shared<SpecMeas>();
      if( !tmp->load_file(options.background_subtract_file, SpecUtils::ParserType::Auto) )
      {
        result.m_error_msg = "Could not load background '" + options.background_subtract_file + "'.";
        result.m_result_code = BatchActivityFitResult::ResultCode::CouldntOpenBackgroundFile;
        return result;
      }
      
      backfile = tmp;
      result.m_background_file = backfile;
    }//if( options.background_subtract_file == filename ) / else
  }//if( !options.background_subtract_file.empty() )
  
  
  // Find the sample number to use for either the foreground, or background measurement
  auto find_sample = []( shared_ptr<const SpecMeas> meas, const SpecUtils::SourceType wanted ) -> int {
    // This logic is probably repeated in a number of places throughout InterSpec now...
    assert( (wanted == SpecUtils::SourceType::Foreground)
           || (wanted == SpecUtils::SourceType::Background) );
    
    if( (wanted != SpecUtils::SourceType::Foreground)
           && (wanted != SpecUtils::SourceType::Background) )
    {
      throw std::logic_error( "Invalid src type" );
    }
    
    if( !meas )
      throw runtime_error( "No SpecMeas passed in." );
    
    const vector<string> &detectors = meas->detector_names();
    const set<int> &sample_nums = meas->sample_numbers();
    
    set<int> foreground_samples, background_samples;
    
    for( const int sample : sample_nums )
    {
      bool classified_sample_num = false;
      for( const string &det : detectors )
      {
        auto m = meas->measurement( sample, det );
        if( !m )
          continue;
        
        switch( m->source_type() )
        {
          case SpecUtils::SourceType::IntrinsicActivity:
          case SpecUtils::SourceType::Calibration:
            break;
            
          case SpecUtils::SourceType::Background:
            classified_sample_num = true;
            background_samples.insert(sample);
            break;
            
          case SpecUtils::SourceType::Foreground:
          case SpecUtils::SourceType::Unknown:
            classified_sample_num = true;
            foreground_samples.insert(sample);
            break;
        }//switch( m->source_type() )
        
        if( classified_sample_num )
          break;
      }//for( const string &det : detectors )
    }//for( const int sample : sample_nums )
    
    const set<int> &samples = (wanted == SpecUtils::SourceType::Foreground) ? foreground_samples 
                                                                            : background_samples;
    if( samples.size() != 1 )
      throw runtime_error( "Sample number to use could not be uniquely identified." );
    
    return *begin(samples);
  };//int find_sample(...)
  
  
  shared_ptr<const SpecUtils::Measurement> foreground;
  try
  {
    const int sample_num = find_sample( specfile, SpecUtils::SourceType::Foreground );
    
    try
    {
      vector<shared_ptr<const SpecUtils::Measurement>> meass = specfile->sample_measurements(sample_num);
      if( meass.size() == 1 )
        foreground = meass[0];
      else
        foreground = specfile->sum_measurements( {sample_num}, specfile->detector_names(), nullptr );
      if( !foreground )
        throw runtime_error( "Missing measurement." );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Could getting foreground measurement: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundSampleNumberUnderSpecified;
      return result;
    }//Try / catch get
    
    result.m_foreground = foreground;
    result.m_foreground_sample_numbers = set<int>{ sample_num };
    assert( result.m_foreground_file );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not determine foreground: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundSampleNumberUnderSpecified;
    return result;
  }// try / catch (find foreground)

  
  set<int> background_sample_nums;
  shared_ptr<const SpecUtils::Measurement> background;
  try
  {
    if( options.background_subtract_file == exemplar_filename )
    {
      if( options.background_subtract_samples.empty() )
      {
        const int sample_num = find_sample( cached_exemplar_n42, SpecUtils::SourceType::Background );
        background_sample_nums = set<int>{sample_num};
      }else
      {
        background_sample_nums = options.background_subtract_samples;
      }
      
      try
      {
        background = specfile->sum_measurements( background_sample_nums, cached_exemplar_n42->detector_names(), nullptr );
      }catch( std::exception &e )
      {
        result.m_error_msg = "Error summing background measurements from exemplar: " + string(e.what());
        result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
        return result;
      }//Try / catch get
    }else if( backfile )
    {
      if( !options.background_subtract_samples.empty() )
      {
        const SpecUtils::SourceType t = (backfile == specfile)
        ? SpecUtils::SourceType::Background
        : SpecUtils::SourceType::Foreground;
        
        const int sample_num = find_sample( backfile, t );
        
        background_sample_nums = set<int>{sample_num};
        
        try
        {
          vector<shared_ptr<const SpecUtils::Measurement>> meass = backfile->sample_measurements(sample_num);
          if( meass.size() == 1 )
            background = meass[0];
          else
            background = backfile->sum_measurements( background_sample_nums, backfile->detector_names(), nullptr );
          if( !background )
            throw runtime_error( "Missing measurement." );
        }catch( std::exception &e )
        {
          result.m_error_msg = "Could not get background measurement: " + string(e.what());
          result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
          return result;
        }//Try / catch get
      }else
      {
        background_sample_nums = options.background_subtract_samples;
        try
        {
          background = backfile->sum_measurements( background_sample_nums, backfile->detector_names(), nullptr );
          if( !background )
            throw runtime_error( "Missing measurement." );
        }catch( std::exception &e )
        {
          result.m_error_msg = "Could not sum background samples: " + string(e.what());
          result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
          return result;
        }//Try / catch get
      }//if( options.background_subtract_samples.empty() ) / else
    }//if( options.background_subtract_file == exemplar_filename ) / else
    
    result.m_background = background;
    result.m_background_sample_numbers = background_sample_nums;
    assert( (!!result.m_background_file) == (!!result.m_background) );
    assert( options.background_subtract_file.empty() == (!background) );
    assert( options.background_subtract_file.empty() == (!result.m_background_file) );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Could not determine background: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::BackgroundSampleNumberUnderSpecified;
    return result;
  }// try / catch (find foreground)
  
  
  // Now need to figure out the foreground and possibly background sample numbers in the exemplar
  BatchPeak::BatchPeakFitOptions peak_fit_options = options;
  peak_fit_options.background_subtract_file = "";
  

  const BatchPeak::BatchPeakFitResult foreground_peak_fit_result
              = BatchPeak::fit_peaks_in_file( exemplar_filename, exemplar_sample_nums,
                                            cached_exemplar_n42, filename, specfile,
                                            result.m_foreground_sample_numbers, peak_fit_options );
  
  result.m_peak_fit_results = make_shared<BatchPeak::BatchPeakFitResult>( foreground_peak_fit_result );
  
  if( !foreground_peak_fit_result.success )
  {
    result.m_error_msg = "Fitting of foreground peaks failed.";
    result.m_result_code = BatchActivityFitResult::ResultCode::ForegroundPeakFitFailed;
  }//if( !foreground_peaks.success )
  
  
  if( !options.background_subtract_file.empty() )
  {
    //if( options.background_subtract_file == exemplar_filename )
    //First see if peaks are already fit - if so, use them, otherwise try to fit them
    
    // TODO: we are fitting the foreground peaks to the background - we should allow specifying exemplar background peaks - either to refit, or to just use
    cerr << "TODO: we are fitting the foreground peaks to the background - we should allow specifying exemplar background peaks - either to refit, or to just use" << endl;
    const BatchPeak::BatchPeakFitResult background_peaks
                = BatchPeak::fit_peaks_in_file( exemplar_filename, 
                                               exemplar_sample_nums,
                                              cached_exemplar_n42, 
                                               options.background_subtract_file,
                                               backfile, 
                                               result.m_background_sample_numbers,
                                               peak_fit_options );
    
    
    
    result.m_background_peak_fit_results = make_shared<BatchPeak::BatchPeakFitResult>( background_peaks );
  }//if( backfile )
  
  
  const rapidxml::xml_document<char> *model_xml = cached_exemplar_n42->shieldingSourceModel();
  const rapidxml::xml_node<char> *base_node = model_xml ? model_xml->first_node() : nullptr;
  if( !base_node )
  {
    result.m_error_msg = "Exemplar file did not have a Shielding/Source fit model in it";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoInputSrcShieldModel;
    return result;
  }
  
  
  deque<shared_ptr<const PeakDef>> foreground_peaks = foreground_peak_fit_result.fit_peaks;
  if( foreground_peaks.empty() )
  {
    result.m_error_msg = "No foreground peaks fit.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoFitForegroundPeaks;
    return result;
  }
  
  shared_ptr<const DetectorPeakResponse> detector = options.drf_override;
  if( !detector )
    detector = cached_exemplar_n42->detector();
  
  if( !detector )
  {
    result.m_error_msg = "No detector efficiency function specified.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoDetEffFnct;
    return result;
  }
  
  
  double distance = 1*PhysicalUnits::meter;
  if( !detector->isFixedGeometry() )
  {
    try
    {
      const rapidxml::xml_node<char> *dist_node = base_node->first_node( "Distance" );
      const string diststr = SpecUtils::xml_value_str( dist_node );
      distance = PhysicalUnits::stringToDistance( diststr );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Failed to get distance.";
      result.m_result_code = BatchActivityFitResult::ResultCode::InvalidDistance;
      return result;
    }//try / catch to get distance
  }//if( !detector->isFixedGeometry() )
  
  GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  
  const rapidxml::xml_node<char> *geom_node = base_node->first_node( "Geometry" );
  const string geomstr = SpecUtils::xml_value_str( geom_node );
  
  for( GammaInteractionCalc::GeometryType type = GammaInteractionCalc::GeometryType(0);
      type != GammaInteractionCalc::GeometryType::NumGeometryType;
      type = GammaInteractionCalc::GeometryType(static_cast<int>(type) + 1) )
  {
    if( SpecUtils::iequals_ascii(geomstr, GammaInteractionCalc::to_str(type)) )
    {
      geometry = type;
      break;
    }
  }//for( loop over geometry types )
  
  if( geometry == GammaInteractionCalc::GeometryType::NumGeometryType )
  {
    result.m_error_msg = "Geometry type not specified.";
    result.m_result_code = BatchActivityFitResult::ResultCode::InvalidGeometry;
    return result;
  }//if( geometry == GammaInteractionCalc::GeometryType::NumGeometryType )
  
  
  vector<ShieldingSourceFitCalc::ShieldingInfo> shield_definitions;
  vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  ShieldingSourceFitCalc::ShieldingSourceFitOptions fit_options;
  try
  {
    fit_options.deSerialize( base_node );
  }catch( std::exception &e )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' Options invalid: " + string(e.what());
    result.m_result_code = BatchActivityFitResult::ResultCode::InvalidFitOptions;
    return result;
  }

  // Check that if options said to use background-subtraction, that we actually can.
  //  I'm a little mixed here - perhaps if no background is specified, we should just ignore this,
  //  or just create a warning???
  //  For now - lets just warn.
  shared_ptr<const deque<std::shared_ptr<const PeakDef>>> background_peaks;
  if( background && result.m_background_peak_fit_results && result.m_background_peak_fit_results->success )
  {
    background_peaks = make_shared<deque<shared_ptr<const PeakDef>>>(result.m_background_peak_fit_results->fit_peaks);
  }
  
  
  
  if( fit_options.background_peak_subtract && (!background_peaks || background_peaks->empty()) )
  {
    fit_options.background_peak_subtract = false;
    string msg = "Exemplar file '" + exemplar_filename + "' Options indicated"
    " background-subtract, but either no background was given, or no background peaks fit.";
    result.m_warnings.push_back( msg );
    //result.m_error_msg = msg;
    //result.m_result_code = BatchActivityFitResult::ResultCode::ExemplarUsedBackSubButNoBackground;
    //return result;
  }
                                           
  const rapidxml::xml_node<char> *shieldings_node = base_node->first_node( "Shieldings" );
  if( !shieldings_node )
  {
    result.m_error_msg = "No Shieldings node specified in Shield/Src fit model.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoShieldingsNode;
    return result;
  }
      
  XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
  {
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo info;
      info.deSerialize( shield_node, &matdb );
      shield_definitions.push_back( std::move(info) );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Invalid Shielding node: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::ErrorParsingShielding;
      return result;
    }
  }//XML_FOREACH_CHILD(shield_node, shieldings_node, "Shielding")
      
      
  const rapidxml::xml_node<char> *srcs_node = base_node->first_node( "Nuclides" );
  if( !srcs_node )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' missing Nuclides (sources) node.";
    result.m_result_code = BatchActivityFitResult::ResultCode::MissingNuclidesNode;
    return result;
  }//if( !srcs_node )
      
  XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
  {
    try
    {
      ShieldingSourceFitCalc::SourceFitDef info;
      info.deSerialize( src_node );
      src_definitions.push_back( std::move(info) );
    }catch( std::exception &e )
    {
      result.m_error_msg = "Exemplar file '" + exemplar_filename + "' has invalid Nuclides node: " + string(e.what());
      result.m_result_code = BatchActivityFitResult::ResultCode::InvalidNuclideNode;
      return result;
    }// try / catch
  }//XML_FOREACH_CHILD( src_node, srcs_node, "Nuclide" )
      
  
  if( src_definitions.empty() )
  {
    result.m_error_msg = "Exemplar file '" + exemplar_filename + "' no sources defined.";
    result.m_result_code = BatchActivityFitResult::ResultCode::NoSourceNuclides;
    return result;
  }//if( src_definitions.empty() )
  
  // TODO: We may not have fit peaks for all `src_definitions`.  We should remove these sources, and add warnings about it

  // We have all the parts, lets do the computation:
  pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
  GammaInteractionCalc::ShieldingSourceChi2Fcn::create( distance, geometry,
                                                           shield_definitions, src_definitions, detector,
                                                           foreground, background, foreground_peaks, background_peaks, fit_options );
      
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  *inputPrams = fcn_pars.second;
      
  auto progress = make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto fit_results = make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  result.m_fit_results = fit_results;
  
  auto progress_fcn = [progress](){
    // We dont really care too much probably right now
  };
      
  bool finished_fit_called = false;
  auto finished_fcn = [fit_results, &finished_fit_called](){
    finished_fit_called = true;
  };
      
  ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, fit_results, finished_fcn );
      
  assert( finished_fit_called );
  
  result.m_result_code = BatchActivityFitResult::ResultCode::Success;
  
  if( fit_results->successful != ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final )
  {
    result.m_error_msg = "Fit not successful.";
    result.m_result_code = BatchActivityFitResult::ResultCode::FitNotSuccessful;
    return result;
  }
  
      
  for( const ShieldingSourceFitCalc::SourceFitDef insrc : src_definitions )
  {
    assert( insrc.nuclide );
    bool found_in_output = false;
    for( const ShieldingSourceFitCalc::IsoFitStruct &fitsrc : fit_results->fit_src_info )
    {
      assert( fitsrc.nuclide );
      found_in_output = (fitsrc.nuclide == insrc.nuclide);
      if( found_in_output )
        break;
    }//for( const ShieldingSourceFitCalc::IsoFitStruct &fitsrc : fit_results->fit_src_info )
    
    if( !found_in_output )
    {
      result.m_warnings.push_back( "Fit results does not include " + insrc.nuclide->symbol );
      result.m_result_code = BatchActivityFitResult::ResultCode::DidNotFitAllSources;
    }//if( !found_in_output )
  }//for( const ShieldingSourceFitCalc::SourceFitDef insrc : src_definitions )
  
  cout << "Finished fitting for '" << filename << "'" << endl;
  
  return result;
}//fit_peaks_in_file(...)
  
}//namespace BatchPeak

