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
#include "InterSpec/SpecMeas.h"
#include "InterSpec/BatchPeak.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/BatchActivity.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"


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
        return;
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
  
  throw runtime_error( "Not implemented yet" );
  /*
  std::shared_ptr<const SpecMeas> cached_exemplar_n42;
  for( const string filename : files )
  {
    const BatchActivityFitResult fit_results
                 = fit_activities_in_file( exemplar_filename, exemplar_sample_nums,
                                     cached_exemplar_n42, filename, options );
    
    if( !cached_exemplar_n42 )
      cached_exemplar_n42 = fit_results.exemplar;
    warnings.insert(end(warnings), begin(fit_results.warnings), end(fit_results.warnings) );
    
    if( !fit_results.success )
      continue;
    
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
  }//for( const string filename : files )
    
  if( !warnings.empty() )
    cerr << endl << endl;
  for( const auto warn : warnings )
    cerr << warn << endl;
   */
}//fit_peaks_in_files(...)
  
  

BatchActivityFitResult fit_activities_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          const BatchActivityFitOptions &options )
{
  throw runtime_error( "fit_activities_in_file not implemented - see fit_peaks_in_file for start of implementation" );
  
  return {};
}//fit_peaks_in_file(...)
  
}//namespace BatchPeak

