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

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>


#include "QR-Code-generator/cpp/qrcodegen.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/UriSpectrum.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/QrCode.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/QRSpectrum.h"


// Headers only needed for `int dev_code();`
#include <mutex>
#include "SpecUtils/SpecUtilsAsync.h"


using namespace std;
using SpecUtils::EncodeOptions;
using SpecUtils::UrlSpectrum;
using SpecUtils::to_url_spectra;

namespace
{
  const int sm_qr_quite_area = 3;
  

  void make_example_qr_codes()
  {
   // This function recursively looks for N42 files in a given directory, and then encodes
   //  a randomly selected 100 of them into QR codes, both as individual SVG files, as well
   //  as them all into a single HTML file, that has a hyperlink and encodeing information as well.
   //  The foreground is encoded alone into QR, then background, then if possible
   //  foreground+background.
   //
   //  Assumes "`base_dir`/svgs/[Low|Medium|Quartile|High]/" directory exists.
   //
   //  Currently set up to pull the foreground and background "Derived" spectra from Verifinder
   //  files and write those out.
   //
    using namespace SpecUtils;
    using namespace QRSpectrum;
    
    const char *base_dir = "...";
    const size_t wanted_num_files = 100;
    const auto wanted_ecl = QrErrorCorrection::Low;
    
    auto ecl_str = []( QrErrorCorrection ecl ) -> string {
      switch ( ecl )
      {
        case QrErrorCorrection::Low: return "Low";
        case QrErrorCorrection::Medium: return "Medium";
        case QrErrorCorrection::Quartile: return "Quartile";
        case QrErrorCorrection::High: return "High";
      }
      return "invalid";
    };
    
    const vector<string> files = SpecUtils::recursive_ls(base_dir, ".n42");
    
    set<int> files_tried;
    size_t num_succes_files = 0;
    
    ofstream output_html( SpecUtils::append_path(base_dir, "qr_codes.html") );
    output_html <<
    R"delim(<!DOCTYPE html>
    <html lang="en">
      <head>
        <meta charset="utf-8">
        <title>Example QR Codes for Verifinder SN20</title>
      </head>
      <body>
    )delim";
    
    output_html << "Target error correction level: " << ecl_str(wanted_ecl) << "<br />";
    output_html << "Quite area, num elements: " << sm_qr_quite_area << "<br />";
    
    int index = 0;
    while( (files_tried.size() != files.size()) && (num_succes_files < wanted_num_files) )
    {
      if( wanted_num_files > files.size() )
      {
        index += 1;
      }else
      {
        index = rand() % files.size();
        while( files_tried.count(index) )
          index = rand() % files.size();
      }
      
      files_tried.insert( index );
      const string filename = files[index];
      
      SpecUtils::SpecFile spec;
      if( !spec.load_file( filename, SpecUtils::ParserType::Auto, filename) )
        continue;
      
      shared_ptr<const SpecUtils::Measurement> foreground, background;
      for( const auto &m : spec.measurements() )
      {
        if( !m || !m->derived_data_properties() || m->num_gamma_channels() < 32 )
          continue;
        
        const string title = m->title();
        if( SpecUtils::icontains(title, "Analysis0Background") )
          background = m;
        else if( SpecUtils::icontains(title, "Analysis0Raw") )
          foreground = m;
      }//for( const auto &m : spec->measurements() )
      
      if( !foreground || !background )
      {
        cout << "Failed to get foreground or background\n";
        continue;
      }
      
      num_succes_files += 1;
      
      output_html << "<br /><br /><h1>" << SpecUtils::filename(filename) << "</h1>\n";
      
      vector<SpecUtils::UrlSpectrum> fore_url_spec = SpecUtils::to_url_spectra( {foreground}, spec.instrument_model() );
      vector<SpecUtils::UrlSpectrum> back_url_spec = SpecUtils::to_url_spectra( {background}, spec.instrument_model() );
      vector<SpecUtils::UrlSpectrum> both_url_spec = SpecUtils::to_url_spectra( {foreground, background}, spec.instrument_model() );
      
      auto doEncode = [wanted_ecl,filename,ecl_str]( const vector<SpecUtils::UrlSpectrum> &url_spec ) -> vector<QrCodeEncodedSpec> {

        return qr_code_encode_spectra( url_spec, wanted_ecl, 0 );
        
        /*
#warning "Doing sanity check on what encoded best"
        static int ntimes = 0;
        if( ntimes++ < 4 )
          cout << "Doing sanity check on what encoded best" << endl;
        
        const uint8_t encode_options[] {
          0x0,
          //static_cast<uint8_t>(EncodeOptions::NoBaseXEncoding),
          //static_cast<uint8_t>(EncodeOptions::CsvChannelData),
          //static_cast<uint8_t>(EncodeOptions::CsvChannelData | EncodeOptions::NoBaseXEncoding),
          //static_cast<uint8_t>(EncodeOptions::NoZeroCompressCounts),
          //static_cast<uint8_t>(EncodeOptions::NoZeroCompressCounts | EncodeOptions::NoBaseXEncoding),

          //static_cast<uint8_t>(EncodeOptions::CsvChannelData | EncodeOptions::NoBaseXEncoding),
          //static_cast<uint8_t>(EncodeOptions::NoZeroCompressCounts),
          //static_cast<uint8_t>(EncodeOptions::CsvChannelData | EncodeOptions::UseUrlSafeBase64),
          //static_cast<uint8_t>(EncodeOptions::NoZeroCompressCounts | EncodeOptions::UseUrlSafeBase64),
          static_cast<uint8_t>(EncodeOptions::UseUrlSafeBase64 | EncodeOptions::AsMailToUri),
        };
        
        const auto option_to_str = []( uint8_t options ) -> string {
          string answer;
          if( options & EncodeOptions::NoDeflate )
            answer += (answer.empty() ? "" : " | ") + string("NoDeflate");
          if( options & EncodeOptions::NoBaseXEncoding )
            answer += (answer.empty() ? "" : " | ") + string("NoBaseXEncoding");
          if( options & EncodeOptions::CsvChannelData )
            answer += (answer.empty() ? "" : " | ") + string("CsvChannelData");
          if( options & EncodeOptions::NoZeroCompressCounts )
            answer += (answer.empty() ? "" : " | ") + string("NoZeroCompressCounts");
          if( options & EncodeOptions::UseUrlSafeBase64 )
            answer += (answer.empty() ? "" : " | ") + string("UseUrlSafeBase64");
          if( options & EncodeOptions::AsMailToUri )
            answer += (answer.empty() ? "" : " | ") + string("AsMailToUri");
          return answer;
        };//option_to_str lamda
        
        
        uint8_t used_option = 0x0;
        int min_qr_size = 10000000;
        size_t min_uri_size = 100000000;
        vector<QrCodeEncodedSpec> encoded;
        for( uint8_t option : encode_options )
        {
          try
          {
            vector<QrCodeEncodedSpec> this_try = url_encode_spectra( url_spec, wanted_ecl, option );
            if( this_try.size() == 1 )
            {
              if( this_try[0].m_qr_size < min_qr_size )
              {
                encoded = this_try;
                used_option = option;
                min_qr_size = this_try[0].m_qr_size;
                min_uri_size = this_try[0].m_url.size();
              }else if( (this_try[0].m_qr_size == min_qr_size)
                       && (this_try[0].m_url.size() < min_uri_size) )
              {
                encoded = this_try;
                used_option = option;
                min_qr_size = this_try[0].m_qr_size;
                min_uri_size = this_try[0].m_url.size();
              }else
              {
                //cout << "  It was " << this_try[0].m_qr_size << " vs best of " << min_qr_size
                //<< " (URI size " << this_try[0].m_url.size() << " vs " << min_uri_size << ")" << endl;
              }
            }else if( encoded.empty() )
            {
              encoded = this_try;
              used_option = option;
            }
          }catch( std::exception &e )
          {
            cerr << "Failed: " << e.what() << endl;
          }
        }//for( uint8_t option : encode_options )
        
        if( encoded.size() == 1 )
        {
          cout << "For " << SpecUtils::filename(filename) << " encode_options="
          << static_cast<int>(used_option)
          << " (" << option_to_str(used_option) << ")"
          << " with size QR " << encoded[0].m_qr_size  //m_qr_size in range [21, 177]
          << " and EC=" << ecl_str(encoded[0].m_error_level)
          << ", len(url)=" << encoded[0].m_url.size()
          //<< ",\nURL=\"" << encoded[0].m_url << "\"" << endl << endl
          << endl;
        }else
        {
          cout << "For " << SpecUtils::filename(filename)
          << " (" << option_to_str(used_option) << ") - failed." << endl;
        }
        
        return encoded;
         */
      };//doEncode
      
      auto toHtml = [&output_html,ecl_str,filename,base_dir,wanted_ecl]( const QrCodeEncodedSpec &urlspec, string tag ){
        output_html << "<div style=\"width: 350px; height: 350px;\">\n" << urlspec.m_qr_svg << "\n</div>\n";
        output_html << "Num elements: " << urlspec.m_qr_size << "x" << urlspec.m_qr_size << " (version " << urlspec.m_qr_version << ").<br />\n";
        output_html << "Error Correction Level: " << ecl_str(urlspec.m_error_level) << "<br />\n";
        output_html << "Link: <a href=\"" <<  urlspec.m_url << "\">here</a>" << "<br />\n";
        
        const string leaf = SpecUtils::filename(filename);
        string svg_name = SpecUtils::append_path(
                          SpecUtils::append_path(
                          SpecUtils::append_path( base_dir, "svgs"), ecl_str(wanted_ecl) ), leaf + "." + tag + ".svg" );
        ofstream svg( svg_name.c_str(), ios::out | ios:: binary );
        svg << urlspec.m_qr_svg << endl;
      };
    
      /*
      try
      {
        output_html << "<h3>Foreground</h3>\n";
        //const vector<QrCodeEncodedSpec> encoded = url_encode_spectra( fore_url_spec, wanted_ecl, 0 );
        const vector<QrCodeEncodedSpec> encoded = doEncode( fore_url_spec );
        
        
        if( encoded.size() != 1 )
        {
          cout << "Foreground larger than 1" << endl;
          output_html << "Encoded to " << encoded.size() << " URLS - skipping\n<br />";
        }else
        {
          toHtml( encoded[0], "fore" );
        }
      }catch( std::exception &e )
      {
        cerr << "Failed to encode foreground to UrlSpectrum" << endl;
        output_html << "Failed to encode foreground to UrlSpectrum: " << e.what() << "\n<br />";
      }
      
      try
      {
        output_html << "<h3>Background</h3>\n";
        //const vector<QrCodeEncodedSpec> encoded = url_encode_spectra( back_url_spec, wanted_ecl, 0 );
        const vector<QrCodeEncodedSpec> encoded = doEncode( back_url_spec );
        
        if( encoded.size() != 1 )
        {
          cout << "Background larger than 1" << endl;
          output_html << "Encoded to " << encoded.size() << " URLS - skipping\n<br />";
        }else
        {
          toHtml( encoded[0], "back" );
        }
      }catch( std::exception &e )
      {
        cerr << "Failed to encode background to UrlSpectrum" << endl;
        output_html << "Failed to encode background to UrlSpectrum: " << e.what() << "\n<br />";
      }
      */
      
      try
      {
        output_html << "<h3>Foreground + Background</h3>\n";
        //const vector<QrCodeEncodedSpec> encoded = qr_code_encode_spectra( both_url_spec, wanted_ecl, 0 );
        const vector<QrCodeEncodedSpec> encoded = doEncode( both_url_spec );
        
        if( encoded.size() != 1 )
        {
          cout << "For " << SpecUtils::filename(filename) << ", Foreground + Background larger than 1" << endl;
          output_html << "Encoded to " << encoded.size() << " URLS - skipping\n<br />";
        }else
        {
          toHtml( encoded[0], "fore_back" );
        }
      }catch( std::exception &e )
      {
        cerr << "Failed to encode foreground+background to UrlSpectrum" << endl;
        output_html << "Failed to encode foreground+background to UrlSpectrum: " << e.what() << "\n<br />";
      }
      
      //cout << "From '" << filename << "', got\n\tNumForeground: " << foregrounds.size() << ", \tNumBack: " << backgrounds.size() << endl;
    }//for( const string filename : files )
    
    output_html <<
    R"delim(
    </body>
    </html>
    )delim";
    
    cout << "Done" << endl;
  }//void make_example_qr_codes()
}//namespace


namespace QRSpectrum
{

std::vector<QrCodeEncodedSpec> qr_code_encode_spectra( const std::vector<SpecUtils::UrlSpectrum> &measurements,
                                           const QrErrorCorrection minErrorCorrection,
                                           const uint8_t encode_options
                                           )
{
// It looks like qrcodegen::QrCode::encodeText(...) will default to binary if the input isnt text
//  const bool alpha_num_qr_encode = ((!use_url_safe_base64) && (use_baseX_encoding || (!use_deflate && !use_bin_chan_data)));
  
  if( measurements.empty() )
    throw runtime_error( "qr_code_encode_spectra: no input spectra." );
  
  qrcodegen::QrCode::Ecc ecc = qrcodegen::QrCode::Ecc::LOW;
  switch( minErrorCorrection )
  {
    case QrErrorCorrection::Low:      ecc = qrcodegen::QrCode::Ecc::LOW;      break;
    case QrErrorCorrection::Medium:   ecc = qrcodegen::QrCode::Ecc::MEDIUM;   break;
    case QrErrorCorrection::Quartile: ecc = qrcodegen::QrCode::Ecc::QUARTILE; break;
    case QrErrorCorrection::High:     ecc = qrcodegen::QrCode::Ecc::HIGH;     break;
  }//switch( minErrorCorrection )
  
  vector<string> urls;
  vector<qrcodegen::QrCode> qrs;

  
  //Now need to check that it will encode into a QR
  bool success_encoding = false;
  const size_t max_parts = (measurements.size() == 1) ? 9 : 1;
  for( size_t num_parts = 1; !success_encoding && (num_parts <= max_parts); ++num_parts )
  {
    urls = SpecUtils::url_encode_spectra( measurements, encode_options, num_parts );
    assert( urls.size() == num_parts );
    if( urls.size() != num_parts )
      throw std::logic_error( "Unexpected number of URLs" );
    
    qrs.clear();
    
    bool all_urls_small_enough = true;
    
    for( size_t url_num = 0; url_num < urls.size(); ++url_num )
    {
      try
      {
        // It looks like qrcodegen::QrCode::encodeText(...) will default to binary if the input isnt text
        //if( alpha_num_qr_encode )
        //{
        qrs.push_back( qrcodegen::QrCode::encodeText( urls[url_num].c_str(), ecc ) );
        //}else
        //{
        // qrcodegen::QrCode::encodeBinary: The maximum number of bytes allowed is 2953
        //  vector<uint8_t> data( url.size() );
        //  memcpy( &(data[0]), &(url[0]), url.size() );
        //  qrs.push_back( qrcodegen::QrCode::encodeBinary( data, ecc ) );
        //}
      }catch( ... )
      {
        all_urls_small_enough = false;
        break;
      }
    }//for( size_t url_num = 0; url_num < urls.size(); ++url_num )
    
    success_encoding = all_urls_small_enough;
  }//for( size_t num_parts = 0; num_parts < 10; ++num_parts )
  
  if( !success_encoding )
  {
    vector<string> urlv;
    try
    {
      urlv = url_encode_spectra( measurements, encode_options, 1 );
    }catch( std::exception & )
    {
      assert( 0 );
    }
    
    throw runtime_error( "qr_code_encode_spectra: Failed to encode "
                        + std::to_string(measurements.size())+ " spectra into a max of "
                        + std::to_string(max_parts) + "QR codes at the desired error"
                        " correction level (len(url)="
                        + std::to_string( urlv.empty() ? size_t(0) : urlv[0].size() ) + ")" );
  }//if( !success_encoding )
  
  // See how much data we are actually getting into a URL
  //  (we do get about the max of 4296 bytes expected)
  //if( urls.size() == 1 )
  //{
  //  static size_t max_size_encoded = 0;
  //  if( urls[0].size() > max_size_encoded )
  //  {
  //    max_size_encoded = urls[0].size();
  //    cout << "\nEncoded " << max_size_encoded << " bytes into URL.\n" << endl;
  //  }
  //}
  
  
  assert( (urls.size() == qrs.size()) && (qrs.size() >= 1) );
  
  std::vector<QrCodeEncodedSpec> answer;
  for( size_t i = 0; (i < urls.size()) && (i < qrs.size()); ++i )
  {
    QrCodeEncodedSpec spec;
    spec.m_url = urls[i];
    spec.m_qr_svg = QrCode::to_svg_string( qrs[i], sm_qr_quite_area );
    
    spec.m_qr_size = qrs[i].getSize();
    spec.m_qr_version = qrs[i].getVersion();
    
    switch( qrs[i].getErrorCorrectionLevel() )
    {
      case qrcodegen::QrCode::Ecc::LOW:      spec.m_error_level = QrErrorCorrection::Low; break;
      case qrcodegen::QrCode::Ecc::MEDIUM:   spec.m_error_level = QrErrorCorrection::Medium; break;
      case qrcodegen::QrCode::Ecc::QUARTILE: spec.m_error_level = QrErrorCorrection::Quartile; break;
      case qrcodegen::QrCode::Ecc::HIGH:     spec.m_error_level = QrErrorCorrection::High; break;
    }//switch( qrs[1].getErrorCorrectionLevel() )
    
    answer.push_back( spec );
  }//for( size_t i = 0; (i < urls.size()) && (i < qrs.size()); ++i )
  
  return answer;
}//qr_code_encode_spectra(...)




int dev_code()
{
  make_example_qr_codes();
  return 0;
  
  // Now that URL encoding has solidified - this function should be cleaned up to use that code
  //  to compute statistics of how many QR codes it typically takes and such.
#define DELETE_UNWANTED_FILES 0

  const char *base_dir = "/Users/wcjohns/rad_ana/qrspec_test_data";
  
  const vector<string> files = SpecUtils::recursive_ls(base_dir);
  
  std::mutex result_mutex;
  map<pair<size_t,string>,SpecUtils::UrlSpectrum> prev_spec;
  
  map<pair<size_t,string>, vector<size_t>>
  data_sizes_bin_zlib_b45, data_sizes_bin_nozero, data_sizes_bin, data_sizes_csv, data_sizes_bin_zlib, data_sizes_csv_zlib_b45,
  num_qr_code_single_spec, num_qr_code_two_spec, qr_code_level_single_spec, qr_code_version_single_spec, qr_code_error_level_single_spec;
  
  
  SpecUtilsAsync::ThreadPool pool;
  
  // Processing to many spectra at a time is super inefficient, so we'll only do just enough to fill up the CPU.
  const int batch_size = 6;
  int thread_num = 0;
  
  for( string filename : files )
  {
    auto do_work = [&]( string filename ){
      
      if( SpecUtils::likely_not_spec_file(filename) )
      {
#if( DELETE_UNWANTED_FILES )
        SpecUtils::remove_file( filename );
#endif
        return;
      }
      
      SpecUtils::SpecFile spec;
      if( !spec.load_file( filename, SpecUtils::ParserType::Auto, filename) )
      {
#if( DELETE_UNWANTED_FILES )
        SpecUtils::remove_file( filename );
#endif
        return;
      }
      
      vector<shared_ptr<const SpecUtils::Measurement>> foreground, background, usable_spectra;
      for( shared_ptr<const SpecUtils::Measurement> m : spec.measurements() )
      {
        if( (m->gamma_count_sum() < 100) || (m->num_gamma_channels() < 255)
           || ((m->live_time() <= 0.0) && (m->real_time() <= 0.0)) )
          continue;
        
        shared_ptr<const SpecUtils::EnergyCalibration> cal = m->energy_calibration();
        if( !cal || !cal->valid() || (cal->type() != SpecUtils::EnergyCalType::Polynomial) )
          continue;
        
        switch( m->source_type() )
        {
          case SpecUtils::SourceType::IntrinsicActivity:
          case SpecUtils::SourceType::Calibration:
            break;
            
          case SpecUtils::SourceType::Background:
            background.push_back( m );
            usable_spectra.push_back( m );
            break;
            
          case SpecUtils::SourceType::Foreground:
            foreground.push_back( m );
            usable_spectra.push_back( m );
            break;
            
          case SpecUtils::SourceType::Unknown:
            if( spec.detector_type() != SpecUtils::DetectorType::KromekD3S )
            {
              foreground.push_back( m );
              usable_spectra.push_back( m );
            }
            break;
        }//switch( m->source_type() )
      }//for( shared_ptr<const SpecUtils::Measurement> m : spec.measurements() )
      
      if( (foreground.size() > 1) || (background.size() > 1) || (usable_spectra.size() < 1) )
      {
#if( DELETE_UNWANTED_FILES )
        SpecUtils::remove_file( filename );
#endif
        return;
      }
      
      string model;
      if( spec.detector_type() != SpecUtils::DetectorType::Unknown )
        model = detectorTypeToString( spec.detector_type() );
      if( model.empty() )
        model = spec.instrument_model();
      
      // Remove equal signs and quotes
      SpecUtils::ireplace_all( model, ":", "" );
      SpecUtils::ireplace_all( model, "\"", "" );
      SpecUtils::ireplace_all( model, "'", "" );
      
      
      const auto key = pair<size_t,string>( spec.num_gamma_channels(), model );
      
      if( (foreground.size() == 1) && (background.size() == 1) )
      {
        try
        {
          vector<SpecUtils::UrlSpectrum> urlspec = SpecUtils::to_url_spectra( {foreground[0], background[0]}, model );
          vector<QrCodeEncodedSpec> encspecs = qr_code_encode_spectra( urlspec, QrErrorCorrection::Low, 0 );
          
          {
            std::lock_guard<std::mutex> lock( result_mutex );
            num_qr_code_two_spec[key].push_back( encspecs.size() );
          }
          
          assert( encspecs.size() == 1 );
          //cout << "Fit two spectra in URL:\n\t" << encspecs[0].m_url << endl << endl << encspecs[0].m_qr_svg << endl << endl << endl;
          
          try
          {
            vector<SpecUtils::UrlSpectrum> decoded = SpecUtils::decode_spectrum_urls( { SpecUtils::url_decode(encspecs[0].m_url) } );
            assert( decoded.size() == 2 );
          }catch( std:: exception &e )
          {
            cerr << "Error decoding multiple URLs: " << e.what() << endl;
            cerr << endl;
          }
          
        }catch( std::exception &e )
        {
          std::lock_guard<std::mutex> lock( result_mutex );
          num_qr_code_two_spec[key].push_back( 0 );
        }
      }else if( usable_spectra.size() == 1 )
      {
        vector<SpecUtils::UrlSpectrum> urlspec = SpecUtils::to_url_spectra( {usable_spectra[0]}, model );
        assert( urlspec.size() == 1 );
        
        bool had_prev = false;
        SpecUtils::UrlSpectrum prev;
        
        {// begin lock on result_mutex
          std::lock_guard<std::mutex> lock( result_mutex );
          if( prev_spec.count(key) )
          {
            had_prev = true;
            prev = prev_spec[key];
          }
          
          prev_spec[key] = urlspec[0];
        }// end lock on result_mutex
        
        if( had_prev )
        {
          try
          {
            vector<SpecUtils::UrlSpectrum> urlsspec{ urlspec[0], prev };
            vector<QrCodeEncodedSpec> encspecs = qr_code_encode_spectra( urlsspec, QrErrorCorrection::Low, 0 );
            
            std::lock_guard<std::mutex> lock( result_mutex );
            num_qr_code_two_spec[key].push_back( encspecs.size() );
          }catch( std::exception &e )
          {
            std::lock_guard<std::mutex> lock( result_mutex );
            num_qr_code_two_spec[key].push_back( 0 );
          }
        }//
      }//if( (foreground.size() == 1) && (background.size() == 1) )
      
      
      for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
      {
        try
        {
          vector<SpecUtils::UrlSpectrum> urlspec = SpecUtils::to_url_spectra( {m}, model );
          vector<QrCodeEncodedSpec> encspecs = qr_code_encode_spectra( urlspec, QrErrorCorrection::Low, 0 );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          num_qr_code_single_spec[key].push_back( encspecs.size() );
          
          if( encspecs.size() == 1 )
          {
            qr_code_level_single_spec[key].push_back( encspecs[0].m_qr_size );
            qr_code_version_single_spec[key].push_back( encspecs[0].m_qr_version );
            qr_code_error_level_single_spec[key].push_back( static_cast<int>(encspecs[0].m_error_level) );
          }//if( encspecs.size() == 1 )
        }catch( std::exception &e )
        {
          std::lock_guard<std::mutex> lock( result_mutex );
          num_qr_code_single_spec[key].push_back( 0 );
        }
      }//for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
      
      
      for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
      {

        const vector<SpecUtils::UrlSpectrum> url_specs = SpecUtils::to_url_spectra( {m}, model );
        assert( url_specs.size() == 1 );
        
        uint8_t encode_options = 0;
        const size_t num_urls = 1;
        const unsigned int skip_encode_options = 0;
        try
        {
          const vector<string> bin_zlib_b45_urls = url_encode_spectrum( url_specs[0], encode_options, num_urls, skip_encode_options );
          assert( bin_zlib_b45_urls.size() == num_urls );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          data_sizes_bin_zlib_b45[key].push_back( bin_zlib_b45_urls[0].size() );
        }catch( std::exception & )
        {
          assert(0);
        }
        
        try
        {
          encode_options = SpecUtils::EncodeOptions::NoZeroCompressCounts;
          const vector<string> bin_nozero_urls = url_encode_spectrum( url_specs[0], encode_options, num_urls, skip_encode_options );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          data_sizes_bin_nozero[key].push_back( bin_nozero_urls[0].size() );
        }catch( std::exception & )
        {
          assert(0);
        }
        
        try
        {
          encode_options = (EncodeOptions::NoDeflate | EncodeOptions::NoBaseXEncoding);
          const vector<string> bin_urls = url_encode_spectrum( url_specs[0], encode_options, num_urls, skip_encode_options );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          data_sizes_bin[key].push_back( bin_urls[0].size() );
        }catch( std::exception & )
        {
          assert(0);
        }
        
        try
        {
          encode_options = (EncodeOptions::NoDeflate | EncodeOptions::NoBaseXEncoding | EncodeOptions::CsvChannelData);
          const vector<string> csv_urls = url_encode_spectrum( url_specs[0], encode_options, num_urls, skip_encode_options );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          data_sizes_csv[key].push_back( csv_urls[0].size() );
        }catch( std::exception & )
        {
          assert(0);
        }
        
        try
        {
          encode_options = EncodeOptions::NoBaseXEncoding;
          const vector<string> bin_zlib_urls = url_encode_spectrum( url_specs[0], encode_options, num_urls, skip_encode_options );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          data_sizes_bin_zlib[key].push_back( bin_zlib_urls[0].size() );
        }catch( std::exception & )
        {
          assert(0);
        }
        
        try
        {
          encode_options = EncodeOptions::CsvChannelData;
          const vector<string> csv_zlib_b45_urls = url_encode_spectrum( url_specs[0], encode_options, num_urls, skip_encode_options );
          
          std::lock_guard<std::mutex> lock( result_mutex );
          data_sizes_csv_zlib_b45[key].push_back( csv_zlib_b45_urls[0].size() );
        }catch( std::exception & )
        {
          assert(0);
        }
      
        
        try
        {
          const vector<SpecUtils::UrlSpectrum> start_meas = SpecUtils::to_url_spectra( {m}, model );
          
          uint8_t encode_options = 0;
          const QrErrorCorrection ecc = QrErrorCorrection::Low;
          
          vector<QrCodeEncodedSpec> encoded = QRSpectrum::qr_code_encode_spectra( start_meas, ecc, 0 );
          assert( encoded.size() >= 1 );
          
          string eccstr;
          switch( encoded[0].m_error_level )
          {
            case QrErrorCorrection::Low:      eccstr = "Low";      break;
            case QrErrorCorrection::Medium:   eccstr = "Medium";   break;
            case QrErrorCorrection::Quartile: eccstr = "Quartile"; break;
            case QrErrorCorrection::High:     eccstr = "High";     break;
          }//switch( encoded[0].m_error_level )
          
          //cout << "Encoded model: " << model << " nchannel: " << m->num_gamma_channels() << " to "
          //<< encoded.size() << " URLS, with size=" << encoded[0].m_qr_size
          //<< ", version=" << encoded[0].m_qr_version << " and ECC=" << eccstr << endl;
          
          vector<string> urls;
          for( const auto &g : encoded )
            urls.push_back( SpecUtils::url_decode( g.m_url ) );
          
#define TEST_EQUAL_ENOUGH(lhs,rhs) \
assert( (fabs((lhs) - (rhs)) < 1.0E-12) \
|| (fabs((lhs) - (rhs)) < 1.0E-4*std::max(fabs(lhs), fabs(rhs))) );
          
          try
          {
            const vector<SpecUtils::UrlSpectrum> decoded_specs = SpecUtils::decode_spectrum_urls( urls );
            //cout << "Decoded URL!" << endl;
            
            assert( decoded_specs.size() == start_meas.size() );
            for( size_t i = 0; i < decoded_specs.size(); ++i )
            {
              const auto &orig = start_meas[i];
              const auto &decoded = decoded_specs[i];
              
              assert( orig.m_source_type == decoded.m_source_type );
              assert( orig.m_energy_cal_coeffs.size() == decoded.m_energy_cal_coeffs.size() );
              for( size_t cal_index = 0; cal_index < orig.m_energy_cal_coeffs.size(); ++cal_index )
              {
                TEST_EQUAL_ENOUGH( orig.m_energy_cal_coeffs[cal_index], decoded.m_energy_cal_coeffs[cal_index] );
              }
              
              assert( orig.m_dev_pairs.size() == decoded.m_dev_pairs.size() );
              for( size_t cal_index = 0; cal_index < orig.m_dev_pairs.size(); ++cal_index )
              {
                TEST_EQUAL_ENOUGH( orig.m_dev_pairs[cal_index].first, decoded.m_dev_pairs[cal_index].first );
                TEST_EQUAL_ENOUGH( orig.m_dev_pairs[cal_index].second, decoded.m_dev_pairs[cal_index].second );
              }
              
              //Temporarily displae below check - need to limit lengths and such
              //assert( orig.m_model == decoded.m_model );
              //assert( orig.m_title == decoded.m_title );
              
              const auto tdiff = orig.m_start_time - decoded.m_start_time;
              assert( (tdiff < std::chrono::seconds(2)) && (tdiff > std::chrono::seconds(-2)) );
              
              TEST_EQUAL_ENOUGH( orig.m_latitude, decoded.m_latitude );
              TEST_EQUAL_ENOUGH( orig.m_longitude, decoded.m_longitude );
              
              assert( orig.m_neut_sum == decoded.m_neut_sum );
              
              if( orig.m_live_time > 0 )
              {
                TEST_EQUAL_ENOUGH( orig.m_live_time, decoded.m_live_time );
              }
              
              if( orig.m_real_time > 0 )
              {
                TEST_EQUAL_ENOUGH( orig.m_real_time, decoded.m_real_time );
              }
              
              assert( orig.m_channel_data.size() == decoded.m_channel_data.size() );
              for( size_t chan_index = 0; chan_index < decoded.m_channel_data.size(); ++chan_index )
              {
                assert( orig.m_channel_data[chan_index] == decoded.m_channel_data[chan_index] );
              }
            }//for( size_t i = 0; i < decoded_specs.size(); ++i )
            
            
            // Implement testing blah blah blah
            
          }catch( std::exception &e )
          {
            cerr << "Failed to decode UTR: " << e.what() << endl;
          }
          
        }catch( std::exception &e )
        {
          cerr << "Failed to encode model: " << model << " nchannel: " << m->num_gamma_channels() << " reason: " << e.what() << endl;
        }//try / catch
        
      }//for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
    };// do_work lamda
    
    pool.post( [=]{ do_work( filename ); } );
    
    thread_num += 1;
    if( thread_num >= batch_size )
    {
      pool.join();
      thread_num = 0;
    }
  }//for( string filename : files )
  
  pool.join();
  
  auto avrg_size = []( const vector<size_t> &sizes ) -> size_t {
    if( sizes.empty() )
      return 0;
    const size_t totalbytes = std::accumulate( begin(sizes), end(sizes), size_t(0) );
    return totalbytes / sizes.size();
  };
  
  // https://en.wikipedia.org/wiki/QR_code#Storage
  const size_t max_ascii_chars_v25 = 1269;
  const size_t max_ascii_chars_v40 = 1852;
  const size_t max_ascii_chars_v40L = 4296; //Low error correction (7% of data bytes can be restored)
  const size_t max_binary_v40L = 2953;
  
  const size_t header_size_single = 17;  // "RADDATA://G0/111/"
  const size_t header_size_multi = 23;   // "RADDATA://G0/111/[CRC-16 ex. 65535]/"
  
  // TODO: print out percentage that fit into max_ascii_chars_v25, max_ascii_chars_v40, max_ascii_chars_v40L
  
  for( const auto &key : data_sizes_bin_zlib_b45 )
  {
    if( key.second.size() < 20 )
      continue;
    
    const size_t nchannel = key.first.first;
    const string model = (key.first.second.size() > 15) ? key.first.second.substr(0,15) : key.first.second;
    
    auto print_stats = [=]( vector<size_t> sizes ){
      std::sort( begin(sizes), end(sizes) );
      
      size_t num_max_ascii_chars_v25 = 0, num_max_ascii_chars_v40 = 0,
             num_max_ascii_chars_v40L = 0, num_max_binary_v40L = 0;
      for( const auto i : sizes )
      {
        num_max_ascii_chars_v25 += ((i+header_size_single) < max_ascii_chars_v25);
        num_max_ascii_chars_v40 += ((i+header_size_single) < max_ascii_chars_v40);
        num_max_ascii_chars_v40L += ((i+header_size_single) < max_ascii_chars_v40L);
        num_max_binary_v40L += ((i+header_size_single) < max_binary_v40L);
      }
      
      cout
      << "," << setw(11) << sizes.front()
      << "," << setw(11) << avrg_size(sizes)
      << "," << setw(11) << sizes.back()
      << "," << setw(11) << (100.0*num_max_ascii_chars_v25) / sizes.size()
      << "," << setw(11) << (100.0*num_max_ascii_chars_v40) / sizes.size()
      << "," << setw(11) << (100.0*num_max_ascii_chars_v40L) / sizes.size()
      << "," << setw(11) << (100.0*num_max_binary_v40L) / sizes.size()
      << endl;
    };//auto print_stats = [=]( vector<size_t> sizes ){
    
    cout
    << setw(15) << "Model" << ","
    << setw(7) << "#Chan" << ","
    << setw(7) << "#Files" << ","
    << setw(17) << "Rep." << ","
    << setw(11) << "MinByte" << ","
    << setw(11) << "AvgByte" << ","
    << setw(11) << "MaxByte" << ","
    << setw(11) << "%<v25" << ","
    << setw(11) << "%<v40" << ","
    << setw(11) << "%<v40L" << ","
    << setw(11) << "%<bin-v40L"
    << endl;
    
    
    cout << setw(15) << model << "," << setw(7) << nchannel << "," << setw(7) << data_sizes_bin_zlib_b45[key.first].size() << endl;
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN+DEF+ZC+B45";
    print_stats( data_sizes_bin_zlib_b45[key.first] );
    
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN+DEF+B45";
    print_stats( data_sizes_bin_nozero[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN";
    print_stats( data_sizes_bin[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "CSV+DEF+ZC";
    print_stats( data_sizes_csv[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN+DEF+ZC";
    print_stats( data_sizes_bin_zlib[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "CSV+DEF+ZC+B45";
    print_stats( data_sizes_csv_zlib_b45[key.first] );
    
    cout << endl;
    
    // Print out stats on how many spec fit in URL, or how many QR codes it took
    auto print_num_qr_stats = [=]( vector<size_t> sizes ){
      map<size_t, size_t> nreq;
      for( auto i : sizes )
        nreq[i] += 1;
      for( auto i : nreq )
        cout << "\t" << i.first << ": " << setw(10) << ((1.0*i.second)/sizes.size());
      cout << endl;
    };
    
    cout << "\nFraction of spectra to fit within a number of QR codes:\n";
    print_num_qr_stats(num_qr_code_single_spec[key.first]);
    
    cout << "\nFraction of foreground+background in a single QR codes (" << num_qr_code_two_spec[key.first].size() << " files)" << ":\n";
    print_num_qr_stats(num_qr_code_two_spec[key.first]);
    
    
    cout << "\nFraction of single-spec-single-QR size: (" << qr_code_level_single_spec[key.first].size() << " files)" << ":\n";
    print_num_qr_stats( qr_code_level_single_spec[key.first]);
    
    cout << "\nFraction of single-spec-single-QR version: (" << qr_code_level_single_spec[key.first].size() << " files)" << ":\n";
    print_num_qr_stats( qr_code_version_single_spec[key.first]);
    
    cout << "\nFraction of single-spec-single-QR error-level: (" << qr_code_level_single_spec[key.first].size() << " files)" << ":\n";
    print_num_qr_stats( qr_code_error_level_single_spec[key.first]);
    
    
    cout << endl;
  }//for( const auto &key : data_sizes_ascii )

  return 1;
}//int dev_code()

}//namespace QRSpectrum
