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

extern "C"{
#include <zlib.h>
}

#include <boost/crc.hpp>
#include <boost/endian/conversion.hpp>

#include <Wt/Utils>

#include "QR-Code-generator/cpp/qrcodegen.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/QrCode.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/QRSpectrum.h"
#include "InterSpec/PhysicalUnits.h"


// Files only needed for `int dev_code();`
#include <mutex>
#include "SpecUtils/SpecUtilsAsync.h"

/** Most for expository purposes, a compile time option to printout
 the steps of creating a spectrum URL
 */
#define LOG_URL_ENCODING 0

#if( LOG_URL_ENCODING && !defined(_WIN32) )
#warning "Explicitly logging URL encoding steps to stdout"
#endif //#if( LOG_URL_ENCODING && !defined(_WIN32) )


using namespace std;

//Make sure assignment in UrlSpectrum definition is valid
static_assert( static_cast<int>(SpecUtils::SourceType::Unknown) == 4,
              "if SpecUtils::SourceType::Unknown is not four, please modify class definition" );

namespace
{
  const int sm_qr_quite_area = 3;
  
  const char * const sm_hex_digits = "0123456789ABCDEF";

  // From: https://datatracker.ietf.org/doc/rfc9285/ , table 1
  const char sm_base45_chars[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";

  // Here we replace the '+' and '/' characters with '-' and '_', respectively.
  //  See https://datatracker.ietf.org/doc/html/rfc4648#section-5
  //  Note: if we want to avoid URL encoding the padding character, we could replace the padding character of '=' with '.'
  const char sm_base64url_padding = '=';
  const char sm_base64url_chars[64 + 1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";

/*
  const char sm_base64url_padding = '=';
  const char sm_base64url_chars[64 + 1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
 */

  // Implement Table 1 in rfc9285 as a switch; should maybe just switch to using a lookup table
  uint8_t b45_to_dec( const char i )
  {
    switch( i )
    {
      case '0': return 0;
      case '1': return 1;
      case '2': return 2;
      case '3': return 3;
      case '4': return 4;
      case '5': return 5;
      case '6': return 6;
      case '7': return 7;
      case '8': return 8;
      case '9': return 9;
     
      // I think the letters should always be uppercase, but we'll allow lower case, jic, for the moment
      case 'A': case 'a': return 10;
      case 'B': case 'b': return 11;
      case 'C': case 'c': return 12;
      case 'D': case 'd': return 13;
      case 'E': case 'e': return 14;
      case 'F': case 'f': return 15;
      case 'G': case 'g': return 16;
      case 'H': case 'h': return 17;
      case 'I': case 'i': return 18;
      case 'J': case 'j': return 19;
      case 'K': case 'k': return 20;
      case 'L': case 'l': return 21;
      case 'M': case 'm': return 22;
      case 'N': case 'n': return 23;
      case 'O': case 'o': return 24;
      case 'P': case 'p': return 25;
      case 'Q': case 'q': return 26;
      case 'R': case 'r': return 27;
      case 'S': case 's': return 28;
      case 'T': case 't': return 29;
      case 'U': case 'u': return 30;
      case 'V': case 'v': return 31;
      case 'W': case 'w': return 32;
      case 'X': case 'x': return 33;
      case 'Y': case 'y': return 34;
      case 'Z': case 'z': return 35;
      case ' ': return 36;
      case '$': return 37;
      case '%': return 38;
      case '*': return 39;
      case '+': return 40;
      case '-': return 41;
      case '.': return 42;
      case '/': return 43;
      case ':': return 44;
      
      default:
        throw std::runtime_error( "Invalid base-45 character with decimal value "
                               + std::to_string( (int)reinterpret_cast<const uint8_t &>(i) ) );
    }//switch( i )
  
    assert( 0 );
    return 255;
  }//int b45_to_dec( char )

 uint8_t hex_to_dec( char v )
 {
   if( v >= '0' && v <= '9' )
     return static_cast<uint8_t>( v - '0' );
   
   switch( v )
   {
     case 'A': case 'a': return 0x0A;
     case 'B': case 'b': return 0x0B;
     case 'C': case 'c': return 0x0C;
     case 'D': case 'd': return 0x0D;
     case 'E': case 'e': return 0x0E;
     case 'F': case 'f': return 0x0F;
   }
   
   throw runtime_error( string("Invalid hex-digit '") + v + string("'") );
   return 0;
 }//uint8_t hex_to_dec( char v )


// We cant just use Wt::Utils::urlEncode(...) because it puts hex-values into lower-case, where
//  we need them in upper case, since QR codes alpha-numeric are upper-case only
template<class T>
string url_or_email_encode( const T &url, const string &invalid_chars )
{
  static_assert( sizeof(typename T::value_type) == 1, "Must be byte-based container" );
  
  string answer;
  answer.reserve( url.size()*2 ); //A large guess, to keep from copying a lot during resizes
  
  for( const typename T::value_type val : url )
  {
    unsigned char c = (unsigned char)val;
    
    if( (c <= 31) || (c >= 127) || (invalid_chars.find(c) != std::string::npos) )
    {
      answer += '%';
      answer += sm_hex_digits[ ((c >> 4) & 0x0F) ];
      answer += sm_hex_digits[ (c & 0x0F) ];
    }else
    {
      answer += val;
    }
  }//for( const string::value_type val : url )
  
  return answer;
}//url_encode(...)

template<class T>
string url_encode( const T &url )
{
  const std::string invalid_chars = " $&+,:;=?@'\"<>#%{}|\\^~[]`/";
  return url_or_email_encode( url, invalid_chars );
}

template<class T>
string email_encode( const T &url )
{
  const std::string invalid_chars = "%&;=/?#[]";
  return url_or_email_encode( url, invalid_chars );
}



// If we arent gzipping, and not base-45 encoding, and using channel data as ascii, then
//  we will make sure this URL can be encoded as a QR in ASCII mode, we will URL-encode
//  all non-base-45 characters.
// Note: The result of this encoding, will also get URL encoded, which is less than ideal
//       but its just a few extra bytes.
string url_encode_non_base45( const string &input )
{
  string answer;
  for( const char val : input )
  {
    if( std::find( begin(sm_base45_chars), end(sm_base45_chars), val ) == end(sm_base45_chars) )
    {
      unsigned char c = (unsigned char)val;
      answer += '%';
      answer += sm_hex_digits[ ((c >> 4) & 0x0F) ];
      answer += sm_hex_digits[ (c & 0x0F) ];
    }else
    {
      answer += val;
    }
  }//for( const char val : operator_notes )
  
  return answer;
};//url_encode_non_base45


template<class T>
std::string base45_encode_bytes( const T &input )
{
  static_assert( sizeof(typename T::value_type) == 1, "Must be byte-based container" );
  
  //From rfc9285:
  // """For encoding, two bytes [a, b] MUST be interpreted as a number n in
  // base 256, i.e. as an unsigned integer over 16 bits so that the number
  // n = (a * 256) + b.
  // This number n is converted to base 45 [c, d, e] so that n = c + (d *
  // 45) + (e * 45 * 45).  Note the order of c, d and e which are chosen
  // so that the left-most [c] is the least significant."""
  
  const size_t input_size = input.size();
  const size_t dest_bytes = 3 * (input_size / 2) + ((input_size % 2) ? 2 : 0);
  string answer( dest_bytes, ' ' );
  
  size_t out_pos = 0;
  for( size_t i = 0; i < input_size; i += 2 )
  {
    if( (i + 1) < input_size )
    {
      // We will process two bytes, storing them into three base-45 letters
      // n = c + (d * 45) + (e * 45 * 45)
      const uint16_t val_0 = reinterpret_cast<const uint8_t &>( input[i] );
      const uint16_t val_1 = reinterpret_cast<const uint8_t &>( input[i+1] );
      
      uint16_t n = (val_0 << 8) + val_1;
      
      //n may be 65535.   65535/(45 * 45)=32
      
      const uint8_t e = n / (45 * 45);
      n %= (45 * 45);
      const uint8_t d = n / 45;
      const uint8_t c = n % 45;
      
      assert( e < sizeof(sm_base45_chars) );
      assert( c < sizeof(sm_base45_chars) );
      assert( d < sizeof(sm_base45_chars) );
      
      answer[out_pos++] = reinterpret_cast<const typename T::value_type &>( sm_base45_chars[c] );
      answer[out_pos++] = reinterpret_cast<const typename T::value_type &>( sm_base45_chars[d] );
      answer[out_pos++] = reinterpret_cast<const typename T::value_type &>( sm_base45_chars[e] );
      assert( out_pos <= dest_bytes );
    }else
    {
      // We have one last dangling byte
      // a = c + (45 * d)
      const uint8_t a = reinterpret_cast<const uint8_t &>( input[i] );
      const uint8_t d = a / 45;
      const uint8_t c = a % 45;
      
      assert( c < sizeof(sm_base45_chars) );
      assert( d < sizeof(sm_base45_chars) );
      
      answer[out_pos++] = reinterpret_cast<const typename T::value_type &>( sm_base45_chars[c] );
      answer[out_pos++] = reinterpret_cast<const typename T::value_type &>( sm_base45_chars[d] );
      assert( out_pos <= dest_bytes );
    }
  }//for( size_t i = 0; i < input_size; i += 2 )
  
  assert( out_pos == dest_bytes );
  
#ifndef NDEBUG
  for( const auto val : answer )
  {
    const char c = static_cast<char>( val );
    assert( std::find( begin(sm_base45_chars), end(sm_base45_chars), c ) != end(sm_base45_chars) );
  }
#endif
  
  return answer;
}//std::string base45_encode_bytes( const vector<uint8_t> &input )


template<class T>
std::string base64url_encode_bytes( const T &input, const bool use_padding )
{
  // The bitshifting logic of this function taken from
  //  https://github.com/emweb/wt/blob/3.7-release/src/web/base64.h
  const size_t num_out_bytes = 4 * ((input.size()+2) / 3);
  
  size_t out_pos = 0;
  string answer( num_out_bytes, '\0' );
  
  for( size_t index = 0; index < input.size(); )
  {
    uint32_t value = 0;

    const size_t bytes_read = std::min( size_t(3), (input.size() - index) );
    
    for( size_t i = 0; i < bytes_read; ++i )
    {
      value <<= 8;
      value += static_cast<uint8_t>( input[index+i] );
    }
    
    index += bytes_read;
    
    // convert to base64
    int bits = static_cast<int>( bytes_read*8 );
    while( bits > 0 )
    {
      bits -= 6;
      const uint8_t index = ((bits < 0) ? value << -bits : value >> bits) & 0x3F;
      answer[out_pos++] = sm_base64url_chars[index];
    }
  }

  // add pad characters
  if( use_padding )
  {
    while( (out_pos % 4) != 0 )
      answer[out_pos++] = sm_base64url_padding;
    assert( out_pos == num_out_bytes );
  }else
  {
    assert( out_pos <= num_out_bytes );
    answer.resize( out_pos );
  }
  
  return answer;
}//std::string base64url_encode_bytes( const T &input )
  

template<class T>
void deflate_compress_internal( const void *in_data, size_t in_data_size, T &out_data)
{
  static_assert( sizeof(typename T::value_type) == 1, "Must be byte-based container" );
  
  uLongf out_len = compressBound( in_data_size );
  T buffer( out_len, 0x0 );
  
  const int rval = compress2( (Bytef *)buffer.data(),  &out_len, (const Bytef *)in_data, in_data_size, Z_BEST_COMPRESSION );
  
  if( (rval != Z_OK) || (out_len == 0) )  //Other possible values: Z_MEM_ERROR, Z_BUF_ERROR
    throw runtime_error( "Error compressing data" );
  
  buffer.resize( out_len );
  out_data.swap( buffer );
}//void deflate_compress_internal( const void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data)



template<class T>
void deflate_decompress_internal( void *in_data, size_t in_data_size, T &out_data )
{
  static_assert( sizeof(typename T::value_type) == 1, "Must be byte-based container" );
  
  z_stream zs;
  memset(&zs, 0, sizeof(zs));
  
  if( inflateInit(&zs) != Z_OK )
    throw(std::runtime_error("deflate_decompress: error from inflateInit while de-compressing."));
  
  zs.next_in = (Bytef*)in_data;
  zs.avail_in = static_cast<unsigned int>( in_data_size );
  
  int ret = Z_OK;
  typename T::value_type buffer[1024*16];
  
  T result;
  
  do
  {
    zs.next_out = reinterpret_cast<Bytef *>(buffer);
    zs.avail_out = sizeof(buffer);
    
    ret = inflate( &zs, 0 );
    
    if( result.size() < zs.total_out)
      result.insert( end(result), buffer, buffer + (zs.total_out - result.size()) );
  }while( ret == Z_OK );
  
  inflateEnd( &zs );
  
  if( ret != Z_STREAM_END )
    throw runtime_error( "deflate_decompress: Error decompressing : ("
                        + std::to_string(ret) + ") "
                        + (zs.msg ? string(zs.msg) : string()) );
  
  out_data.swap( result );
}//void deflate_decompress_internal( const void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data)

vector<uint32_t> compress_to_counted_zeros( const vector<uint32_t> &input )
{
  vector<uint32_t> results;
  
  const size_t nBin = input.size();
  
  for( size_t bin = 0; bin < nBin; ++bin )
  {
    const bool isZero = (input[bin] == 0);
    results.push_back( input[bin] );
    
    if( isZero )
    {
      uint32_t nBinZeroes = 0;
      while( ( bin < nBin ) && (input[bin] == 0) )
      {
        ++nBinZeroes;
        ++bin;
      }//while more zero bins
      
      results.push_back( nBinZeroes );
      
      if( bin != nBin )
        --bin;
    }//if( input[bin] == 0.0 )
  }//for( size_t bin = 0; bin < input.size(); ++bin )
  
  return results;
}//void compress_to_counted_zeros(...)


vector<uint32_t> zero_compress_expand( const vector<uint32_t> &input )
{
  vector<uint32_t> expanded;
  const auto dstart = begin(input);
  const auto dend = end(input);
  
  for( auto iter = dstart; iter != dend; ++iter)
  {
    if( ((*iter) != 0) || ((iter+1)==dend) )
    {
      expanded.push_back(*iter);
    }else
    {
      iter++;
      if( (*iter) == 0 )
        throw runtime_error( "Invalid counted zeros: less than one number of zeros" );
      
      const uint32_t nZeroes = *iter;
      
      if( (expanded.size() + nZeroes) > 131072 )
        throw runtime_error( "Invalid counted zeros: too many total elements" );
      
      for( uint32_t k = 0; k < nZeroes; ++k )
        expanded.push_back( 0 );
    }//if( at a non-zero value, the last value, or the next value is zero) / else
  }//for( iterate over data, iter )
  
  return expanded;
};//zero_compress_expand


string to_hex_bytes_str( const string &input )
{
  string answer;
  for( size_t i = 0; i < input.size(); ++i )
  {
    if( i )
      answer += " ";
    uint8_t v = static_cast<uint8_t>(input[i]);
    answer += sm_hex_digits[((v >> 4) & 0x0F)];
    answer += sm_hex_digits[(v & 0x0F)];
  }
  
  return answer;
}

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
      
      std::vector<UrlSpectrum> fore_url_spec = to_url_spectra( {foreground}, spec.instrument_model() );
      std::vector<UrlSpectrum> back_url_spec = to_url_spectra( {background}, spec.instrument_model() );
      std::vector<UrlSpectrum> both_url_spec = to_url_spectra( {foreground, background}, spec.instrument_model() );
      
      auto doEncode = [wanted_ecl,filename,ecl_str]( const vector<UrlSpectrum> &url_spec ) -> vector<UrlEncodedSpec> {

        return url_encode_spectra( url_spec, wanted_ecl, 0 );
        
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
        vector<UrlEncodedSpec> encoded;
        for( uint8_t option : encode_options )
        {
          try
          {
            vector<UrlEncodedSpec> this_try = url_encode_spectra( url_spec, wanted_ecl, option );
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
      
      auto toHtml = [&output_html,ecl_str,filename,base_dir,wanted_ecl]( const UrlEncodedSpec &urlspec, string tag ){
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
        //const vector<UrlEncodedSpec> encoded = url_encode_spectra( fore_url_spec, wanted_ecl, 0 );
        const vector<UrlEncodedSpec> encoded = doEncode( fore_url_spec );
        
        
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
        //const vector<UrlEncodedSpec> encoded = url_encode_spectra( back_url_spec, wanted_ecl, 0 );
        const vector<UrlEncodedSpec> encoded = doEncode( back_url_spec );
        
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
        //const vector<UrlEncodedSpec> encoded = url_encode_spectra( both_url_spec, wanted_ecl, 0 );
        const vector<UrlEncodedSpec> encoded = doEncode( both_url_spec );
        
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
std::string base45_encode( const std::vector<uint8_t> &input )
{
  return base45_encode_bytes( input );
}

std::string base45_encode( const std::string &input )
{
  return base45_encode_bytes( input );
}


vector<uint8_t> base45_decode( const string &input )
{
  const size_t input_size = input.size();
  
  if( (input_size == 1) || ((input_size % 3) && ((input_size - 2) % 3)) )
    throw runtime_error( "base45_decode: invalid input size (" + std::to_string(input_size) + ")" );
  
  const size_t output_size = ( 2*(input_size / 3) + ((input_size % 3) ? 1 : 0) );
  vector<uint8_t> answer( output_size );
  
  for( size_t i = 0, output_pos = 0; i < input_size; i += 3 )
  {
    assert( (i + 2) <= input_size );
    
    const uint32_t c = b45_to_dec( input[i] );
    const uint32_t d = b45_to_dec( input[i+1] );
    uint32_t n = c + 45 * d;
    
    if( (i + 2) < input_size )
    {
      uint32_t e = b45_to_dec( input[i+2] );
      
      n += e * 45 * 45;
      
      if( n >= 65536 )
        throw runtime_error( "base45_decode: Invalid three character sequence ("
                             + input.substr(i,3) + " -> " + std::to_string(n) + ")" );
      
      assert( (n / 256) <= 255 );
      
      answer[output_pos++] = static_cast<uint8_t>( n / 256 );
      n %= 256;
    }//if( (i + 2) < input_size )
    
    assert( n <= 255 );
    answer[output_pos++] = static_cast<uint8_t>( n );
  }//for( size_t i = 0, output_pos = 0; i < input_size; i+=3 )
  
  return answer;
}//vector<uint8_t> base45_decode( const string &input )


std::string base64url_encode( const std::string &input, const bool use_padding )
{
  return base64url_encode_bytes( input, use_padding );
}

std::string base64url_encode( const std::vector<uint8_t> &input, const bool use_padding )
{
  return base64url_encode_bytes( input, use_padding );
}

vector<uint8_t> base64url_decode( const std::string &input )
{
  // The bitshifting logic of this function taken from
  //  https://github.com/emweb/wt/blob/3.7-release/src/web/base64.h
  
  const char * const alphabet_start = sm_base64url_chars;
  const char * const alphabet_end   = sm_base64url_chars + 64;
  
  size_t output_len = 0;
  
  // If we require padding, I think the length would be `1 + 3*(input.size()/4)`,
  //  but we wont do this, so its `2 + 3*(input.size()/4)`.  For example
  //  (ex. "MTIzNDU" --decode--> "12345", but with padding the input would have
  //  been "MTIzNDU=")
  vector<uint8_t> answer( 2 + 3*(input.size()/4), 0x0 );
  
  for( size_t index = 0; index < input.size(); )
  {
    uint8_t input_bytes[4] = {0, 0, 0, 0};

    size_t num_chars = 0;
    while( (num_chars < 4) && (index < input.size()) )
    {
      const char character = static_cast<char>( input[index++] );
      if( character == sm_base64url_padding )
      {
        index = input.size();
        break;
      }
      
      const auto pos = std::find( alphabet_start, alphabet_end, character );
      if( pos != alphabet_end )
      {
        input_bytes[num_chars++] = (pos - alphabet_start);
        assert( (input_bytes[num_chars-1] < 128) && (input_bytes[num_chars-1] >= 0) );
      }else
      {
        // rfc4648 says we can choose to ignore characters out of the alphabet
        // as long as our specification says to do this, e.g.
        //  "be liberal in what you accept"
      }
    }//while( (num_chars < 4) && (index < input.size()) )

    // output the binary data
    if( num_chars >= 2 )
    {
      assert( output_len < answer.size() );
      answer[output_len++] = static_cast<uint8_t>((input_bytes[0] << 2) + (input_bytes[1] >> 4));
      if( num_chars >= 3 )
      {
        assert( output_len < answer.size() );
        answer[output_len++] = static_cast<uint8_t>((input_bytes[1] << 4) + (input_bytes[2] >> 2));
        if( num_chars >= 4 )
        {
          assert( output_len < answer.size() );
          answer[output_len++] = static_cast<uint8_t>((input_bytes[2] << 6) + input_bytes[3]);
        }//if( num_chars >= 4 )
      }//if( num_chars >= 3 )
    }else
    {
      // In principle we should maybe throw exception here; each char represents 6 bits,
      //  so we need at least two chars to represent a byte.
      // But I guess we'll just ignore when this happens
    }//if( num_chars >= 2 ) / else
  }//for( size_t index = 0; index < input.size(); )

  assert( output_len <= answer.size() );
  answer.resize( output_len );
  
  return answer;
}//std::vector<uint8_t> base64url_decode( const std::string &input )

  
/** Performs the same encoding as `streamvbyte_encode` from https://github.com/lemire/streamvbyte,
 but pre-pended with a uint16_t to give number of integer entries, has a c++ interface, and is way,
 way slower.
 */
vector<uint8_t> encode_stream_vbyte( const vector<uint32_t> &input )
{
  // Performs the same encoding as `streamvbyte_encode` from https://github.com/lemire/streamvbyte,
  //  but prepends answer with a uint16_t to give number of integer entries.
  //  This niave implementation is based on README of that project.
  
  // I this function might be okay on big-endian machines, but untested, so we'll leave a
  //  compile time assert here
  static_assert( boost::endian::order::native == boost::endian::order::little, "This function not tested in big-endian" );
  
  const size_t count = input.size();
  
  // We'll limit the size, for our implementation, because we should never be seeing greater than
  //  64k channel spectra, I think.  The leading uint16_t is the only part of this implementation
  //  that is size-limited.
  if( count > std::numeric_limits<uint16_t>::max() )
    throw runtime_error( "encode_stream_vbyte: input too large" );
  
  // Encoded data starts with ((count + 3) / 4), two bit ints
  const size_t num_cntl_bytes = (count + 3) / 4;
  
  vector<uint8_t> answer;
  answer.reserve( 2 * input.size() ); // Assumes not-super-high-statistics spectra, and is probably a bit larger than needed
  answer.resize( num_cntl_bytes + 2, 0x00 );
  
  const uint16_t ncount = static_cast<uint16_t>( count );
  answer[0] = static_cast<uint8_t>( ncount & 0x00FF );
  answer[1] = static_cast<uint8_t>( (ncount & 0xFF00) >> 8 );
  
  uint16_t test_ncount;
  memcpy( &test_ncount, &(answer[0]), 2 );
  assert( test_ncount == ncount );
  
  if( input.empty() )
    return answer;
  
  for( size_t i = 0; i < count; ++i )
  {
    const uint32_t val = input[i];
    
    uint8_t ctrl_val = 0;
    answer.push_back( static_cast<uint8_t>( val & 0x000000FF ) );
    if( val < 256 )
    {
      ctrl_val = 0;
    }else if( val < 65536 )
    {
      ctrl_val = 1;
      answer.push_back( static_cast<uint8_t>( (val & 0x0000FF00 ) >> 8  ) );
    }else if( val < 16777216 )
    {
      ctrl_val = 2;
      answer.push_back( static_cast<uint8_t>( (val & 0x0000FF00 ) >> 8  ) );
      answer.push_back( static_cast<uint8_t>( (val & 0x00FF0000 ) >> 16 ) );
    }else
    {
      ctrl_val = 3;
      answer.push_back( static_cast<uint8_t>( (val & 0x0000FF00 ) >> 8  ) );
      answer.push_back( static_cast<uint8_t>( (val & 0x00FF0000 ) >> 16 ) );
      answer.push_back( static_cast<uint8_t>( (val & 0xFF000000 ) >> 24 ) );
    }//if( we can represent in 1 byte ) / else 2 byte / 3 / 4
    
    const size_t ctrl_byte = i / 4;
    const uint8_t ctrl_shift = 2 * (i % 4);
    
    assert( ctrl_byte < num_cntl_bytes );
    
    answer[2 + ctrl_byte] |= (ctrl_val << ctrl_shift);
  }//for( size_t i = 0; i < count; ++i )
  
  return answer;
}//vector<uint8_t> encode_stream_vbyte( const vector<uint32_t> &input )


template<class T>
size_t decode_stream_vbyte( const T * const input_begin, const size_t nbytes, vector<uint32_t> &answer )
{
  //Performs the same encoding as `streamvbyte_decode` from https://github.com/lemire/streamvbyte,
  //  except assumes first two bytes gives a uint16_t for the number of integer entries.
  //  This niave implementation is based on README of that project.
  
  static_assert( sizeof(T) == 1, "Must be byte-based container" );
  // Maybe fine on big-endian, but untested
  static_assert( boost::endian::order::native == boost::endian::order::little, "This function not tested in big-endian" );
  
  if( nbytes < 2 )
    throw runtime_error( "decode_stream_vbyte: input isnt long enough to give num integers." );
  
  size_t nints = reinterpret_cast<const uint8_t &>( input_begin[1] );
  nints = nints << 8;
  nints += reinterpret_cast<const uint8_t &>( input_begin[0] );
  
  answer.resize( nints );
  
  const size_t num_cntl_bytes = (nints + 3) / 4;
  
  if( nbytes < (num_cntl_bytes + 2) )
    throw runtime_error( "decode_stream_vbyte: input isnt long enough for control bytes." );
  
  const uint8_t *begin_input = (const uint8_t *)input_begin;
  const uint8_t *ctrl_begin = begin_input + 2;
  const uint8_t *data_pos = ctrl_begin + num_cntl_bytes;
  const uint8_t * const data_end = (const uint8_t *)(input_begin + nbytes);
  
  for( size_t i = 0; i < nints; ++i )
  {
    uint32_t &val = answer[i];
    const size_t ctrl_byte = i / 4;
    const uint8_t ctrl_shift = 2 * (i % 4);
    const uint8_t ctrl_val = ((ctrl_begin[ctrl_byte] >> ctrl_shift) & 0x00003);
    assert( (ctrl_val == 0) || (ctrl_val == 1) || (ctrl_val == 2) || (ctrl_val == 3) );
    
    if( (data_pos + ctrl_val + 1) > data_end )
      throw runtime_error( "decode_stream_vbyte: input shorter than needed for specified number of integers" );
    
    val = *data_pos++;
    if( ctrl_val >= 1 )
      val += (static_cast<uint32_t>( *(data_pos++) ) << 8);
    if( ctrl_val >= 2 )
      val += (static_cast<uint32_t>( *(data_pos++) ) << 16);
    if( ctrl_val == 3 )
      val += (static_cast<uint32_t>( *(data_pos++) ) << 24);
  }//for( size_t i = 0; i < nints; ++i )
  
  return data_pos - begin_input;
}//size_t decode_stream_vbyte( const T input, size_t inlen, vector<uint32_t> &answer )


size_t decode_stream_vbyte( const std::vector<uint8_t> &inbuff, std::vector<uint32_t> &answer )
{
  return decode_stream_vbyte( inbuff.data(), inbuff.size(), answer );
}



/** Creates the "data" portion of the, split into `num_parts` separate URLs.
 */
vector<string> url_encode_spectrum( const UrlSpectrum &m,
                                    const uint8_t encode_options,
                                    const size_t num_parts,
                                    const unsigned int skip_encode_options )
{
  if( (num_parts == 0) || (num_parts > 9) )
    throw runtime_error( "url_encode_spectrum: invalid number of URL portions." );
  
  if( m.m_channel_data.size() < 1 )
    throw runtime_error( "url_encode_spectrum: invalid Measurement passed in." );
  
  if( m.m_channel_data.size() < num_parts )
    throw runtime_error( "url_encode_spectrum: more parts requested than channels." );
  
  assert( encode_options <= 0x80 );
  
  // Break out the encoding options.
  // The options given by `encode_options` are for the final message - however, we may be
  //  creating a multi-spectrum QR code, which means we dont want to do deflate, base45, or URL
  //  encoding yet
  const bool use_deflate = !(encode_options & EncodeOptions::NoDeflate);
  if( (encode_options & EncodeOptions::NoBaseXEncoding) && (encode_options & EncodeOptions::UseUrlSafeBase64) )
    throw runtime_error( "url_encode_spectrum: cant specify NoBaseXEncoding and UseUrlSafeBase64" );
  
  const bool use_baseX_encoding = (!(encode_options & EncodeOptions::NoBaseXEncoding)
                           && !(encode_options & EncodeOptions::UseUrlSafeBase64) );
  const bool use_url_safe_base64 = (encode_options & EncodeOptions::UseUrlSafeBase64);
  
  const bool skip_url_encode = (skip_encode_options & SkipForEncoding::UrlEncoding);

  const bool use_bin_chan_data = !(encode_options & EncodeOptions::CsvChannelData);
  const bool zero_compress = !(encode_options & EncodeOptions::NoZeroCompressCounts);
  
  const bool skip_encoding = (skip_encode_options & SkipForEncoding::Encoding);
  const bool skip_energy   = (skip_encode_options & SkipForEncoding::EnergyCal);
  const bool skip_model    = (skip_encode_options & SkipForEncoding::DetectorModel);
  const bool skip_gps      = (skip_encode_options & SkipForEncoding::Gps);
  const bool skip_title    = (skip_encode_options & SkipForEncoding::Title);
  
  assert( !skip_encoding || (num_parts == 1) );
  const char *comma = (use_deflate || use_url_safe_base64 || use_baseX_encoding || use_bin_chan_data) ? "," : "$";
  vector<string> answer( num_parts );
  
  string &first_url = answer[0];
  first_url.reserve( 100 + 3*m.m_channel_data.size() ); //An absolute guess, has not been checked
  
  switch( m.m_source_type )
  {
    case SpecUtils::SourceType::IntrinsicActivity: first_url += "I:I "; break;
    case SpecUtils::SourceType::Calibration:       first_url += "I:C "; break;
    case SpecUtils::SourceType::Background:        first_url += "I:B "; break;
    case SpecUtils::SourceType::Foreground:        first_url += "I:F "; break;
    case SpecUtils::SourceType::Unknown:                           break;
  }//switch( m->source_type() )
  
  const float lt = (m.m_live_time > 0.0) ? m.m_live_time : m.m_real_time;
  const float rt = (m.m_real_time > 0.0) ? m.m_real_time : m.m_live_time;
  first_url += "T:" + PhysicalUnits::printCompact(rt,6) + "," + PhysicalUnits::printCompact(lt,6) + " ";
  
  if( !skip_energy )
  {
    
    // Lower channel energy not currently implemented
    if( !m.m_energy_cal_coeffs.empty() )
    {
      first_url += "C:";
      for( size_t i = 0; i < m.m_energy_cal_coeffs.size(); ++i )
        first_url += (i ? comma : "") + PhysicalUnits::printCompact(m.m_energy_cal_coeffs[i], 7 );
      first_url += " ";
      
      const vector<pair<float,float>> &dev_pairs = m.m_dev_pairs;
      if( !dev_pairs.empty() )
      {
        first_url += "D:";
        for( size_t i = 0; i < dev_pairs.size(); ++i )
        {
          first_url += (i ? comma : "") + PhysicalUnits::printCompact(dev_pairs[i].first, 5 )
          + comma + PhysicalUnits::printCompact(dev_pairs[i].second, 5 );
        }
        first_url += " ";
      }//if( !dev_pairs.empty() )
    }//if( we have energy calibration info - that we can use )
  }//if( !skip_energy )
  
  // Lamdba to remove parts of a string that would look like a field delimiter
  //  i.e., a space, followed by a colon, followed by an upper-case letter
  //  TODO: Just replaces colon with a space, for now - perhaps we should define some 'escape' sequence
  auto remove_field_delimiters = []( string input ) -> string {
    size_t pos = 0;
    while( (pos + 1) < input.size() )
    {
      pos = input.find(":", pos + 1);
      if( pos == string::npos)
        break;
      
      const char prev_char = pos ? input[pos-1] : ' ';
      const char next_char = ((pos+1) < input.size()) ? input[pos+1] : 'A';
      
      if( (prev_char == ' ') && (next_char >= 'A') && (next_char <= 'Z') )
        input[pos] = ' ';
    }//while( (pos + 1) < input.size() )
    
    return input;
  };//remove_field_delimiters lambda
  
  if( !skip_model )
  {
    // Remove equal signs and quotes
    string det_model = remove_field_delimiters( m.m_model );
    if( det_model.size() > 30 )
      det_model = det_model.substr(0,30);
    
    SpecUtils::trim( det_model );
    if( !det_model.empty() )
    {
      if( !use_deflate && !use_baseX_encoding && !use_url_safe_base64 && !use_bin_chan_data )
        det_model = url_encode_non_base45( det_model );
      first_url += "M:" + det_model + " ";
    }
  }//if( !skip_model )
  
  if( !SpecUtils::is_special(m.m_start_time) )
  {
    // TODO: ISO string takes up 15 characters - could represent as a double, a la Microsoft time
    std::string t = SpecUtils::to_iso_string( m.m_start_time );
    const size_t dec_pos = t.find( "." );
    if( dec_pos != string::npos )
      t = t.substr(0, dec_pos);
    first_url += "P:" + t + " ";
  }//if( !SpecUtils::is_special(m->start_time()) )
  
  if( !skip_gps
     && SpecUtils::valid_longitude(m.m_longitude)
     && SpecUtils::valid_latitude(m.m_latitude) )
  {
    first_url += "G:" + PhysicalUnits::printCompact(m.m_latitude, 7)
    + comma + PhysicalUnits::printCompact(m.m_longitude, 7) + " ";
  }
  
  if( m.m_neut_sum >= 0 )
    first_url += "N:" + std::to_string(m.m_neut_sum) + " ";
  
  if( !skip_title && !m.m_title.empty() )
  {
    string operator_notes = remove_field_delimiters( m.m_title );
    if( operator_notes.size() > 60 )
      operator_notes.substr(0, 60);
    
    string remark;
    remark = operator_notes;
    
    if( !use_deflate && !use_baseX_encoding && !use_url_safe_base64 && !use_bin_chan_data )
      remark = url_encode_non_base45( operator_notes );
    
    first_url += "O:" + remark + " ";
  }//if( !skip_title && !m->title().empty() )
  
  first_url += "S:";

#if( LOG_URL_ENCODING )
  const std::string url_up_to_spectrum = first_url;
#endif //#if( LOG_URL_ENCODING )
  
  const vector<uint32_t> channel_counts = zero_compress ? compress_to_counted_zeros( m.m_channel_data )
                                                        : m.m_channel_data;
  
  // TODO: we could better split up the channel data, so the first URL gets fewer channels, since it carries the meta-information, so then the remaining QR codes could be better filled up.  Or we could try different lengths of channel data to come up with what is best to evenly distribute.
  for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
  {
    string &url_data = answer[msg_num];
    
    //size_t num_chan_first = channel_counts.size() / num_parts;
    //if( num_parts > 1 && (num_chan_first > 10) )
    //{
    //  //Assume two bytes per channel, and subtract off the number of channels all the front-matter will take ups
    //  num_chan_first -= std::min( num_chan_first - 1, first_url.size() / 2 );
    //  blah blah blah, split rest up evenly...
    //}
    
    const size_t num_channel_per_part = channel_counts.size() / num_parts;
    
    if( !num_channel_per_part )
      throw runtime_error( "url_encode_spectrum: more url-parts requested than channels after zero-compressing." );
    
    const size_t start_int_index = num_channel_per_part * msg_num;
    const size_t end_int_index = ((msg_num + 1) == num_parts)
                                        ? channel_counts.size()
                                        : start_int_index + num_channel_per_part;
    
    if( use_bin_chan_data )
    {
      vector<uint32_t> integral_counts;
      for( size_t i = start_int_index; i < end_int_index; ++i )
        integral_counts.push_back( channel_counts[i] );
      
      const uint32_t num_counts = static_cast<uint32_t>( integral_counts.size() );
      if( num_counts > 65535 )
        throw runtime_error( "url_encode_spectrum: a max of 65535 channels are supported." );
      
      const vector<uint8_t> encoded_bytes = encode_stream_vbyte( integral_counts );
      const size_t start_size = url_data.size();
      url_data.resize( start_size + encoded_bytes.size() );
      memcpy( &(url_data[start_size]), (void *)encoded_bytes.data(), encoded_bytes.size() );
      
#ifndef NDEBUG
      { // Begin quick check we can get back the original data
        vector<uint32_t> decoded;
        const size_t nbytdecoded = decode_stream_vbyte( encoded_bytes, decoded );
        assert( nbytdecoded == encoded_bytes.size() );
        assert( decoded == integral_counts );
      } // End quick check we can get back the original data
#endif
    }else
    {
      for( size_t i = start_int_index; i < end_int_index; ++i )
        url_data += ((i == start_int_index) ? "" : comma) + std::to_string( channel_counts[i] );
    }//if( use_bin_chan_data ) / else
  }//for( size_t msg_num = 0; msg_num < num_parts. ++msg_num )
  
  
  if( use_deflate && !skip_encoding )
  {
    for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
    {
      string &url_data = answer[msg_num];
      
      //cout << "During encoding, before compression: " << to_hex_bytes_str(url_data.substr(0,60)) << endl;
#ifndef NDEBUG
      const string orig = url_data;
      deflate_compress( &(url_data[0]), url_data.size(), url_data );
      
      string decompressed;
      deflate_decompress( &(url_data[0]), url_data.size(), decompressed );
      assert( orig == decompressed );
#else
      deflate_compress( &(url_data[0]), url_data.size(), url_data );
#endif
      //cout << "During encoding, after compression, before base-45: '" << to_hex_bytes_str(url_data.substr(0,60)) << "'" << endl << endl;
    
      assert( !url_data.empty() );
    }//for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
  }//if( use_deflate )
  
#if( LOG_URL_ENCODING )
  const vector<string> after_deflate = answer;
#endif //#if( LOG_URL_ENCODING )
  
  // We have just gzipped the URIs, if wanted (e.g., we have what comes "RADDATA://G0/xxx/")
  if( use_baseX_encoding && !skip_encoding )
  {
    for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
    {
#ifndef NDEBUG
      const string orig = answer[msg_num];

      if( use_url_safe_base64 )
      {
        const string sub64_encoded = base64url_encode( answer[msg_num], false );
        
        answer[msg_num] = sub64_encoded;
        
        vector<uint8_t> sb64_decoded_bytes = base64url_decode( sub64_encoded );
        assert( !sb64_decoded_bytes.empty() );
        string sub64_decoded( sb64_decoded_bytes.size(), 0x0 );
        for( size_t i = 0; i < sb64_decoded_bytes.size(); ++i )
          sub64_decoded[i] = sb64_decoded_bytes[i];
        
        if( sub64_decoded != orig )
          cerr << "\n\nSafe-Url-Base-64 Not matching:\n\tOrig='" << to_hex_bytes_str(orig) << "'" << endl
               << "\tDecs='" << to_hex_bytes_str(sub64_decoded) << "'\n\n";
        assert( sub64_decoded == orig );
      }else
      {
        answer[msg_num] = base45_encode( answer[msg_num] );
        
        vector<uint8_t> decoded_bytes = base45_decode( answer[msg_num] );
        assert( !decoded_bytes.empty() );
        string decoded( decoded_bytes.size(), 0x0 );
        //memcpy( &(decoded[0]), &(decoded_bytes[0]), decoded_bytes.size() );
        for( size_t i = 0; i < decoded_bytes.size(); ++i )
          decoded[i] = decoded_bytes[i];
        
        if( decoded != orig )
          cerr << "\n\nNot matching:\n\tOrig='" << to_hex_bytes_str(orig) << "'" << endl
          << "\tDecs='" << to_hex_bytes_str(decoded) << "'\n\n";
        assert( decoded == orig );
      }

#else
      if( use_url_safe_base64 )
      {
        answer[msg_num] = base64url_encode( answer[msg_num], false );
      }else
      {
        answer[msg_num] = base45_encode( answer[msg_num], false );
      }
#endif
      //cout << "During encoding, after  base-45: '" << to_hex_bytes_str(answer[msg_num].substr(0,60)) << "'" << endl << endl;
    }
  }//if( use_baseX_encoding )
  
  

  if( use_url_safe_base64 && !skip_encoding )
  {
    for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
    {
#ifndef NDEBUG
      const string orig = answer[msg_num];
      const string sub64_encoded = base64url_encode( answer[msg_num], false );
      
      answer[msg_num] = sub64_encoded;
      
      vector<uint8_t> sb64_decoded_bytes = base64url_decode( sub64_encoded );
      assert( !sb64_decoded_bytes.empty() );
      string sub64_decoded( sb64_decoded_bytes.size(), 0x0 );
      for( size_t i = 0; i < sb64_decoded_bytes.size(); ++i )
        sub64_decoded[i] = sb64_decoded_bytes[i];
      
      if( sub64_decoded != orig )
        cerr << "\n\nSafe-Url-Base-64 Not matching:\n\tOrig='" << to_hex_bytes_str(orig) << "'" << endl
        << "\tDecs='" << to_hex_bytes_str(sub64_decoded) << "'\n\n";
      assert( sub64_decoded == orig );
#else
      answer[msg_num] = base64url_encode( answer[msg_num], false );
#endif
    }//for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
  }//if( use_url_safe_base64 )

  
#if( LOG_URL_ENCODING )
  const vector<string> after_base45 = answer;
#endif //#if( LOG_URL_ENCODING )
  
  if( !skip_encoding && !skip_url_encode )
  {
    for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
    {
#ifndef NDEBUG
      string urlencoded = url_encode( answer[msg_num] );
      assert( Wt::Utils::urlDecode(urlencoded) == answer[msg_num] );
      answer[msg_num].swap( urlencoded );
#else
      answer[msg_num] = url_encode( answer[msg_num] );
#endif
      //cout << "During encoding, after  URL encode: '" << to_hex_bytes_str(answer[msg_num].substr(0,60)) << "'" << endl << endl;
    }//for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
  }//if( !skip_encoding )
  
#if( LOG_URL_ENCODING )
  const vector<string> after_urlencode = answer;
#endif //#if( LOG_URL_ENCODING )
  
#ifndef NDEBUG
  if( !skip_encoding && (use_baseX_encoding || (!use_deflate && !use_bin_chan_data)) )
  {
    for( const string &url : answer )
    {
      for( const char val : url )
      {
        assert( std::find( begin(sm_base45_chars), end(sm_base45_chars), val ) != end(sm_base45_chars) );
      }
    }
  }//if( use_baseX_encoding || (!use_deflate && !use_bin_chan_data) )
  

  if( !skip_encoding && use_url_safe_base64 )
  {
    for( const string &url : answer )
    {
      for( const char val : url )
      {
        assert( std::find( begin(sm_base64url_chars), end(sm_base64url_chars), val ) != end(sm_base64url_chars)
               || (val == sm_base64url_padding) );
      }
    }
  }//if( use_url_safe_base64 || (!use_deflate && !use_bin_chan_data) )
  
#endif  //#ifndef NDEBUG
  
  
#if( LOG_URL_ENCODING )
  const size_t max_nchannel_printout = 16;
  const size_t max_bytes_printout = 48;
  const size_t max_data_chars = 96;
  
  auto to_hex = []( uint8_t i ) -> std::string {
    std::stringstream stream;
    stream << "0x" << std::setfill ('0') << std::setw(2) << std::hex << static_cast<unsigned int>(i);
    return stream.str();
  };//to_hex
  
  cout << "--------------------------------------------------------------------------------\n";
  cout << "URL start:\n\t";
  const string url_start = string("RADDATA://G0/") + sm_hex_digits[encode_options & 0x0F];
  cout << "\t" << url_start << endl;
  cout << "--------------------------------------------------------------------------------\n";
  
  cout << "--------------------------------------------------------------------------------\n";
  cout << "Input spectrum info:\n";
  cout << "\tItem Type: ";
  switch( m.m_source_type )
  {
    case SpecUtils::SourceType::IntrinsicActivity: cout << "Intrinsic\n"; break;
    case SpecUtils::SourceType::Calibration:       cout << "Calibration\n"; break;
    case SpecUtils::SourceType::Background:        cout << "Background\n"; break;
    case SpecUtils::SourceType::Foreground:        cout << "Foreground\n"; break;
    case SpecUtils::SourceType::Unknown:           cout << "Unspecified\n"; break;
  }//switch( m.m_source_type )
  
  if( !skip_energy )
  {
    cout << "\tEnergy Calibration Coefficients: ";
    for( size_t i = 0; i < m.m_energy_cal_coeffs.size(); ++i )
      cout << (i ? ", " : "") << m.m_energy_cal_coeffs[i];
    cout << endl;
    
    if( !m.m_dev_pairs.empty() )
    {
      cout << "\tDeviation Pairs: ";
      for( size_t i = 0; i < m.m_dev_pairs.size(); ++i )
        cout << (i ? ", (" : "(") << m.m_dev_pairs[i].first << ", " << m.m_dev_pairs[i].second << ")";
      cout << endl;
    }//if( !m.m_dev_pairs.empty() )
  }//if( !skip_energy )
  
  if( !m.m_model.empty() && !skip_model )
    cout << "\tDetector Model: " << m.m_model << endl;
  
  if( !m.m_title.empty() && !skip_title )
    cout << "\tOperator Notes: " << m.m_title << endl;
  
  if( !SpecUtils::is_special(m.m_start_time) )
  {
    std::string t = SpecUtils::to_iso_string( m.m_start_time );
    const size_t dec_pos = t.find( "." );
    if( dec_pos != string::npos )
      t = t.substr(0, dec_pos);
    cout << "\tMeasurement Start Time: " << t << endl;
  }//if( !SpecUtils::is_special(m.m_start_time) )
  
  if( !skip_gps
     && SpecUtils::valid_longitude(m.m_longitude)
     && SpecUtils::valid_latitude(m.m_latitude) )
  {
    cout << "\tGPS coordinates: " << m.m_latitude << ", " << m.m_longitude << endl;
  }
  
  if( m.m_neut_sum >= 0 )
    cout << "\tNeutron Counts: " << m.m_neut_sum << endl;
  cout << "\tLive/Real Time: " << m.m_live_time << ", " << m.m_real_time << endl;
  
  cout << "\tChannel Counts: ";
  std::vector<uint32_t> m_channel_data;
  for( size_t i = 0; i < m.m_channel_data.size() && i < max_nchannel_printout; ++i )
    cout << (i ? "," : "") << m.m_channel_data[i];
  if( m.m_channel_data.size() > max_nchannel_printout )
    cout << "... + " << (m.m_channel_data.size() - max_nchannel_printout) << " more";
  cout << endl;
  cout << "--------------------------------------------------------------------------------\n";
  
  cout << endl << endl;
  
  cout << "--------------------------------------------------------------------------------\n";
  cout << "Initial textual representation:\n\t" << url_up_to_spectrum;
  for( size_t i = 0; i < m.m_channel_data.size() && (i < max_nchannel_printout); ++i )
    cout << (i ? "," : "") << m.m_channel_data[i];
  if( m.m_channel_data.size() > max_nchannel_printout )
    cout << "... + " << (m.m_channel_data.size() - max_nchannel_printout) << " more";
  cout << "\n";
  cout << "--------------------------------------------------------------------------------\n";
  
  if( zero_compress )
  {
    cout << endl << endl;
    cout << "--------------------------------------------------------------------------------\n";
    cout << "After zero compression:\n\t" << url_up_to_spectrum << endl;
    for( size_t i = 0; i < channel_counts.size() && (i < max_nchannel_printout); ++i )
      cout << (i ? "," : "") << channel_counts[i];
    if( channel_counts.size() > max_nchannel_printout )
      cout << "... + " << (channel_counts.size() - max_nchannel_printout) << " more";
    cout << "\n";
    cout << "--------------------------------------------------------------------------------\n";
  }//if( zero_compress )
  
  
  if( use_bin_chan_data )
  {
    const vector<uint8_t> encoded_bytes = encode_stream_vbyte( channel_counts );
    
    cout << endl << endl;
    cout << "--------------------------------------------------------------------------------\n";
    cout << "After Stream VByte bit-packing of channel counts:\n\t" << url_up_to_spectrum << "[";
    for( size_t i = 0; i < encoded_bytes.size() && i < max_bytes_printout; ++i )
      cout << (i ? "," : "") << to_hex( encoded_bytes[i] );
    if( encoded_bytes.size() > max_bytes_printout )
      cout << "... + " << (encoded_bytes.size() - max_bytes_printout) << " more bytes";
    cout << "]\n";
    cout << "--------------------------------------------------------------------------------\n";
  }//if( use_bin_chan_data )
  
  
  if( use_deflate )
  {
    cout << endl << endl;
    cout << "--------------------------------------------------------------------------------\n";
    cout << "After DEFLATE compression:" << endl;
    for( size_t url_index = 0; url_index < after_deflate.size(); ++url_index )
    {
      const string &url = after_deflate[url_index];
      
      cout << "\t";
      if( after_deflate.size() > 1 )
        cout << "Part " << url_index << ": ";
      cout << "[";
      for( size_t i = 0; i < url.size() && i < max_data_chars; ++i )
        cout << (i ? "," : "") << to_hex( static_cast<const uint8_t &>( url[i] ) );
      
      if( url.size() > max_data_chars )
        cout << "... + " << (url.size() - max_data_chars) << " more bytes";
      cout << "]\n";
    }//for( size_t url_index = 0; url_index < after_deflate.size(); ++url_index )
    cout << "--------------------------------------------------------------------------------\n";
  }//if( use_deflate )
  
  if( use_baseX_encoding || use_url_safe_base64 )
  {
    cout << endl << endl;
    cout << "--------------------------------------------------------------------------------\n";
    cout << "After base-45/url-safe-64 encoding:\n\t";

    for( size_t url_index = 0; url_index < after_base45.size(); ++url_index )
    {
      const string &url = after_base45[url_index];
      
      cout << "\t";
      if( after_base45.size() > 1 )
        cout << "Part " << url_index << ": ";
      for( size_t i = 0; i < url.size() && i < max_data_chars; ++i )
        cout << url[i];
      
      if( url.size() > max_data_chars )
        cout << "... + " << (url.size() - max_data_chars) << " more characters";
      cout << endl;
    }//for( size_t url_index = 0; url_index < after_base45.size(); ++url_index )
    cout << "--------------------------------------------------------------------------------\n";
  }//if( use_baseX_encoding )
  
  if( !skip_encoding )
  {
    cout << endl << endl;
    cout << "--------------------------------------------------------------------------------\n";
    cout << "After URL encoding:\n\t" << endl;
    for( size_t url_index = 0; url_index < after_urlencode.size(); ++url_index )
    {
      const string &url = after_urlencode[url_index];
      
      cout << "\t";
      if( after_urlencode.size() > 1 )
        cout << "Part " << url_index << ": ";
      for( size_t i = 0; i < url.size() && i < max_data_chars; ++i )
        cout << url[i];
      
      if( url.size() > max_data_chars )
        cout << "... + " << (url.size() - max_data_chars) << " more characters";
      cout << endl;
    }//for( size_t url_index = 0; url_index < after_urlencode.size(); ++url_index )
  }//if( !skip_encoding )
  cout << "--------------------------------------------------------------------------------\n";
#endif //#if( LOG_URL_ENCODING )
  
  
  return answer;
}//vector<string> url_encode_spectrum(...)


std::vector<UrlEncodedSpec> url_encode_spectra( const std::vector<UrlSpectrum> &measurements,
                                           const QrErrorCorrection minErrorCorrection,
                                           const uint8_t encode_options
                                           )
{
  if( measurements.empty() )
    throw runtime_error( "url_encode_spectra: no measurements passed in." );
  if( measurements.size() > 9 )
    throw runtime_error( "url_encode_spectra: to many measurements passed in." );
  

  if( encode_options & ~(EncodeOptions::NoDeflate | EncodeOptions::NoBaseXEncoding | EncodeOptions::CsvChannelData
                         | EncodeOptions::NoZeroCompressCounts | EncodeOptions::UseUrlSafeBase64 | EncodeOptions::AsMailToUri) )
    throw runtime_error( "url_encode_spectra: invalid option passed in - see EncodeOptions." );
  
  assert( encode_options < 0x80 );
  
  const bool use_deflate = !(encode_options & EncodeOptions::NoDeflate);
  const bool use_baseX_encoding = !(encode_options & EncodeOptions::NoBaseXEncoding);
  const bool use_bin_chan_data = !(encode_options & EncodeOptions::CsvChannelData);
  //const bool zero_compress = !(encode_options & EncodeOptions::NoZeroCompressCounts);
  const bool use_url_safe_base64 = (encode_options & EncodeOptions::UseUrlSafeBase64);
  const bool use_mailto = (encode_options & EncodeOptions::AsMailToUri);
  
  if( (encode_options & EncodeOptions::NoBaseXEncoding) && (encode_options & EncodeOptions::UseUrlSafeBase64) )
    throw runtime_error( "url_encode_spectra: invalid option passed in - cant specify NoBaseXEncoding and UseUrlSafeBase64." );
  
// It looks like qrcodegen::QrCode::encodeText(...) will default to binary if the input isnt text
//  const bool alpha_num_qr_encode = ((!use_url_safe_base64) && (use_baseX_encoding || (!use_deflate && !use_bin_chan_data)));
  
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
  

  auto make_url_start = [=]( const size_t url_num, const size_t num_parts ) -> string {
    string url_start;
    if( use_mailto )
    {
      string title_postfix;
      if( num_parts > 1 )
        title_postfix += "%20" + std::to_string(url_num) + "-" + std::to_string(num_parts);
      
      url_start = string("mailto:user@example.com?subject=spectrum" + title_postfix
                           + "&body=Spectrum%20URI%0D%0Araddata://G0/");
      //static int warn_int = 0;
      //if( warn_int++ < 4 )
      //  cerr << "Warning: not using correct email uri stuff." << endl;
      //url_start = string("RADDATA://G0/");
    }else
    {
      url_start = string("RADDATA://G0/");
    }
    
    const char more_sig_char = sm_hex_digits[(encode_options >> 4) & 0x0F];
    const char less_sig_char = sm_hex_digits[encode_options & 0x0F];
    
    url_start += ((more_sig_char == '0') ? string() : (string() + more_sig_char))
                 + (string() + less_sig_char)
                 + std::to_string(num_parts - 1) + std::to_string(url_num) + "/";
    
    return url_start;
  };//auto make_url_start lamda

  
  if( measurements.size() == 1 )
  {
    const unsigned int skip_encode_options = 0x00;
    
    const UrlSpectrum &m = measurements[0];
    
    //Now need to check that it will encode into a QR
    bool success_encoding = false;
    for( size_t num_parts = 1; num_parts < 10; ++num_parts )
    {
      vector<string> trial_urls = url_encode_spectrum( m, encode_options, num_parts, skip_encode_options );
      assert( trial_urls.size() == num_parts );
      
      urls.resize( trial_urls.size() );
      qrs.clear();
      
      string crc16;
      if( trial_urls.size() > 1 )
      {
        boost::crc_16_type crc_computer;
        
        for( const string &v : trial_urls )
        {
          const string preUrlDecode = Wt::Utils::urlDecode( v );
          crc_computer.process_bytes( (void const *)preUrlDecode.data(), preUrlDecode.size() );
        }
        
        const uint16_t crc = crc_computer.checksum();
        crc16 = std::to_string( static_cast<unsigned int>(crc) );
        crc16 += "/";
      }//if( trial_urls.size() > 1 )
      
      bool all_urls_small_enough = true;
      
      for( size_t url_num = 0; url_num < trial_urls.size(); ++url_num )
      {
        string &url = urls[url_num];
        
        url = make_url_start(url_num,num_parts) + crc16;
        
        if( use_mailto )
          url += email_encode( trial_urls[url_num] );
        else
          url += trial_urls[url_num];
        
        try
        {
          // It looks like qrcodegen::QrCode::encodeText(...) will default to binary if the input isnt text
          //if( alpha_num_qr_encode )
          //{
            qrs.push_back( qrcodegen::QrCode::encodeText( url.c_str(), ecc ) );
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
      }//for( size_t url_num = 0; url_num < trial_urls.size(); ++url_num )
      
      if( all_urls_small_enough )
      {
        success_encoding = true;
        break;
      }//if( all_urls_small_enough )
    }//for( size_t num_parts = 0; num_parts < 10; ++num_parts )
    
    if( !success_encoding )
    {
      const vector<string> urls = url_encode_spectrum( m, encode_options, 1, 0 );

      throw runtime_error( "url_encode_spectra: Failed to encode spectrum into less than 10 QR"
                           " codes at the desired error correction level (len(url)="
                           + std::to_string( urls.empty() ? size_t(0) : urls[0].size() ) + ")" );
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
  }else //
  {
    // Multiple measurements to put in single URL
    urls.resize( 1 );
    string &url = urls[0];
    url = "";
    
    qrs.clear();
    
    for( size_t meas_num = 0; meas_num < measurements.size(); ++meas_num )
    {
      const UrlSpectrum &m = measurements[meas_num];
      const UrlSpectrum &m_0 = measurements[0];
      
      unsigned int skip_encode_options = SkipForEncoding::Encoding;
      if( meas_num )
        skip_encode_options |= SkipForEncoding::DetectorModel;
      
      if( meas_num && ((m.m_energy_cal_coeffs == m_0.m_energy_cal_coeffs) && (m.m_dev_pairs == m_0.m_dev_pairs)) )
        skip_encode_options |= SkipForEncoding::EnergyCal;
      
      if( meas_num && ((m.m_latitude == m_0.m_latitude) && (m.m_longitude == m.m_longitude)) )
        skip_encode_options |= SkipForEncoding::Gps;
      
      if( meas_num && (m.m_title == m_0.m_title) )
        skip_encode_options |= SkipForEncoding::Title;
      
      vector<string> spec_url = url_encode_spectrum( m, encode_options, 1, skip_encode_options );
      assert( spec_url.size() == 1 );
      
      if( meas_num )
        url += ":0A:"; //"/G" + std::to_string(meas_num)+ "/"; //Arbitrary
      url += spec_url[0];
    }//for( size_t meas_num = 0; meas_num < measurements.size(); ++meas_num )
    
    if( use_deflate )
      deflate_compress( &(url[0]), url.size(), url );
    
    if( use_baseX_encoding )
    {
      //cout << "During encoding, before base-45 encoded: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
      
      if( use_url_safe_base64 )
      {
#ifndef NDEBUG
        const auto orig_url = url;
#endif
        //cout << "Pre url_safe_base64:\n" << url << endl << endl << endl;
        
        //printf( "Pre  url_safe_base64:\n\"" );
        //for(size_t j = 0; j < url.size(); j++)
        //{
        //  unsigned int v = reinterpret_cast<uint8_t &>( url[j] );
        //  printf( "\\x%02X", v );
        //}
        //printf( "\"\n\n" );
        
        url = base64url_encode( url, false );
        
    
        //cout << "Post UseUrlSafeBase64:\n\"" << url << "\"" << endl << endl;
        
#ifndef NDEBUG
        vector<uint8_t> sb64_decoded_bytes = base64url_decode( url );
        assert( !sb64_decoded_bytes.empty() );
        string sub64_decoded( sb64_decoded_bytes.size(), 0x0 );
        for( size_t i = 0; i < sb64_decoded_bytes.size(); ++i )
          sub64_decoded[i] = sb64_decoded_bytes[i];
        const size_t orig_size = orig_url.size();
        const size_t decoded_size = orig_url.size();
        assert( orig_size == decoded_size );
        for( size_t i = 0; i < orig_size && i < decoded_size; ++i )
        {
          uint8_t decoded = sub64_decoded[i];
          uint8_t orig = orig_url[i];
          assert( decoded == orig );
        }
        assert( sub64_decoded == orig_url );
#endif
      }else
      {
        url = base45_encode( url );
      }
      
      //cout << "During encoding, after base-45, before url-encoding: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
      
//#ifndef NDEBUG
//      const vector<uint8_t> raw = base45_decode( url );
//      assert( !raw.empty() );
//      string decoded( raw.size(), 0x0 );
//      memcpy( &(decoded[0]), &(raw[0]), raw.size() );
//      assert( decoded == url );
//#endif
    }//if( use_baseX_encoding )
    
    if( use_mailto )
    {
      //cout << "Base-x encoded:\n" << url << "\n\n\n";
      //cout << "Pre(url_encode)=" << url.size();
      //cout << "\n\t" << url << endl << endl;
      url = url_encode( url );
      //cout << ", Post(url_encode)=" << url.size();
      //cout << "\n\t" << url << endl << endl;
     
      url = email_encode( url );
      //cout << ", Post email encode=" << url.size() << endl;
    }else
    {
      url = url_encode( url );
    }
        
    // cout << "During encoding, after URL encode: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl;
    url = make_url_start(0,1) + url;
    
    //cout << "During encoding len(url)=" << url.size() << ":\n\t" << url << endl << endl << endl;
    
    try
    {
      // It looks like qrcodegen::QrCode::encodeText(...) will default to binary if the input isnt text
      //if( alpha_num_qr_encode )
      //{
        qrs.push_back( qrcodegen::QrCode::encodeText( url.c_str(), ecc ) );
      //}else
      //{
      //  vector<uint8_t> data( url.size() );
      //  memcpy( &(data[0]), &(url[0]), url.size() );
      //  qrs.push_back( qrcodegen::QrCode::encodeBinary( data, ecc ) );
      //}
    }catch( std::exception &e )
    {
      //cout << "URL Failed (size=" << url.size() << "): " /*<< url*/ << endl << endl;
      throw runtime_error( "url_encode_spectra: could not fit all the spectra into a QR code - "
                          + std::string(e.what()) );
    }
  }//if( measurements.size() == 1 ) / else
  
  assert( (urls.size() == qrs.size()) && (qrs.size() >= 1) );
  
  std::vector<UrlEncodedSpec> answer;
  for( size_t i = 0; (i < urls.size()) && (i < qrs.size()); ++i )
  {
    UrlEncodedSpec spec;
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
}//encode_spectra(...)


EncodedSpectraInfo get_spectrum_url_info( std::string url )
{
  EncodedSpectraInfo answer;
  answer.m_orig_url = url;
  
  //"RADDATA://G0/111/[CRC-16 ex. 65535]/"
  if( SpecUtils::istarts_with(url, "RADDATA://G0/") )
    url = url.substr(13);
  else if( SpecUtils::istarts_with(url, "INTERSPEC://G0/") )
    url = url.substr(15);
  else if( SpecUtils::istarts_with(url, "RADDATA:G0/") )
    url = url.substr(11);
  else
    throw runtime_error( "get_spectrum_url_info: URL does not start with 'RADDATA://G0/'" );
  
  if( url.size() < 4 )
    throw runtime_error( "get_spectrum_url_info: URL too short" );
  
  try
  {
    const bool has_email_opt = (url[3] != '/');
    if( has_email_opt )
    {
      if( (url.size() < 5) || (url[4] != '/') )
        throw runtime_error( "options not terminated with a '/' character." );
      
      answer.m_encode_options = hex_to_dec( url[0] );
      answer.m_encode_options = (answer.m_encode_options << 4);
      answer.m_encode_options += hex_to_dec( url[1] );
    }else
    {
      assert( url[3] == '/' );
      
      answer.m_encode_options = hex_to_dec( url[0] );
    }
    
    if( answer.m_encode_options
       & ~(EncodeOptions::NoDeflate | EncodeOptions::NoBaseXEncoding
           | EncodeOptions::CsvChannelData | EncodeOptions::NoZeroCompressCounts
           | EncodeOptions::UseUrlSafeBase64 | EncodeOptions::AsMailToUri) )
    {
      throw runtime_error( string("Encoding option had invalid bit set (hex digit ")
                           + url[0] + (has_email_opt ? (string() + url[1]) : string()) + ")" );
    }
   
    answer.m_number_urls = hex_to_dec( url[(has_email_opt ? 2 : 1)] ) + 1;
    if( answer.m_number_urls > 10 )
      throw std::runtime_error( "Invalid number of total URLs specified" );
    
    if( answer.m_number_urls > 1 )
    {
      answer.m_spectrum_number = hex_to_dec( url[(has_email_opt ? 3 : 2)] );
      if( answer.m_spectrum_number >= answer.m_number_urls )
        throw runtime_error( "Spectrum number larger than total number URLs" );
    }else
    {
      answer.m_num_spectra = hex_to_dec( url[(has_email_opt ? 3 : 2)] ) + 1;
      if( answer.m_num_spectra > 10 )
        throw std::runtime_error( "Invalid number of spectra in URL." );
    }
    
    url = url.substr( has_email_opt ? 5 : 4 );
  }catch( std::exception &e )
  {
    throw runtime_error( "get_spectrum_url_info: options portion (three hex digits after"
                        " the //G/) of url is invalid: " + string(e.what()) );
  }//try / catch
  
  if( answer.m_number_urls > 1 )
  {
    const size_t pos = url.find( '/' );
    if( (pos == string::npos) || (pos > 5) )
      throw runtime_error( "get_spectrum_url_info: missing or invalid CRC-16 for multi-url spectrum" );
    
    string crcstr = url.substr(0, pos);
    int crc_val;
    if( !SpecUtils::parse_int(crcstr.c_str(), crcstr.size(), crc_val) )
      throw runtime_error( "get_spectrum_url_info: CRC-16 portion of URL ('" + crcstr + "') was not int" );
      
    if( (crc_val < 0) || (crc_val > std::numeric_limits<uint16_t>::max()) )
      throw runtime_error( "get_spectrum_url_info: CRC-16 portion of URL ('" + crcstr + "') was not not a uint16_t" );
      
    answer.m_crc = static_cast<uint16_t>( crc_val );
    
    url = url.substr( pos + 1 );
  }//if( answer.m_number_urls > 1 )
  
  answer.m_raw_data = url;
  
  //cout << "Before being URL decoded: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
  //url = Wt::Utils::urlDecode( url );
  
  if( !(answer.m_encode_options & EncodeOptions::NoBaseXEncoding) )
  {
    //cout << "Before being base-45 decoded: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
    vector<uint8_t> raw;
    if( answer.m_encode_options & EncodeOptions::UseUrlSafeBase64 )
      raw = base64url_decode( url );
    else
      raw = base45_decode( url );

    assert( !raw.empty() );
    url.resize( raw.size() );
    memcpy( &(url[0]), &(raw[0]), raw.size() );
  }//if( use_baseX_encoding )
  
  if( !(answer.m_encode_options & EncodeOptions::NoDeflate) )
  {
    // cout << "Going into deflate decompress: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
    deflate_decompress( &(url[0]), url.size(), url );
  }
  
  // cout << "After deflate decompress: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
  
  answer.m_data = url;
  
  return answer;
}//EncodedSpectraInfo get_spectrum_url_info( const std::string &url )


std::vector<UrlSpectrum> spectrum_decode_first_url( const std::string &url, const EncodedSpectraInfo &info )
{
  size_t pos = url.find( " S:" );
  
  if( pos == string::npos )
    throw runtime_error( "spectrum_decode_first_url: No ' S:' marker (i.e. the channel counts) found" );
  
  const string metainfo = " " + url.substr(0, pos);
  string counts_str = url.substr( pos + 3 ); //may have additional meas after the counts
  
  // Check to make sure the metainfo doesnt have the between spectrum delimiter
  if( metainfo.find(":0A:") != string::npos ) //"/G\\d/"
    throw runtime_error( "Found the measurement delimiter in pre-information - probably means a missing channel counts." );
  
  //We'll go through and get the non-spectrum info first
  UrlSpectrum spec;
  spec.m_source_type = SpecUtils::SourceType::Unknown;
  
  pos = metainfo.find( " I:" );
  if( (pos != string::npos) && ((pos + 4) < metainfo.size()) )
  {
    switch( metainfo[pos + 3] )
    {
      case 'I': spec.m_source_type = SpecUtils::SourceType::IntrinsicActivity; break;
      case 'C': spec.m_source_type = SpecUtils::SourceType::Calibration;       break;
      case 'B': spec.m_source_type = SpecUtils::SourceType::Background;        break;
      case 'F': spec.m_source_type = SpecUtils::SourceType::Foreground;        break;
        
      default:
        throw runtime_error( "Invalid source type character following ' I:'" );
    }//switch( v )
  }//if( (pos != string::npos) && ((pos + 4) < metainfo.size()) )
  
  auto get_str_field = [&metainfo]( char val ) -> string {
    const string key = string(" ") + val + string(":");
    const size_t index = metainfo.find( key );
    if( index == string::npos )
      return "";
    
    const size_t start_pos = index + 3;
    size_t end_pos = start_pos;
    
    // Look for the next string like " X:", where X can be any upper-case letter
    while( (end_pos < metainfo.size())
          && ((metainfo[end_pos] != ' ')
              || ((end_pos >= metainfo.size()) || ( (metainfo[end_pos+1] < 'A') && (metainfo[end_pos+1] > 'Z') ))
              || (((end_pos+1) >= metainfo.size()) ||(metainfo[end_pos + 2] != ':')) ) )
    {
      ++end_pos;
    }
    
    return metainfo.substr( start_pos, end_pos - start_pos);
  };//get_str_field(...)
  
  
  string cal_str = get_str_field( 'C' );
  if( !cal_str.empty() )
  {
    SpecUtils::ireplace_all( cal_str, "$", "," );
    
    if( !SpecUtils::split_to_floats( cal_str, spec.m_energy_cal_coeffs ) )
      throw runtime_error( "Invalid CSV for energy calibration coefficients." );
    
    if( spec.m_energy_cal_coeffs.size() < 2 )
      throw runtime_error( "Not enough energy calibration coefficients." );
  }//if( pos != string::npos )
  
  string dev_pair_str = get_str_field( 'D' );
  if( !dev_pair_str.empty() )
  {
    SpecUtils::ireplace_all( dev_pair_str, "$", "," );
    
    vector<float> dev_pairs;
    if( !SpecUtils::split_to_floats( dev_pair_str, dev_pairs ) )
      throw runtime_error( "Invalid CSV for deviation pairs." );
    
    if( (dev_pairs.size() % 2) != 0 )
      throw runtime_error( "Not an even number of deviation pairs." );
    
    for( size_t i = 0; (i + 1) < dev_pairs.size(); i += 2 )
      spec.m_dev_pairs.push_back( {dev_pairs[i], dev_pairs[i+1]} );
  }//if( pos != string::npos )
  
  
  spec.m_model = get_str_field( 'M' );
  spec.m_title = get_str_field( 'O' );
  
  if( (info.m_encode_options & EncodeOptions::NoDeflate)
     && (info.m_encode_options & EncodeOptions::NoBaseXEncoding)
     && (info.m_encode_options & EncodeOptions::CsvChannelData) )
  {
    spec.m_model = Wt::Utils::urlDecode( spec.m_model );
    spec.m_title = Wt::Utils::urlDecode( spec.m_title );
  }
  
  const string neut_str = get_str_field( 'N' );
  if( !neut_str.empty() )
  {
    if( !SpecUtils::parse_int( neut_str.c_str(), neut_str.size(), spec.m_neut_sum ) )
      throw runtime_error( "Failed to parse neutron count." );
  }
  
  string lt_rt_str = get_str_field( 'T' );
  if( lt_rt_str.empty() )
    throw runtime_error( "Real and Live times not given." );
  
  SpecUtils::ireplace_all( lt_rt_str, "$", "," );
  
  vector<float> rt_lt;
  if( !SpecUtils::split_to_floats( lt_rt_str, rt_lt ) )
    throw runtime_error( "Could not parse real and live times." );
  
  if( rt_lt.size() != 2 )
    throw runtime_error( "Did not get exactly two times in Real/Live time field." );
  
  spec.m_real_time = rt_lt[0];
  spec.m_live_time = rt_lt[1];
  
  const string starttime_str = get_str_field( 'P' ); //Formatted by to_iso_string
  if( !starttime_str.empty() )
  {
    spec.m_start_time = SpecUtils::time_from_string( starttime_str );
    if( SpecUtils::is_special(spec.m_start_time) )
      throw runtime_error( "Invalid start time given (" + starttime_str + ")" );
  }//if( !starttime_str.empty() )
  
  string gps_str = get_str_field( 'G' );
  if( !gps_str.empty() )
  {
    //latitude,longitude
    vector<string> parts;
    SpecUtils::split( parts, gps_str, ",$" );
    if( parts.size() != 2 )
      throw runtime_error( "GPS does not have exactly two comma separated doubles" );
    
    if( !SpecUtils::parse_double( parts[0].c_str(), parts[0].size(), spec.m_latitude ) )
      throw runtime_error( "Invalid latitude given" );
    
    if( !SpecUtils::parse_double( parts[1].c_str(), parts[1].size(), spec.m_longitude ) )
      throw runtime_error( "Invalid longitude given" );
  }//if( !gps_str.empty() )
  

  //UrlSpectrum
  string next_spec_info;
  if( info.m_encode_options & EncodeOptions::CsvChannelData )
  {
    //if( info.m_num_spectra
    const size_t end_pos = counts_str.find(":0A:");
    if( end_pos != string::npos )
    {
      next_spec_info = counts_str.substr( end_pos + 4 );
      counts_str = counts_str.substr( end_pos );
    }
    
    SpecUtils::ireplace_all( counts_str, "$", "," );
    
    vector<long long> counts;
    if( !SpecUtils::split_to_long_longs( counts_str.c_str(), counts_str.size(), counts )
       || counts.empty() )
      throw runtime_error( "Failed to parse CSV channel counts." );
    
    spec.m_channel_data.resize( counts.size() );
    for( size_t i = 0; i < counts.size(); ++i )
      spec.m_channel_data[i] = static_cast<uint32_t>( counts[i] );
  }else
  {
    if( counts_str.size() < 2 )
      throw runtime_error( "Missing binary channel data info." );
    
    const size_t nread = decode_stream_vbyte( counts_str.data(), counts_str.size(), spec.m_channel_data );
    next_spec_info = counts_str.substr( nread );
    
    if( next_spec_info.size() > 4 )
    {
      if( !SpecUtils::istarts_with(next_spec_info, ":0A:") )
        throw runtime_error( "Left-over bytes after channel data." );
      
      if( !info.m_num_spectra )
        throw runtime_error( "Multiple spectra in URI when not specified as such" );
      
      next_spec_info = next_spec_info.substr(4);
    }//if( next_spec_info.size() > 4 )
    
    //cout << "next_spec_info.len=" << next_spec_info.length() << endl;
  }//if( info.m_encode_options & EncodeOptions::CsvChannelData ) / else
  
  
  if( !(info.m_encode_options & EncodeOptions::NoZeroCompressCounts)
     && (info.m_number_urls == 1) )
  {
    spec.m_channel_data = zero_compress_expand( spec.m_channel_data );
  }//if( !(info.m_encode_options & EncodeOptions::NoZeroCompressCounts) )
  
  vector<UrlSpectrum> answer;
  answer.push_back( spec );
  
  if( next_spec_info.length() > 1 )
  {
    vector<UrlSpectrum> the_rest = spectrum_decode_first_url( next_spec_info, info );
    answer.insert( end(answer), begin(the_rest), end(the_rest) );
  }
  
  return answer;
}//std::vector<UrlSpectrum> spectrum_decode_first_url( string url, const EncodedSpectraInfo & )


std::vector<UrlSpectrum> spectrum_decode_first_url( const std::string &url )
{
  EncodedSpectraInfo info = get_spectrum_url_info( url );
  
  if( info.m_spectrum_number != 0 )
    throw runtime_error( "spectrum_decode_first_url: URL indicates this is not first URL" );
  
  return spectrum_decode_first_url( info.m_data, info );
}//std::vector<UrlSpectrum> spectrum_decode_first_url( std::string url )


std::vector<uint32_t> spectrum_decode_not_first_url( std::string url )
{
  EncodedSpectraInfo info = get_spectrum_url_info( url );
  
  if( info.m_spectrum_number == 0 )
    throw runtime_error( "spectrum_decode_not_first_url: URL indicates it is first URL" );
  
  if( info.m_data.size() < 4 )
    throw runtime_error( "spectrum_decode_not_first_url: data too short" );
  
  vector<uint32_t> answer;
  if( info.m_encode_options & EncodeOptions::CsvChannelData )
  {
    SpecUtils::ireplace_all( info.m_data, "$", "," );
    vector<long long> cd_ll;
    if( !SpecUtils::split_to_long_longs( info.m_data.data(), info.m_data.size(), cd_ll ) )
       throw runtime_error( "spectrum_decode_not_first_url: error splitting into integer channel counts" );
    
    answer.resize( cd_ll.size() );
    for( size_t i = 0; i < cd_ll.size(); ++i )
      answer[i] = ((cd_ll[i] >= 0) ? static_cast<uint32_t>( cd_ll[i] ) : uint32_t(0));
  }else
  {
    const size_t nread = decode_stream_vbyte( info.m_data.data(), info.m_data.size(), answer );
    assert( nread == info.m_data.size() );
    
    if( nread != info.m_data.size() )
      throw runtime_error( "spectrum_decode_not_first_url: Extra unrecognized information in not-first-url" );
  }//if( CSV ) / else
  
  return answer;
}//std::vector<uint32_t> spectrum_decode_not_first_url( std::string url )


std::vector<UrlSpectrum> decode_spectrum_urls( vector<string> urls )
{
  if( urls.empty() )
    throw runtime_error( "decode_spectrum_urls: no input" );
  
  const EncodedSpectraInfo info = get_spectrum_url_info( urls[0] );
  
  if( info.m_spectrum_number != 0 )
    throw runtime_error( "decode_spectrum_urls: URL indicates this is not first URL" );
  
  vector<UrlSpectrum> spec_infos = spectrum_decode_first_url( info.m_data, info );
  
  if( (urls.size() > 1) && (spec_infos.size() > 1) )
    throw runtime_error( "decode_spectrum_urls: Multiple spectra were in first URL, but multiple URLs passed in." );
  
  assert( !spec_infos.empty() );
  if( spec_infos.empty() )
    throw logic_error( "decode_spectrum_urls: no spectra in URL." );
  
  if( urls.size() == 1 )
  {
    const UrlSpectrum &first_spec = spec_infos[0];
    
    for( size_t i = 1; i < spec_infos.size(); ++i )
    {
      UrlSpectrum &spec = spec_infos[i];
      if( spec.m_model.empty() )
        spec.m_model = first_spec.m_model;
      
      if( spec.m_energy_cal_coeffs.empty()
         && (spec.m_channel_data.size() == first_spec.m_channel_data.size()) )
        spec.m_energy_cal_coeffs = first_spec.m_energy_cal_coeffs;
      
      if( spec.m_dev_pairs.empty()
         && (spec.m_channel_data.size() == first_spec.m_channel_data.size()) )
        spec.m_dev_pairs = first_spec.m_dev_pairs;
      
      if( SpecUtils::valid_latitude(first_spec.m_latitude) && !SpecUtils::valid_latitude(spec.m_latitude) )
        spec.m_latitude = first_spec.m_latitude;
      if( SpecUtils::valid_longitude(first_spec.m_longitude) && !SpecUtils::valid_longitude(spec.m_longitude) )
        spec.m_longitude = first_spec.m_longitude;
    }//for( size_t i = 1; i < spec_infos.size(); ++i )
  }else
  {
    UrlSpectrum &spec = spec_infos[0];
    for( size_t i = 1; i < urls.size(); ++i )
    {
      const vector<uint32_t> more_counts = spectrum_decode_not_first_url( urls[i] );
      spec.m_channel_data.insert( end(spec.m_channel_data), begin(more_counts), end(more_counts) );
    }
    
    if( !(info.m_encode_options & EncodeOptions::NoZeroCompressCounts) )
      spec.m_channel_data = zero_compress_expand( spec.m_channel_data );
  }//if( urls.size() == 1 ) / else
  
  
  return spec_infos;
}//shared_ptr<SpecUtils::SpecFile> decode_spectrum_urls( vector<string> urls )



vector<UrlSpectrum> to_url_spectra( vector<shared_ptr<const SpecUtils::Measurement>> specs, string det_model )
{
  vector<UrlSpectrum> answer;
  for( size_t i = 0; i < specs.size(); ++i )
  {
    shared_ptr<const SpecUtils::Measurement> spec = specs[i];
    assert( spec );
    if( !spec || (spec->num_gamma_channels() < 1) )
      throw runtime_error( "to_url_spectra: invalid Measurement" );
    
    UrlSpectrum urlspec;
    urlspec.m_source_type = spec->source_type();
    
    shared_ptr<const SpecUtils::EnergyCalibration> cal = spec->energy_calibration();
    assert( cal );
    // Lower channel energy not currently implemented
    if( (cal->type() != SpecUtils::EnergyCalType::LowerChannelEdge)
       && (cal->type() != SpecUtils::EnergyCalType::InvalidEquationType) )
    {
      vector<float> coefs = cal->coefficients();
      switch( cal->type() )
      {
        case SpecUtils::EnergyCalType::FullRangeFraction:
          coefs = SpecUtils::fullrangefraction_coef_to_polynomial( coefs, cal->num_channels() );
          break;
          
        case SpecUtils::EnergyCalType::LowerChannelEdge:
        case SpecUtils::EnergyCalType::InvalidEquationType:
          assert( 0 );
          break;
          
        case SpecUtils::EnergyCalType::Polynomial:
        case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
          break;
      }//switch( cal->type() )
      
      urlspec.m_energy_cal_coeffs = coefs;
      urlspec.m_dev_pairs = cal->deviation_pairs();
    }//if( poly or FWF )
      
    
    urlspec.m_model = det_model;
    
    string operator_notes = spec->title();
    SpecUtils::ireplace_all( operator_notes, ":", " " );
    if( operator_notes.size() > 60 )
      operator_notes.substr(0, 60);
    
    // TODO: look for this info in the "Remarks" - I think Ortec Detectives and some other models will end up getting user input to there maybe?
    urlspec.m_title = operator_notes;
    urlspec.m_start_time = spec->start_time();
    
    if( spec->has_gps_info() )
    {
      urlspec.m_latitude = spec->latitude();
      urlspec.m_longitude = spec->longitude();
    }
    
    if( spec->contained_neutron() )
    {
      const float neut_sum = static_cast<float>( spec->neutron_counts_sum() );
      urlspec.m_neut_sum = SpecUtils::float_to_integral<int>( neut_sum );
    }
    
    
    urlspec.m_live_time = spec->live_time();
    urlspec.m_real_time = spec->real_time();
    
    const vector<float> &counts = *spec->gamma_counts();
    urlspec.m_channel_data.resize( counts.size() );
    for( size_t i = 0; i < counts.size(); ++i )
      urlspec.m_channel_data[i] = static_cast<uint32_t>( counts[i] );
    
    answer.push_back( urlspec );
  }//for( size_t i = 0; i < specs.size(); ++i )
  
  return answer;
}//vector<UrlSpectrum> to_url_spectra( vector<shared_ptr<const SpecUtils::Measurement>> specs, string det_model )


std::shared_ptr<SpecMeas> to_spec_file( const std::vector<UrlSpectrum> &spec_infos )
{
  auto specfile = make_shared<SpecMeas>();
  specfile->set_instrument_model( spec_infos[0].m_model );
  
  shared_ptr<SpecUtils::EnergyCalibration> first_cal;
  for( size_t spec_index = 0; spec_index < spec_infos.size(); ++spec_index )
  {
    const UrlSpectrum &spec = spec_infos[spec_index];
    
    auto m = make_shared<SpecUtils::Measurement>();
    m->set_source_type( spec.m_source_type );
    m->set_start_time( spec.m_start_time );
    m->set_position( spec.m_longitude, spec.m_latitude, {} );
    m->set_title( spec.m_title );
    
    if( spec.m_neut_sum >= 0 )
      m->set_neutron_counts( { static_cast<float>(spec.m_neut_sum) } );
    
    auto counts = make_shared<vector<float>>();
    counts->insert( end(*counts), begin(spec.m_channel_data), end(spec.m_channel_data) );
    
    m->set_gamma_counts( counts, spec.m_live_time, spec.m_real_time );
    
    if( !spec.m_energy_cal_coeffs.empty() )
    {
      // Check to see if we can re-use energy cal from previous measurement
      if( spec_index
         && (spec.m_energy_cal_coeffs == first_cal->coefficients())
         && (spec.m_dev_pairs == first_cal->deviation_pairs())
         && (spec.m_channel_data.size() == first_cal->num_channels()) )
      {
        m->set_energy_calibration( first_cal );
      }else
      {
        first_cal = make_shared<SpecUtils::EnergyCalibration>();
        
        try
        {
          first_cal->set_polynomial( spec.m_channel_data.size(), spec.m_energy_cal_coeffs, spec.m_dev_pairs );
        }catch( std::exception &e )
        {
          throw runtime_error( "Energy cal given is invalid: " + string(e.what()) );
        }
        
        m->set_energy_calibration( first_cal );
      }//if( we can re-use energy cal ) / else
    }//if( !spec.m_energy_cal_coeffs.empty() )
    
    
    specfile->add_measurement( m, false );
  }//for( const UrlSpectrum &spec : spec_infos )
  
  specfile->cleanup_after_load();
  
  return specfile;
}//std::shared_ptr<SpecUtils::SpecFile> to_spec_file( const std::vector<UrlSpectrum> &meas );



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
  map<pair<size_t,string>,UrlSpectrum> prev_spec;
  
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
          vector<UrlSpectrum> urlspec = to_url_spectra( {foreground[0], background[0]}, model );
          vector<UrlEncodedSpec> encspecs = url_encode_spectra( urlspec, QrErrorCorrection::Low, 0 );
          
          {
            std::lock_guard<std::mutex> lock( result_mutex );
            num_qr_code_two_spec[key].push_back( encspecs.size() );
          }
          
          assert( encspecs.size() == 1 );
          //cout << "Fit two spectra in URL:\n\t" << encspecs[0].m_url << endl << endl << encspecs[0].m_qr_svg << endl << endl << endl;
          
          try
          {
            vector<UrlSpectrum> decoded = decode_spectrum_urls( { Wt::Utils::urlDecode(encspecs[0].m_url) } );
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
        vector<UrlSpectrum> urlspec = to_url_spectra( {usable_spectra[0]}, model );
        assert( urlspec.size() == 1 );
        
        bool had_prev = false;
        UrlSpectrum prev;
        
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
            vector<UrlSpectrum> urlsspec{ urlspec[0], prev };
            vector<UrlEncodedSpec> encspecs = url_encode_spectra( urlsspec, QrErrorCorrection::Low, 0 );
            
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
          vector<UrlSpectrum> urlspec = to_url_spectra( {m}, model );
          vector<UrlEncodedSpec> encspecs = url_encode_spectra( urlspec, QrErrorCorrection::Low, 0 );
          
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

        const vector<UrlSpectrum> url_specs = to_url_spectra( {m}, model );
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
          encode_options = EncodeOptions::NoZeroCompressCounts;
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
          const vector<UrlSpectrum> start_meas = to_url_spectra( {m}, model );
          
          uint8_t encode_options = 0;
          const QrErrorCorrection ecc = QrErrorCorrection::Low;
          
          vector<UrlEncodedSpec> encoded = QRSpectrum::url_encode_spectra( start_meas, ecc, 0 );
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
            urls.push_back( Wt::Utils::urlDecode( g.m_url ) );
          
#define TEST_EQUAL_ENOUGH(lhs,rhs) \
assert( (fabs((lhs) - (rhs)) < 1.0E-12) \
|| (fabs((lhs) - (rhs)) < 1.0E-4*std::max(fabs(lhs), fabs(rhs))) );
          
          try
          {
            const vector<UrlSpectrum> decoded_specs = decode_spectrum_urls( urls );
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



void deflate_compress( const void *in_data, size_t in_data_size, std::string &out_data )
{
  deflate_compress_internal( in_data, in_data_size, out_data );
}

void deflate_compress( const void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data )
{
  deflate_compress_internal( in_data, in_data_size, out_data );
}

void deflate_decompress( void *in_data, size_t in_data_size, std::string &out_data )
{
  deflate_decompress_internal( in_data, in_data_size, out_data );
}

void deflate_decompress( void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data )
{
  deflate_decompress_internal( in_data, in_data_size, out_data );
}


}//namespace QRSpectrum
