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



using namespace std;
namespace
{
  const char * const sm_hex_digits = "0123456789ABCDEF";

  // From: https://datatracker.ietf.org/doc/rfc9285/ , table 1
  const char sm_base45_chars[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
  
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
string url_encode( const T &url )
{
  static_assert( sizeof(typename T::value_type) == 1, "Must be byte-based container" );
  
  const std::string invalid_chars = " $&+,:;=?@'\"<>#%{}|\\^~[]`/";
  
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
  
  assert( encode_options < 16 );
  
  // Break out the encoding options.
  // The options given by `encode_options` are for the final message - however, we may be
  //  creating a multi-spectrum QR code, which means we dont want to do deflate, base45, or URL
  //  encoding yet
  const bool use_deflate = !(encode_options & EncodeOptions::NoDeflate);
  const bool use_base45 = !(encode_options & EncodeOptions::NoBase45);
  const bool use_bin_chan_data = !(encode_options & EncodeOptions::CsvChannelData);
  const bool zero_compress = !(encode_options & EncodeOptions::NoZeroCompressCounts);
  
  const bool skip_encoding = (skip_encode_options & SkipForEncoding::Encoding);
  const bool skip_energy   = (skip_encode_options & SkipForEncoding::EnergyCal);
  const bool skip_model    = (skip_encode_options & SkipForEncoding::DetectorModel);
  const bool skip_gps      = (skip_encode_options & SkipForEncoding::Gps);
  const bool skip_title    = (skip_encode_options & SkipForEncoding::Title);
  
  assert( !skip_encoding || (num_parts == 1) );
  const char *comma = (use_deflate || use_base45 || use_bin_chan_data) ? "," : "$";
  
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
  
  
  if( !skip_model )
  {
    // Remove equal signs and quotes
    string det_model = m.m_model;
    SpecUtils::ireplace_all( det_model, ":", "" );
    if( det_model.size() > 30 )
      det_model = det_model.substr(0,30);
    
    SpecUtils::trim( det_model );
    if( !det_model.empty() )
    {
      if( !use_deflate && !use_base45 && !use_bin_chan_data )
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
    string operator_notes = m.m_title;
    SpecUtils::ireplace_all( operator_notes, ":", " " );
    if( operator_notes.size() > 60 )
      operator_notes.substr(0, 60);
    
    string remark;
    remark = operator_notes;
    
    if( !use_deflate && !use_base45 && !use_bin_chan_data )
      remark = url_encode_non_base45( operator_notes );
    
    first_url += "O:" + remark + " ";
  }//if( !skip_title && !m->title().empty() )
  
  first_url += "S:";
  
  vector<uint32_t> channel_counts;
  if( zero_compress )
    channel_counts = compress_to_counted_zeros( m.m_channel_data );
  else
    channel_counts = m.m_channel_data;
  
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
  
  if( use_base45 && !skip_encoding )
  {
    for( size_t msg_num = 0; msg_num < num_parts; ++msg_num )
    {
#ifndef NDEBUG
      const string orig = answer[msg_num];
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
#else
      answer[msg_num] = base45_encode( answer[msg_num] );
#endif
      //cout << "During encoding, after  base-45: '" << to_hex_bytes_str(answer[msg_num].substr(0,60)) << "'" << endl << endl;
    }
  }//if( use_base45 )
  
  if( !skip_encoding )
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
  
#ifndef NDEBUG
  if( !skip_encoding && (use_base45 || (!use_deflate && !use_bin_chan_data)) )
  {
    for( const string &url : answer )
    {
      for( const char val : url )
      {
        assert( std::find( begin(sm_base45_chars), end(sm_base45_chars), val ) != end(sm_base45_chars) );
      }
    }
  }//if( use_base45 || (!use_deflate && !use_bin_chan_data) )
#endif
  
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
  
  if( encode_options & ~(EncodeOptions::NoDeflate | EncodeOptions::NoBase45 | EncodeOptions::CsvChannelData | EncodeOptions::NoZeroCompressCounts) )
    throw runtime_error( "url_encode_spectra: invalid option passed in - see EncodeOptions." );
  
  assert( encode_options < 16 );
  
  const bool use_deflate = !(encode_options & EncodeOptions::NoDeflate);
  const bool use_base45 = !(encode_options & EncodeOptions::NoBase45);
  const bool use_bin_chan_data = !(encode_options & EncodeOptions::CsvChannelData);
  const bool zero_compress = !(encode_options & EncodeOptions::NoZeroCompressCounts);
  
  const bool alpha_num_qr_encode = (use_base45 || (!use_deflate && !use_bin_chan_data));
  
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
  
  const string url_start = string("RADDATA://G0/") + sm_hex_digits[encode_options & 0x0F];
  
  if( measurements.size() == 1 )
  {
    const UrlSpectrum &m = measurements[0];
    
    //Now need to check that it will encode into a QR
    bool success_encoding = false;
    for( size_t num_parts = 1; num_parts < 10; ++num_parts )
    {
      vector<string> trial_urls = url_encode_spectrum( m, encode_options, num_parts, 0x00 );
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
        
        url = url_start + std::to_string(num_parts - 1)
        + std::to_string(url_num) + "/" + crc16
        + trial_urls[url_num];
        
        try
        {
          if( alpha_num_qr_encode )
          {
            qrs.push_back( qrcodegen::QrCode::encodeText( url.c_str(), ecc ) );
          }else
          {
            vector<uint8_t> data( url.size() );
            memcpy( &(data[0]), &(url[0]), url.size() );
            qrs.push_back( qrcodegen::QrCode::encodeBinary( data, ecc ) );
          }
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
    
    if( use_base45 )
    {
      //cout << "During encoding, before base-45 encoded: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
      
      url = base45_encode( url );
      
      //cout << "During encoding, after base-45, before url-encoding: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
      
//#ifndef NDEBUG
//      const vector<uint8_t> raw = base45_decode( url );
//      assert( !raw.empty() );
//      string decoded( raw.size(), 0x0 );
//      memcpy( &(decoded[0]), &(raw[0]), raw.size() );
//      assert( decoded == url );
//#endif
    }//if( use_base45 )
    
#ifndef NDEBUG
    string urlencoded = url_encode( url );
    assert( Wt::Utils::urlDecode(urlencoded) == url );
    url.swap( urlencoded );
#else
    url = url_encode( url );
#endif
    
    // cout << "During encoding, after URL encode: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl;
    url = url_start + "0" + std::to_string( measurements.size() - 1 ) + "/" + url;
    
    try
    {
      if( alpha_num_qr_encode )
      {
        qrs.push_back( qrcodegen::QrCode::encodeText( url.c_str(), ecc ) );
      }else
      {
        vector<uint8_t> data( url.size() );
        memcpy( &(data[0]), &(url[0]), url.size() );
        qrs.push_back( qrcodegen::QrCode::encodeBinary( data, ecc ) );
      }
    }catch( std::exception &e )
    {
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
    spec.m_qr_svg = QrCode::to_svg_string( qrs[i], 1 );
    
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
  else
    throw runtime_error( "get_spectrum_url_info: URL does not start with 'RADDATA://G0/'" );
  
  if( url.size() < 4 )
    throw runtime_error( "get_spectrum_url_info: URL too short" );
  
  try
  {
    answer.m_encode_options = hex_to_dec( url[0] );
    
    if( answer.m_encode_options
       & ~(EncodeOptions::NoDeflate | EncodeOptions::NoBase45
           | EncodeOptions::CsvChannelData | EncodeOptions::NoZeroCompressCounts) )
    {
      throw runtime_error( string("Encoding option had invalid bit set (hex digit ") + url[0] + ")" );
    }
    
    answer.m_number_urls = hex_to_dec( url[1] ) + 1;
    if( answer.m_number_urls > 10 )
      throw std::runtime_error( "Invalid number of total URLs specified" );
    
    if( answer.m_number_urls > 1 )
    {
      answer.m_spectrum_number = hex_to_dec( url[2] );
      if( answer.m_spectrum_number >= answer.m_number_urls )
        throw runtime_error( "Spectrum number larger than total number URLs" );
    }else
    {
      answer.m_num_spectra = hex_to_dec( url[2] ) + 1;
      if( answer.m_num_spectra > 10 )
        throw std::runtime_error( "Invalid number of spectra in URL." );
    }
    
    if( url[3] != '/' )
      throw runtime_error( "options not followed by a '/' character." );
    
    url = url.substr( 4 );
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
  
  if( !(answer.m_encode_options & EncodeOptions::NoBase45) )
  {
    //cout << "Before being base-45 decoded: '" << to_hex_bytes_str(url.substr(0,60)) << "'" << endl << endl;
    const vector<uint8_t> raw = base45_decode( url );
    assert( !raw.empty() );
    url.resize( raw.size() );
    memcpy( &(url[0]), &(raw[0]), raw.size() );
  }//if( use_base45 )
  
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
     && (info.m_encode_options & EncodeOptions::NoBase45)
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
         && (spec.m_channel_data.size() == first_spec.m_energy_cal_coeffs.size()) )
        spec.m_energy_cal_coeffs = first_spec.m_energy_cal_coeffs;
      
      if( spec.m_dev_pairs.empty()
         && (spec.m_channel_data.size() == first_spec.m_energy_cal_coeffs.size()) )
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
  // Now that URL encoding has solidified - this function should be cleaned up to use that code
  //  to compute statistics of how many QR codes it typically takes and such.
#define DELETE_UNWANTED_FILES 0
#define SAVE_ASCII_OUTPUT 0
#define USE_ZSTDLIB_CL 0

  const char *base_dir = "/Users/wcjohns/rad_ana/qrspec_test_data";
#if( SAVE_ASCII_OUTPUT )
  const char *out_dir = "/Users/wcjohns/rad_ana/processed_qrspec";
#endif
  
  const vector<string> files = SpecUtils::recursive_ls(base_dir);
  
  map<pair<size_t,string>,UrlSpectrum> prev_spec;
  map<pair<size_t,string>, vector<size_t>> data_sizes_ascii, data_sizes_raw_bin,
                                           data_sizes_zlib, data_sizes_zlib_url, data_sizes_ascii_zlib,
                                           data_sizes_ascii_zlib_url, data_sizes_ascii_zlib_base_45_url,
                                           data_sizes_bin_base45, data_sizes_bin_base45_url,
                                           num_qr_code_single_spec, num_qr_code_two_spec;
  
#if( USE_ZSTDLIB_CL )
  map<pair<size_t,string>, vector<size_t>> data_sizes_zstdlib;
#endif
    
  
  for( string filename : files )
  {
    if( SpecUtils::likely_not_spec_file(filename) )
    {
#if( DELETE_UNWANTED_FILES )
      SpecUtils::remove_file( filename );
#endif
      continue;
    }
    
    SpecUtils::SpecFile spec;
    if( !spec.load_file( filename, SpecUtils::ParserType::Auto, filename) )
    {
#if( DELETE_UNWANTED_FILES )
      SpecUtils::remove_file( filename );
#endif
      continue;
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
      continue;
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
        num_qr_code_two_spec[key].push_back( encspecs.size() );
        
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
        num_qr_code_two_spec[key].push_back( 0 );
      }
    }else if( usable_spectra.size() == 1 )
    {
      vector<UrlSpectrum> urlspec = to_url_spectra( {usable_spectra[0]}, model );
      assert( urlspec.size() == 1 );
      
      if( prev_spec.count(key) )
      {
        try
        {
          vector<UrlSpectrum> urlspec{ urlspec[0], prev_spec[key] };
          vector<UrlEncodedSpec> encspecs = url_encode_spectra( urlspec, QrErrorCorrection::Low, 0 );
          num_qr_code_two_spec[key].push_back( encspecs.size() );
        }catch( std::exception &e )
        {
          num_qr_code_two_spec[key].push_back( 0 );
        }
      }
      
      prev_spec[key] = urlspec[0];
    }
    
    for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
    {
      try
      {
        vector<UrlSpectrum> urlspec = to_url_spectra( {m}, model );
        vector<UrlEncodedSpec> encspecs = url_encode_spectra( urlspec, QrErrorCorrection::Low, 0 );
        num_qr_code_single_spec[key].push_back( encspecs.size() );
      }catch( std::exception &e )
      {
        num_qr_code_single_spec[key].push_back( 0 );
      }
    }//for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
    
    
    //EncodeOptions::NoDeflate
    //EncodeOptions::NoBase45
    //EncodeOptions::CsvChannelData
    //EncodeOptions::NoZeroCompressCounts
    //SkipForEncoding::Encoding
    //SkipForEncoding::EnergyCal
    //SkipForEncoding::DetectorModel
    //SkipForEncoding::Gps
    //SkipForEncoding::Title
    
    for( shared_ptr<const SpecUtils::Measurement> m : usable_spectra )
    {
      stringstream strm;
      
      switch( m->source_type() )
      {
        case SpecUtils::SourceType::IntrinsicActivity: strm << "I=I "; break;
        case SpecUtils::SourceType::Calibration:       strm << "I=C "; break;
        case SpecUtils::SourceType::Background:        strm << "I=B "; break;
        case SpecUtils::SourceType::Foreground:        strm << "I=F "; break;
        case SpecUtils::SourceType::Unknown:                           break;
      }//switch( m->source_type() )
      
      const float lt = (m->live_time() > 0.0) ? m->live_time() : m->real_time();
      const float rt = (m->real_time() > 0.0) ? m->real_time() : m->live_time();
      strm << "T:" << PhysicalUnits::printCompact(rt,6) << "," << PhysicalUnits::printCompact(lt,6) << " ";
      
      shared_ptr<const SpecUtils::EnergyCalibration> cal = m->energy_calibration();
      assert( cal->type() == SpecUtils::EnergyCalType::Polynomial );
      strm << "C:";
      for( size_t i = 0; i < cal->coefficients().size(); ++i )
        strm << (i ? "," : "") << PhysicalUnits::printCompact(cal->coefficients()[i], i ? 4 : 7 );
      strm << " ";
      
      const vector<pair<float,float>> &dev_pairs = cal->deviation_pairs();
      if( !dev_pairs.empty() )
      {
        strm << "D:";
        for( size_t i = 0; i < dev_pairs.size(); ++i )
        {
          strm << (i ? "," : "") << PhysicalUnits::printCompact(dev_pairs[i].first, 5 )
          << "," << PhysicalUnits::printCompact(dev_pairs[i].second, 4 );
        }
        strm << " ";
      }//if( !dev_pairs.empty() )
      
      
      if( !model.empty() )
        strm << "M:" << model << " ";
      
      if( !SpecUtils::is_special(m->start_time()) )
      {
        // TODO: ISO string takes up 15 characters - could represent as a double, a la Microsoft time
        std::string t = SpecUtils::to_iso_string(m->start_time());
        const size_t dec_pos = t.find( "." );
        if( dec_pos != string::npos )
          t = t.substr(0, dec_pos);
        strm << "P:" << t << " ";
      }//if( !SpecUtils::is_special(m->start_time()) )
      
      if( m->has_gps_info() )
      {
        strm << "G:" << PhysicalUnits::printCompact(m->latitude(), 7)
             << "," << PhysicalUnits::printCompact(m->longitude(), 7) << " ";
      }
      
      if( m->contained_neutron() )
      {
        const float neut_sum = static_cast<float>( m->neutron_counts_sum() );
        const int ineut_sum = SpecUtils::float_to_integral<int>( neut_sum );
        strm << "N:" << ineut_sum << " ";
      }
      
      if( !m->title().empty() )
      {
        // TODO: look for this info in the "Remarks" - I think Ortec Detectives and some other models will end up getting user input to there maybe?
        string operator_notes = m->title();
        SpecUtils::ireplace_all( operator_notes, ":", " " );
        if( operator_notes.size() > 60 )
          operator_notes.substr(0, 60);
        
        string remark;
        remark = operator_notes;
        
        // If we arent gzipping, and not base-45 encoding, and channel data as ascii, then we should make sure this is ascii so the QR code can be generated in ascii mode - however, URL encoding here is wont work, because we will URL encode things again later on...
        //for( const char val : operator_notes )
        //{
        //  if( std::find( begin(sm_base45_chars), end(sm_base45_chars), val ) == end(sm_base45_chars) )
        //  {
        //    unsigned char c = (unsigned char)val;
        //    remark += '%';
        //    remark += sm_hex_digits[ ((c >> 4) & 0x0F) ];
        //    remark += sm_hex_digits[ (c & 0x0F) ];
        //  }else
        //  {
        //    remark += val;
        //  }
        //}
        
        strm << "O:" << remark << " ";
      }//if( !m->title().empty() )
      
      strm << "S:";
      const string data_up_to_spectrum = strm.str();
      
      vector<float> zero_compressed_counts;
      SpecUtils::compress_to_counted_zeros( *m->gamma_counts(), zero_compressed_counts );
      
      for( size_t i = 0; i < zero_compressed_counts.size(); ++i )
        strm << (i ? "," : "") << SpecUtils::float_to_integral<int>( zero_compressed_counts[i] );
      
      vector<int32_t> signed_compressed_integral_counts;
      vector<uint32_t> compressed_integral_counts;
      for( const float f : zero_compressed_counts )
      {
        compressed_integral_counts.push_back( SpecUtils::float_to_integral<uint32_t>( f ) );
        signed_compressed_integral_counts.push_back( SpecUtils::float_to_integral<int32_t>( f ) );
      }
      
      const string ascii_data = strm.str();
      
      if( model.empty() )
        model = "other";
      
#if( SAVE_ASCII_OUTPUT )
      const string save_dir = SpecUtils::append_path(out_dir, model);
      if( !SpecUtils::is_directory(save_dir) )
      {
        if( !SpecUtils::create_directory(save_dir) )
        {
          cerr << "Failed to make directory '" << save_dir << "'" << endl;
          assert( 0 );
        }
      }//if( !SpecUtils::is_directory(save_dir) )

      const string data_hash = Wt::Utils::hexEncode( Wt::Utils::md5(ascii_data) );
      const string save_file_name = SpecUtils::append_path(save_dir, data_hash + ".spec.txt" );
      
      ofstream output( save_file_name.c_str(), ios::out | ios::binary );
      assert( output.is_open() );
      output.write( ascii_data.c_str(), ascii_data.size() );
      assert( output.good() );
#endif
      
      if( compressed_integral_counts.size() > 65535 )
        throw runtime_error( "A max of 65535 is supported." );
      
      const vector<uint8_t> encoded_bytes = encode_stream_vbyte( compressed_integral_counts );
      
      /*
      auto print_test_case = [=](){
        assert( encoded_bytes == encode_stream_vbyte( compressed_integral_counts ) );
        
        static int test_num = 1;
        const auto inflags = cout.flags();
        
        cout << "\n\n\n  // Test case " << test_num << endl;
        cout << "  const vector<uint32_t> test_" << test_num << "_chan_cnts{ ";
        for( size_t i = 0; i < compressed_integral_counts.size(); ++i )
        {
          cout << (i ? ", " : "");
          if( !(i % 20) )
            cout << "\n    ";
          cout << compressed_integral_counts[i];
        }
        cout << "  };\n";
        cout << "  assert( test_" << test_num << "_chan_cnts.size() == " << compressed_integral_counts.size() << " );\n";
        
        cout << "  const vector<uint8_t> test_" << test_num << "_packed{ ";
        
        for( size_t i = 0; i < encoded_bytes.size(); ++i )
        {
          cout << (i ? ", " : "");
          if( !(i % 50) )
            cout << "\n    ";
          
          unsigned int val = encoded_bytes[i];
          //cout << "0x" << ios::hex << ios::uppercase << val;
          cout << val;
        }
        cout << "\n  };\n";
        cout << "  assert( test_" << test_num << "_packed.size() == " << encoded_bytes.size() << " );\n";
        
        cout << "  const vector<uint8_t> test_" << test_num << "_encoded = QRSpectrum::encode_stream_vbyte( test_" << dec << test_num << "_chan_cnts );\n"
        "  assert( test_" << dec << test_num << "_encoded == test_" << test_num << "_packed );\n"
        "  vector<uint32_t> test_" << test_num << "_dec;\n"
        "  const size_t test_"<< test_num << "_nbytedec = QRSpectrum::decode_stream_vbyte(test_" << test_num << "_encoded,test_" << test_num << "_dec);\n"
        "  assert( test_" << test_num << "_nbytedec == test_" << test_num << "_packed.size() );\n"
        "  assert( test_" << test_num << "_dec == test_" << dec << test_num << "_chan_cnts );\n"
        << endl;
        
        cout.flags(inflags);
        test_num += 1;
      };
      
      if( m->gamma_counts()->size() == 512 )
      {
        static int nprinted512 = 0;
        if( nprinted512++ < 5 )
          print_test_case();
      }else if( m->gamma_counts()->size() == 1024 )
      {
        static int nprinted1k = 0;
        if( nprinted1k++ < 5 )
          print_test_case();
      }else if( m->gamma_counts()->size() == 8192 )
      {
        static int nprinted8k = 0;
        if( nprinted8k++ < 5 )
          print_test_case();
      }else if( m->gamma_counts()->size() == 16384 )
      {
        static int nprinted16k = 0;
        if( nprinted16k++ < 5 )
          print_test_case();
      }
       */
      
      {// Begin quick check we can recover things
        vector<uint32_t> recovdata;
        const size_t nread = decode_stream_vbyte( encoded_bytes.data(), encoded_bytes.size(), recovdata );
        assert( nread == encoded_bytes.size() );
        assert( recovdata == compressed_integral_counts );
      }// End quick check we can recover things
      
      
      
      /*
      //No bitpacking
      std::vector<uint8_t> encoded_bytes( 4*compressed_integral_counts.size() );
      memcpy( (void *)encoded_bytes.data(), (void *)compressed_integral_counts.data(), 4*compressed_integral_counts.size() );
      */
      
      /*
#define ENCODE_BITS 7
#define ENCODE_FREF 0
      using codec = oroch::bitfor_codec<uint32_t>;
      std::vector<uint8_t> encoded_bytes( codec::basic_codec::space(compressed_integral_counts.size(),ENCODE_BITS)  );
      codec::parameters params(ENCODE_FREF, ENCODE_BITS);
      oroch::dst_bytes_t enc_dest = (unsigned char *)encoded_bytes.data();
      codec::encode( enc_dest, begin(compressed_integral_counts), end(compressed_integral_counts), params );
      encoded_bytes.resize( enc_dest - (unsigned char *)encoded_bytes.data() );
      
      {// Begin quick check we can decode
        codec::parameters params(ENCODE_FREF, ENCODE_BITS);
        vector<uint32_t> output( compressed_integral_counts.size() );
        
        oroch::src_bytes_t b_it = encoded_bytes.data();
        codec::decode( begin(output), end(output), b_it, params);
        
        for (int i = 0; i < compressed_integral_counts.size(); i++) {
          if( compressed_integral_counts[i] != output[i] )
       throw runtime_error(
        }
      }// End check we can decode
      */
      
      
      //TODO:
      //  - Try https://github.com/lemire/LittleIntPacker
      //  - try out lzma compression
      vector<uint8_t> raw_bin_data( data_up_to_spectrum.size() + encoded_bytes.size() );
      
      
      
      memcpy( raw_bin_data.data(), data_up_to_spectrum.data(), data_up_to_spectrum.size() );
      memcpy( raw_bin_data.data() + data_up_to_spectrum.size(), encoded_bytes.data(), encoded_bytes.size() );
      
#if( SAVE_ASCII_OUTPUT )
      const string bin_save_file_name = SpecUtils::append_path(save_dir, data_hash + ".spec.bin" );
      
      {// Begin block to write bin output
        ofstream bin_output( bin_save_file_name.c_str(), ios::out | ios::binary );
        assert( bin_output.is_open() );
        bin_output.write( (const char *)raw_bin_data.data(), raw_bin_data.size() );
        assert( bin_output.good() );
      }// End block to write bin output
#endif
      
#if( USE_ZSTDLIB_CL )
      vector<char> plain_zstdlib_data;
      {// Begin block to do zstd compression, and read back in
        const string out_stdzlib_name = bin_save_file_name + ".zstd";
        const string cmd = "/usr/local/bin/zstd --ultra -22 -q -D '/Users/wcjohns/rad_ana/processed_qrspec/zstd_spec.dict' '" + bin_save_file_name + "' -o '" + out_stdzlib_name + "'";
        const int rval = system( cmd.c_str() );
        assert( rval == 0 );
        SpecUtils::load_file_data( out_stdzlib_name.c_str(), plain_zstdlib_data );
        SpecUtils::remove_file( out_stdzlib_name );
      }// End block to do zstd compression, and read back in
#endif
      
      // Need to generate a dictionary from all files
      //vector<char> dict_zstdlib_data;
      //system( ("cd '" + save_dir + "'; zstd --ultra -22 " + data_hash + ".spec.bin").c_str() );
      //SpecUtils::load_file_data( (bin_save_file_name + ".zst").c_str(), dict_zstdlib_data );
      //SpecUtils::remove_file( bin_save_file_name + ".zst" );
      
      
      vector<uint8_t> zlib_data, ascii_zlib_data;
      deflate_compress( raw_bin_data.data(), raw_bin_data.size(), zlib_data );
      deflate_compress( (const void *)&(ascii_data[0]), ascii_data.size(), ascii_zlib_data );
      
      {
        vector<uint8_t> decomp_out_data;
        deflate_decompress( zlib_data.data(), zlib_data.size(), decomp_out_data );
        assert( decomp_out_data == raw_bin_data );
        
        deflate_decompress( ascii_zlib_data.data(), ascii_zlib_data.size(), decomp_out_data );
        vector<uint8_t> ascii_in( (uint8_t *)ascii_data.data(), (uint8_t *)(ascii_data.data() + ascii_data.size()) );
        assert( decomp_out_data == ascii_in );
      }
      
      
      //Now add in base45 encoding, then URL encoding.
      const string base45_data = base45_encode( zlib_data );
      const string url_base45_data = url_encode( base45_data );
      const string zlib_url = url_encode( zlib_data );
      
      const string ascii_base45_data = base45_encode( ascii_zlib_data );
      const string ascii_url_base45_data = url_encode( ascii_base45_data );
      
      const string ascii_url = url_encode( ascii_data );
      
      const size_t nchannel = m->num_gamma_channels();
      
      data_sizes_ascii[key].push_back( ascii_url.size() );
      data_sizes_raw_bin[key].push_back( raw_bin_data.size() );
      data_sizes_zlib[key].push_back( zlib_data.size() );
      data_sizes_zlib_url[key].push_back( zlib_url.size() );
      
#if( USE_ZSTDLIB_CL )
      data_sizes_zstdlib[key].push_back( plain_zstdlib_data.size() );
#endif
      
      data_sizes_ascii_zlib[key].push_back( ascii_zlib_data.size() );
      data_sizes_bin_base45[key].push_back( base45_data.size() );
      data_sizes_bin_base45_url[key].push_back( url_base45_data.size() );
      data_sizes_ascii_zlib_url[key].push_back( ascii_zlib_data.size() );
      data_sizes_ascii_zlib_base_45_url[key].push_back( ascii_url_base45_data.size() );
      
      
      //Add in statistics for how much zero compress helps, then how much encoding helps
      //data_up_to_spectrum
      //vector<float> zero_compressed_counts;
      //vector<uint32_t> compressed_integral_counts;

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
  }//for( string filename : files )
  
  
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
  
  for( const auto &key : data_sizes_ascii )
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
    
    cout << setw(15) << model << "," << setw(7) << nchannel
         << "," << setw(7) << data_sizes_ascii[key.first].size()
         << endl;
    
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "ASCII URL";
    print_stats( data_sizes_ascii[key.first] );
    
    //cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "RAW BIN";
    //print_stats( data_sizes_raw_bin[key.first] );
    
    //cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << ","<< setw(17) << "BIN ZLIB";
    //print_stats( data_sizes_zlib[key.first] );
    
#if( USE_ZSTDLIB_CL )
    //cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << ","<< setw(17) << "ZSTDLIB";
    //print_stats( data_sizes_zstdlib[key.first] );
#endif
    
    //cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "ASCII ZLIB";
    //print_stats( data_sizes_ascii_zlib_url[key.first] );
    
    //cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN B45 ZLIB";
    //print_stats( data_sizes_bin_base45[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN B45URL ZLIB";
    print_stats( data_sizes_bin_base45_url[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "ASCII B45URL ZLIB";
    print_stats( data_sizes_ascii_zlib_base_45_url[key.first] );
    
    cout << setw(15) << "" << "," << setw(7) << "" << "," << setw(7) << "" << "," << setw(17) << "BIN ZLIB URL";
    print_stats( data_sizes_zlib_url[key.first] );
    
    cout << endl << endl;
    
    // Print out stats on how many spec fit in URL, or how many QR codes it took
    auto print_num_qr_stats = [=]( vector<size_t> sizes ){
      map<size_t, size_t> nreq;
      for( auto i : sizes )
        nreq[i] += 1;
      for( auto i : nreq )
        cout << "\t" << i.first << ": " << setw(10) << ((1.0*i.second)/sizes.size());
      cout << endl;
    };
    
    cout << "\nPercentage of spectra to fit within a number of QR codes:\n";
    print_num_qr_stats(num_qr_code_single_spec[key.first]);
    
    cout << "\nPercentage of foreground+background in a single QR codes (" << num_qr_code_two_spec[key.first].size() << " files)" << ":\n";
    print_num_qr_stats(num_qr_code_two_spec[key.first]);
    
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
