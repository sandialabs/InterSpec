#ifndef QRSpecDev_H
#define QRSpecDev_H
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

#include <string>
#include <vector>
#include <cstdint>

// Forward declarations
class SpecMeas;
namespace SpecUtils
{
  class Measurement;
  enum class SourceType : int;
}

/** The \c QRSpectrum namespace provides functions to encode spectra into URLs and QR codes.
 
 An example of using this library could be:
 \code{.cpp}
 QRSpectrum::UrlSpectrum foreground;
 foreground.m_source_type = SpecUtils::SourceType::Foreground;
 foreground.m_energy_cal_coeffs = {0.0f, 3.0f};
 foreground.m_model = "SomeDetector";
 foreground.m_title = "User entered Notes";
 foreground.m_start_time = std::chrono::system_clock::now();
 foreground.m_latitude = 10.1223;
 foreground.m_longitude = -13.232;
 foreground.m_neut_sum = 5;
 foreground.m_live_time = 295.1f;
 foreground.m_real_time = 300.0f;
 foreground.m_channel_data = {0,0,0,5,6,22,15,... };
 
 QRSpectrum::UrlSpectrum background;
 // No need to set energy cal, model, or lat/lon for background - will assume same as foreground
 background.m_source_type = SpecUtils::SourceType::Background;
 background.m_neut_sum = 1;
 background.m_live_time = 299.1;
 background.m_real_time = 300.0f;
 foreground.m_channel_data = {0,0,0,3,4,2,1,...};
 
 const vector<UrlEncodedSpec> encoded = url_encode_spectra( {foreground, background}, QRSpectrum::QrErrorCorrection::Low );
 assert( encoded.size() == 1 );
 cout << "URL: " << encoded[0].m_url << endl;
 cout << "SVG: " << encoded[0].m_qr_svg << endl;
 cout << "Num Elements: " << encoded[0].m_qr_size << endl;
 \endcode
 
 */

namespace QRSpectrum
{
  int dev_code();

  std::string base45_encode( const std::string &input );
  std::string base45_encode( const std::vector<uint8_t> &input );

  std::vector<uint8_t> base45_decode( const std::string &input );

// TODO: need to make tests for url_safe_base64 encode/decode!
  /** base64url is a URL and filename safe variant of base64 encoding.
   See: https://datatracker.ietf.org/doc/html/rfc4648#section-5
   and: https://en.wikipedia.org/wiki/Base64#Implementations_and_history
   */
  std::string base64url_encode( const std::string &input, const bool use_padding );
  std::string base64url_encode( const std::vector<uint8_t> &input, const bool use_padding );
  std::vector<uint8_t> base64url_decode( const std::string &input );
  
  void deflate_compress( const void *in_data, size_t in_data_size, std::string &out_data );
  void deflate_compress( const void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data );
  
  void deflate_decompress( void *in_data, size_t in_data_size, std::string &out_data );
  void deflate_decompress( void *in_data, size_t in_data_size, std::vector<uint8_t> &out_data );

enum class QrErrorCorrection
{
  Low,       // ~7% missing codewords
  Medium,   // ~15% missing codewords
  Quartile, // ~25% missing codewords
  High      // ~30% missing codewords
};//enum class QrErrorCorrection

/** Options for how the data can be encoded as a URL.*/
enum EncodeOptions
{
  /** Do not apply zlib based DEFLATE compression. */
  NoDeflate = 0x01,
  
  /** Do not encode data (after optional DEFLATE) as base-45 or as url-safe-base-64.
   E.g., for clickable URLs you may do this, or if you will encode to a binary QR-code, but
   will generally lead to longer URIs or larger QR codes.
   */
  NoBaseXEncoding = 0x02,
  
  /** Keep channel data as text-based numbers.  Will actually be separated by the '$'
   sign, since this is a base-45 character and commas arent - but commas would be valid
   to.
   */
  CsvChannelData = 0x04,
  
  /** Do not zero-compress channel data. */
  NoZeroCompressCounts = 0x08,
  
  /** Use a modified version of base-64 encoding (that doesnt have and invalid URL characters), instead of base-45.
   
   If this option is specified, then #NoBaseXEncoding must not be specified.
   
   This option is useful if you are generating a binary QR code, but want all printable characters; e.g.,
   when you create a QR code with a 'mailto:'.
   */
  UseUrlSafeBase64 = 0x10,
  
  /** If specified, the URI will be one to generate a email (i.e., a 'mailto:...' URI), instead of a 'raddata://' URI.
   This effects both the begining of the URI starting with 'mailto:', but also the URL data is email encoded
   according to RFC 6068, (i.e., only the "%&;=/?#[]" characters are escaped), instead of total URL
   encoding (i.e., the " $&+,:;=?@'\"<>#%{}|\\^~[]`/" characters escaped).
   
   For example, a returned URI will start with:
   "mailto:user@example.com?subject=spectrum&body=Spectrum%20URI%0D%0Araddata:..."
   and it will not be strictly URL encoded.
   
   Note that this option will not be written/indicated in the URI; however, this value 0x20, was being written
   until 20231202, so we will ignore this bit when it is in URIs, for the time being, and in the future if more
   options are added, consider skipping this bit.
   */
  AsMailToUri = 0x20,
};//enum EncodeOptions


struct UrlEncodedSpec
{
  std::string m_url;
  std::string m_qr_svg;
  int m_qr_size = 0;
  int m_qr_version = 0;
  
  QrErrorCorrection m_error_level;
};//struct UrlEncodedSpec

/** Struct that represents information that can be included in a spectrum URL. */
struct UrlSpectrum
{
  SpecUtils::SourceType m_source_type = SpecUtils::SourceType(4); // There is static_assert in .cpp file to make sure SpecUtils::SourceType::Unknown == 4
  std::vector<float> m_energy_cal_coeffs;
  std::vector<std::pair<float,float>> m_dev_pairs;
  
  std::string m_model;
  std::string m_title;
  
  // Maybe add serial number?
  
  // Use same time definition as SpecUtils.
  using time_point_t = std::chrono::time_point<std::chrono::system_clock,std::chrono::microseconds>;
  time_point_t m_start_time = time_point_t{}; // defaults to 0
  
  double m_latitude = -999.9;
  double m_longitude = -999.9;
  int m_neut_sum = -1; ///< Neutron count; will be negative if not present
  
  float m_live_time = -1;
  float m_real_time = -1;
  std::vector<uint32_t> m_channel_data;
};//struct UrlSpectrum

std::vector<UrlSpectrum> to_url_spectra(
                                  std::vector<std::shared_ptr<const SpecUtils::Measurement>> specs,
                                  std::string detector_model );

std::shared_ptr<SpecMeas> to_spec_file( const std::vector<UrlSpectrum> &meas );

/** Encodes the specified measurements into one or more URL and QR codes.
 
 If a single spectrum is passed in, you may get back multiple urls/qr-codes.
 If multiple spectra are passed in, you will get back a single url/qr-code, if they can all fit
 in a single QR code (otherwise an exception will be thrown).
 
 Note: URLs will be URL-encoded, while #get_spectrum_url_info, #decode_spectrum_urls, and similar
       all expect non-URL-encoded URLs.
 
 @param measurements One or measurement to put in the result
 @param detector_model The detector model (e.x., "Detective-X") to include.  May be blank.
 @param minErrorCorrection The minimum error correction to require.  The QR code size will be
        adjusted based on this.  Or if a higher correction can be achieved at the same level, the
        higher level will be used.
 @param encode_options Options, set by #EncodeOptions bits, for how the encoding should be done.
        Default value of 0x00 is binary channel data, DEFLATE compressed, and BASE-45 and URL
        encoded, so as to be sent as an ASCII QR code.
 
 Throws exception on error.
 */
std::vector<UrlEncodedSpec> url_encode_spectra( const std::vector<UrlSpectrum> &measurements,
                                          const QrErrorCorrection minErrorCorrection = QrErrorCorrection::Low,
                                           const uint8_t encode_options = 0
                                               );

/** Options to adjust how the URL for a spectrum is created.
 For a single spectrum going into one or more QR codes, you will not specify any of these.
 For creating a QR code with multiple spectra, you will specify to skip encoding, and you
 may be able to skip energy cal, detector model, GPS, and title for the second spectrum.
 */
enum SkipForEncoding
{
  Encoding = 0x01,      ///< Skip DEFLATE, base-45/url-safe-base-64, and URL encoding (e.g., for placing multiple spectra in a URL, where you will do this after combination)
  EnergyCal = 0x02,
  DetectorModel = 0x04,
  Gps = 0x08,
  Title = 0x10,
  UrlEncoding = 0x20   ///< Skip just URL encoding (e.g., for when you are making a 'mailto:' uri with an invalid-uri body
};//enum SkipForEncoding

/** Puts the specified spectrum into `num_parts` URLs.
 Will happily return a very long URL that cant fit into a QR code.
 
 @param meas The measurement to be encoded.
 @param det_model The detector model (ex. "Detective-x") to include.
 @param encode_options How the final URL will be encoded, as specified by #EncodeOptions.
        Note that #skip_encode_options may over-ride these options for when you are combining
        multiple spectra into a single URL.
 @param num_parts How many URLs the result should be broken into.  This is done by breaking the
        channel data into
 @param skip_encode_options Options given by #SkipForEncoding for what encoding or information not
        to do/include.  Useful for including multiple spectra in a single URL.
 
 Throws exception on error.
 */
std::vector<std::string> url_encode_spectrum( const UrlSpectrum &meas,
                                   const uint8_t encode_options,
                                   const size_t num_parts,
                                   const unsigned int skip_encode_options );


struct EncodedSpectraInfo
{
  /** See #EncodeOptions for bit meanings. */
  uint8_t m_encode_options = 0;
  uint8_t m_num_spectra = 0;
  /** starts at 0. */
  uint8_t m_spectrum_number = 0;
  uint8_t m_number_urls = 0;
  
  uint16_t m_crc = 0;

  /** The original URL, before any manipulation.
   
   Note: will not be URL encoded.
   */
  std::string m_orig_url;
  
  /** The spectrum relevant data that wont be URL encoded, but may be base-45 encoded, as well
   as Deflate compressed.
   
   This is needed as the CRC-16 is computed from all the parts, after these steps, but before
   URL encoding (because URL encoding is not unique, and the OS may have .
   */
  std::string m_raw_data;
  
  /** The spectrum relevant data, after un-base-45/un-base-85 encoded, and un-Deflated, if applicable.
   */
  std::string m_data;
};//struct EncodedSpectraInfo


/** Breaks out the various information from the url, and un-base-45 encodes, and un-Deflates,
 if necessary.
 
 Expects the input url to have already been url-decoded.
 
 Throws exception on error.
 */
EncodedSpectraInfo get_spectrum_url_info( std::string url );


std::vector<UrlSpectrum> spectrum_decode_first_url( const std::string &url );
std::vector<uint32_t> spectrum_decode_not_first_url( std::string url );

/** Decodes the given urls to a single spectrum (one or more input urls), or
 multiple spectrum (a single input url).
 
 Expects each urls to have already been url-decoded.
 */
std::vector<UrlSpectrum> decode_spectrum_urls( std::vector<std::string> urls );


/** Performs the same encoding as `streamvbyte_encode` from https://github.com/lemire/streamvbyte,
 but pre-pended with a uint16_t to give number of integer entries, has a c++ interface, and is way,
 way slower, but is a safer in terms of buffer overflow.
 */
std::vector<uint8_t> encode_stream_vbyte( const std::vector<uint32_t> &input );

/** Performs the same encoding as `streamvbyte_decode` from https://github.com/lemire/streamvbyte,
but assumes data is prepended with uin16_t that gives the number of integer entries, has a c++
 interface, is way, way slower, but is a safer in terms of buffer overflow.
 */
size_t decode_stream_vbyte( const std::vector<uint8_t> &inbuff, std::vector<uint32_t> &answer );
}//namespace QRSpectrum

#endif  //QRSpecDev_H
