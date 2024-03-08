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

#include "SpecUtils/UriSpectrum.h"

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
 SpecUtils::UrlSpectrum foreground;
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
 
 SpecUtils::UrlSpectrum background;
 // No need to set energy cal, model, or lat/lon for background - will assume same as foreground
 background.m_source_type = SpecUtils::SourceType::Background;
 background.m_neut_sum = 1;
 background.m_live_time = 299.1;
 background.m_real_time = 300.0f;
 foreground.m_channel_data = {0,0,0,3,4,2,1,...};
 
 const vector<QRSpectrum::QrCodeEncodedSpec> encoded
           = qr_code_encode_spectra( {foreground, background}, QRSpectrum::QrErrorCorrection::Low );
 assert( encoded.size() == 1 );
 cout << "URL: " << encoded[0].m_url << endl;
 cout << "SVG: " << encoded[0].m_qr_svg << endl;
 cout << "Num Elements: " << encoded[0].m_qr_size << endl;
 \endcode
 
 */

namespace QRSpectrum
{
  int dev_code();

enum class QrErrorCorrection
{
  Low,       // ~7% missing codewords
  Medium,   // ~15% missing codewords
  Quartile, // ~25% missing codewords
  High      // ~30% missing codewords
};//enum class QrErrorCorrection


struct QrCodeEncodedSpec
{
  std::string m_url;
  std::string m_qr_svg;
  int m_qr_size = 0;
  int m_qr_version = 0;
  
  QrErrorCorrection m_error_level;
};//struct QrCodeEncodedSpec

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
 @param encode_options Options, set by #SpecUtils::EncodeOptions bits, for how the encoding should be done.
        Default value of 0x00 is binary channel data, DEFLATE compressed, and BASE-45 and URL
        encoded, so as to be sent as an ASCII QR code.
 
 Throws exception on error.
 */
std::vector<QrCodeEncodedSpec> qr_code_encode_spectra( const std::vector<SpecUtils::UrlSpectrum> &measurements,
                                          const QrErrorCorrection minErrorCorrection = QrErrorCorrection::Low,
                                           const uint8_t encode_options = 0
                                               );
}//namespace QRSpectrum

#endif  //QRSpecDev_H
