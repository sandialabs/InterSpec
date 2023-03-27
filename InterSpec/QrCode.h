#ifndef QrCode_h
#define QrCode_h
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

#include <tuple>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>

// Forward declarations
class SimpleDialog;
namespace qrcodegen
{
  class QrCode;
}


namespace QrCode
{
  /** The approximate error correction level to try for, or the actual level that was used.
   The name of the enum gives the approcimate number of erroneous codewords the QR
   code can tollerate.
   */
  enum class ErrorCorrLevel
  {
    About7Percent = 0,
    About15Percent,
    About25Percent,
    About30Percent,
  };//enum class ErrorCorrLevel
  
  /** Encordes the given UTF-8 string into a SVG image.
   
   Will succeed for strings that have 2953 or fewer UTF-8 code units.
   
   @param input The URL to encode
   @param prefferedECL Preffered error correction level; encoding at the requested level may fail, in which case lower ECL will be tried.
   @param quietSpace The number of emty (or "quiet") QR code elements to place in the border around the actual QR code; the spec says this should be at least two elements.
   
   Returns a tuple with the SVG inamfe in the string, the integer number of elements (not including the border "quiet" space), as well
   as the actual used error correction level.
   
   Throws exception on error.
   */
  std::tuple<std::string,int,ErrorCorrLevel> utf8_string_to_svg_qr( const std::string &input,
                                                                   ErrorCorrLevel prefferedECL,
                                                                   const int quietSpace = 3 );

  std::pair<std::string,int> binary_to_svg_qr( const std::vector<std::uint8_t> &data );

  /** Turns a QrCode result into a SVG with a viewBox of height/width of `qr.getSize()`.
   
   The border size provides the "quiet zone" around the QR code.
   At least a 4 module wide margin around the QR code is required; see https://www.qrcode.com/en/howto/code.html
   */
  std::string to_svg_string( const qrcodegen::QrCode &qr, int border = 4 );


  /** Returns nullptr if an error is encountered. */
  SimpleDialog *displayTxtAsQrCode( const std::string &url,
                                    const std::string &title,
                                    const std::string &description );
}//namespace QrCode

#endif //QrCode_h
