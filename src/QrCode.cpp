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
#include <sstream>
#include <iostream>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WAnchor>
#include <Wt/WApplication>
#include <Wt/WPushButton>
#include <Wt/WMemoryResource>
#include <Wt/WContainerWidget>
#include <Wt/WJavaScriptPreamble>

#include "QR-Code-generator/cpp/qrcodegen.hpp"

#include "InterSpec/QrCode.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"

using namespace Wt;
using namespace std;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


namespace
{

/** Converts a QR code to a SVG string we can include in the DOM.
 
 Adapted from QrCodeGeneratorDemo.cpp that is part of QR Code generator library.
 */
std::string to_svg_string( const qrcodegen::QrCode &qr, int border )
{
  if (border < 0)
    throw std::domain_error("Border must be non-negative");
  if (border > INT_MAX / 2 || border * 2 > INT_MAX - qr.getSize())
    throw std::overflow_error("Border too large");
  
  std::ostringstream sb;
  //sb << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  //sb << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
  sb << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"0 0 ";
  sb << (qr.getSize() + border * 2) << " " << (qr.getSize() + border * 2) << "\" stroke=\"none\">\n";
  sb << "\t<rect width=\"100%\" height=\"100%\" fill=\"#FFFFFF\"/>\n";
  sb << "\t<path d=\"";
  for (int y = 0; y < qr.getSize(); y++) {
    for (int x = 0; x < qr.getSize(); x++) {
      if (qr.getModule(x, y)) {
        if (x != 0 || y != 0)
          sb << " ";
        sb << "M" << (x + border) << "," << (y + border) << "h1v1h-1z";
      }
    }
  }
  sb << "\" fill=\"#000000\"/>\n";
  sb << "</svg>\n";
  return sb.str();
}//std::string to_svg_string( const qrcodegen::QrCode &qr, int border )


WT_DECLARE_WT_MEMBER
 (CopyUrlToClipboard, Wt::JavaScriptFunction, "CopyUrlToClipboard",
  function( sender, event, id, domSignalId, text )
{
  const success_msg = 'showMsg-info-Copied URL to clipboard';
  const fail_msg = 'showMsg-error-Failed to copy URL to clipboard';
  
  if( !navigator.clipboard )
  {
    let textArea = document.createElement("textarea");
    textArea.value = text;
    temparea.style.position = "fixed";
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    
    let successful = false;
    try
    {
      successful = document.execCommand('copy');
    }catch( ex )
    {
      successful = false;
    }
    
    document.body.removeChild( textArea );
    
    if( successful )
    {
      console.log('Copied text to clipboard via appending textarea to DOM');
      Wt.emit( domSignalId, {name:'miscSignal'}, success_msg );
    }else
    {
      console.warn('Failed to copy to clipboard using textarea');
      Wt.emit( domSignalId, {name:'miscSignal'}, fail_msg );
    }
    
    return successful;
  }//if( !navigator.clipboard )
  
  navigator.clipboard.writeText(text).then(function() {
    console.log('Copying text to clipboard using async method.');
    Wt.emit( domSignalId, {name:'miscSignal'}, success_msg );
  }, function(err) {
    console.warn('Failed copying text to clipboard using async method.: ', err);
    Wt.emit( domSignalId, {name:'miscSignal'}, fail_msg );
  });
}
  );

}//namespace





namespace QrCode
{
pair<string,int> utf8_string_to_svg_qr( const std::string &input )
{
  // LOW (7% erroneous codewords), MEDIUM (15%), QUARTILE (25%), HIGH (30%).
  const qrcodegen::QrCode qr = qrcodegen::QrCode::encodeText( input.c_str(), qrcodegen::QrCode::Ecc::MEDIUM );
  
  cout << "Qr code created has version: " << qr.getVersion() << ", size: " << qr.getSize()
       << ", ErrorCorrectionLevel: " << static_cast<int>(qr.getErrorCorrectionLevel())
       << endl;
  
  return { to_svg_string( qr, 1 ), qr.getSize() };
}//utf8_string_to_svg_qr(...)

pair<string,int> binary_to_svg_qr( const std::vector<std::uint8_t> &data )
{
  const qrcodegen::QrCode qr = qrcodegen::QrCode::encodeBinary( data, qrcodegen::QrCode::Ecc::MEDIUM );
  
  return { to_svg_string( qr, 1 ), qr.getSize() };
}//binary_to_svg_qr(...)


SimpleDialog *displayTxtAsQrCode( const std::string &url,
                                 const std::string &title,
                                 const std::string &description )
{
  try
  {
    const pair<string,int> qr_svg = utf8_string_to_svg_qr( url );
    const string &qr_svg_str = qr_svg.first;
    const int qr_size = qr_svg.second;  //A simple DRF is like 70, or so
    
    int svg_size = 5*(qr_size+1);
    InterSpec *interspec = InterSpec::instance();
    if( interspec && (interspec->renderedWidth() > 100) && (interspec->renderedHeight() > 100) )
    {
      int wdim = std::min( interspec->renderedWidth(), interspec->renderedHeight() ) / 3;
      if( wdim > 2*qr_size )
        svg_size = std::min( wdim, 500 );
    }//if( InterSpec knows the window size )
    
    if( qr_svg_str.empty() )
    {
      passMessage( "Sorry, couldnt create QR code", 3 );
      return nullptr;
    }
    
    SimpleDialog *window = new SimpleDialog( title, "" );
    window->rejectWhenEscapePressed();
    window->addButton( "Close" );
    
    const unsigned char *svg_begin = (unsigned char *) &(qr_svg_str[0]);
    const unsigned char *svg_end = svg_begin + qr_svg_str.size();
    const vector<unsigned char> svg_data( svg_begin, svg_end );
    WMemoryResource *svgResource = new WMemoryResource( "image/svg+xml", svg_data, window );
    if( title.size() )
      svgResource->suggestFileName( title + ".svg", WResource::Attachment );
    else
      svgResource->suggestFileName( "qr.svg", WResource::Attachment );
    
    
    WImage *qrImage = new WImage( WLink(svgResource), window->contents() );
    //WText *qrImage = new WText( qr_svg.first, window->contents() );
    qrImage->setInline( false );
    qrImage->resize( svg_size, svg_size );
    qrImage->setMargin( 15, Wt::Side::Bottom );
    qrImage->setMargin( WLength::Auto, Wt::Side::Left );
    qrImage->setMargin( WLength::Auto, Wt::Side::Right );
    
    
    if( description.length() )
    {
      WText *message = new WText( description, window->contents() );
      message->addStyleClass( "content" );
      message->setInline( false );
    }
    
    WContainerWidget *btndiv = new WContainerWidget( window->contents() );
    
#if( BUILD_AS_OSX_APP || IOS )
    WAnchor *svgDownload = new WAnchor( WLink(svgResource), btndiv );
    svgDownload->setTarget( AnchorTarget::TargetNewWindow );
    svgDownload->setStyleClass( "LinkBtn DownloadLink DrfXmlDownload" );
#else
    WPushButton *svgDownload = new WPushButton( btndiv );
    svgDownload->setIcon( "InterSpec_resources/images/download_small.svg" );
    svgDownload->setLink( WLink(svgResource) );
    svgDownload->setLinkTarget( Wt::TargetNewWindow );
    svgDownload->setStyleClass( "LinkBtn DownloadBtn" );
    
#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    svgDownload->clicked().connect( std::bind([svgResource](){
      android_download_workaround(svgResource, "qr.svg");
    }) );
#endif //ANDROID
    
#endif
    svgDownload->setText( "SVG" );
    svgDownload->setFloatSide( Wt::Side::Right );
    
    LOAD_JAVASCRIPT(wApp, "QrCode.cpp", "QrCode", wtjsCopyUrlToClipboard );
    
    WPushButton *copyBtn = new WPushButton( "Copy url to clipboard", btndiv );
    copyBtn->setStyleClass( "LinkBtn" );
    copyBtn->setFloatSide( Wt::Side::Left );
    
    // TODO: there is probably a way to get wApp->root()->id() directly in the JS
    copyBtn->clicked().connect( "function(s,e){ "
                                 "Wt.WT.CopyUrlToClipboard(s,e,'" + copyBtn->id() + "','" + wApp->root()->id() + "','" + url + "');"
                                 "}" );
    
    
    return window;
  }catch( std::exception &e )
  {
    passMessage( "Sorry, error creating QR code: " + string(e.what()), 3 );
    cerr << "Error creating QR code window: " << e.what() << endl;
  }
  
  return nullptr;
}//displayTxtAsQrCode



}//namespace QrCode


