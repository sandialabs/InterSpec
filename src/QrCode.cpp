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

#include <cctype>
#include <string>
#include <sstream>
#include <iostream>

#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WComboBox>
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

  
bool starts_with( const std::string &line, const char *label )
{
  const size_t len1 = line.size();
  const size_t len2 = strlen(label);
  
  if( len1 < len2 )
    return false;
  
  for( size_t i = 0; i < len2; ++i )
  {
    if( std::toupper( (int)line[i] ) != std::toupper( (int)label[i] ) )
      return false;
  }
  
  return true;
}//bool starts_with( const std::string &line, const char *label )
}//namespace





namespace QrCode
{
/** Converts a QR code to a SVG string we can include in the DOM.
 
 Adapted from QrCodeGeneratorDemo.cpp that is part of QR Code generator library.
 */
std::string to_svg_string( const qrcodegen::QrCode &qr, int border )
{
  if( border < 0 )
    throw std::domain_error("Border must be non-negative");
  
  if( border > 256 )
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


tuple<string,int,ErrorCorrLevel> utf8_string_to_svg_qr( const std::string &input,
                                                        ErrorCorrLevel prefferedECL,
                                                       const int quietSpace )
{
  auto toOurEnum = []( qrcodegen::QrCode::Ecc ecc ) -> ErrorCorrLevel {
    switch( ecc )
    {
      case qrcodegen::QrCode::Ecc::LOW:
        return ErrorCorrLevel::About7Percent;

      case qrcodegen::QrCode::Ecc::MEDIUM:
        return ErrorCorrLevel::About15Percent;
        
      case qrcodegen::QrCode::Ecc::QUARTILE:
        return ErrorCorrLevel::About25Percent;
        
      case qrcodegen::QrCode::Ecc::HIGH:
        return ErrorCorrLevel::About30Percent;
    }//switch( ecc )
    
    assert( 0 );
    return ErrorCorrLevel::About30Percent;
  };//auto toOurEnum
  
  vector<qrcodegen::QrCode::Ecc> eccs{
    qrcodegen::QrCode::Ecc::HIGH,
    qrcodegen::QrCode::Ecc::QUARTILE,
    qrcodegen::QrCode::Ecc::MEDIUM,
    qrcodegen::QrCode::Ecc::LOW
  };
  
  switch( prefferedECL )
  {
    case ErrorCorrLevel::About7Percent:
      eccs = { qrcodegen::QrCode::Ecc::LOW };
      break;
      
    case ErrorCorrLevel::About15Percent:
      eccs = { qrcodegen::QrCode::Ecc::MEDIUM, qrcodegen::QrCode::Ecc::LOW };
      break;
      
    case ErrorCorrLevel::About25Percent:
      eccs = { qrcodegen::QrCode::Ecc::QUARTILE, qrcodegen::QrCode::Ecc::MEDIUM, qrcodegen::QrCode::Ecc::LOW };
      break;
      
    case ErrorCorrLevel::About30Percent:
      break;
  }//switch( prefferedECL )
  
  
  for( const qrcodegen::QrCode::Ecc ecc : eccs )
  {
    try
    {
      // LOW (7% erroneous codewords), MEDIUM (15%), QUARTILE (25%), HIGH (30%).
      qrcodegen::QrCode qr = qrcodegen::QrCode::encodeText( input.c_str(), ecc );
      
      //cout << "Qr code created has version: " << qr.getVersion() << ", size: " << qr.getSize()
      //  << ", ErrorCorrectionLevel: " << static_cast<int>(qr.getErrorCorrectionLevel())
      //  << endl;
      
      return { to_svg_string( qr, quietSpace ), qr.getSize(), toOurEnum(ecc) };
    }catch( std::exception &e )
    {
      cerr << "Failed to encode QR code at level " << static_cast<int>(ecc) << endl;
    }
  }//for( const qrcodegen::QrCode::Ecc ecc : eccs )
  
  throw runtime_error( "Failed to be able to encode URL to QR code" );
  
  return { "", 0, ErrorCorrLevel::About30Percent };
}//utf8_string_to_svg_qr(...)

  
pair<string,int> binary_to_svg_qr( const std::vector<std::uint8_t> &data )
{
  const qrcodegen::QrCode qr = qrcodegen::QrCode::encodeBinary( data, qrcodegen::QrCode::Ecc::MEDIUM );
  
  return { to_svg_string( qr, 5 ), qr.getSize() };
}//binary_to_svg_qr(...)


SimpleDialog *displayTxtAsQrCode( const std::string &url,
                                 const std::string &title,
                                 const std::string &description )
{
  try
  {
    InterSpec *interspec = InterSpec::instance();
    const bool is_phone = (interspec && interspec->isPhone());
    
    // If its a "mailto:" URI, then put error correction to lowest level, to make the QR code
    //  as small as possible, since user will most likely be directly reading the QR code onto
    //  their phone.
    // For phones with a small screen, prefer fewer elements
    QrCode::ErrorCorrLevel wanted_ecl = ErrorCorrLevel::About30Percent;
    if( starts_with( url, "mailto:") || is_phone )
      wanted_ecl = ErrorCorrLevel::About7Percent;
    
    const tuple<string,int,ErrorCorrLevel> qr_svg = utf8_string_to_svg_qr( url, wanted_ecl, 3 );
    const string &qr_svg_str = get<0>(qr_svg);
    const int qr_size = get<1>(qr_svg);  //A simple DRF is like 70, or so
    const ErrorCorrLevel ecl = get<2>(qr_svg);
    
    int svg_size = 5*(qr_size+1);
    if( interspec && (interspec->renderedWidth() > 100) && (interspec->renderedHeight() > 100) )
    {
      const int h = interspec->renderedHeight();
      
      const int otherVertSpace = (is_phone ? ((h/20) + 95) : 205) + (title.empty() ? 0 : 15) + (description.empty() ? 0 : 20);
      
      int wdim = interspec->renderedWidth() / (is_phone ? 2 : 3);
      wdim = std::min( wdim, (h - otherVertSpace) );
      wdim = std::max( wdim, 125 );
      svg_size = std::min( wdim, svg_size );
    }//if( InterSpec knows the window size )
    
    svg_size = std::min( svg_size, 640 );
    
    if( qr_svg_str.empty() )
    {
      passMessage( "Sorry, couldn't create QR code", 3 );
      return nullptr;
    }
    
    SimpleDialog *window = new SimpleDialog( title, "" );
    window->rejectWhenEscapePressed();
    window->addButton( "Close" );
    
    if( is_phone )
      window->setAttributeValue( "style", "max-width: 95vw; " + window->attributeValue("style") );
    
    const unsigned char *svg_begin = (unsigned char *) &(qr_svg_str[0]);
    const unsigned char *svg_end = svg_begin + qr_svg_str.size();
    const vector<unsigned char> svg_data( svg_begin, svg_end );
    WMemoryResource *svgResource = new WMemoryResource( "image/svg+xml", svg_data, window );
    if( title.size() )
      svgResource->suggestFileName( title + ".svg", WResource::Attachment );
    else
      svgResource->suggestFileName( "qr.svg", WResource::Attachment );
    
    if( is_phone )
      window->contents()->setAttributeValue( "style",
                                            "display: flex; "
                                            "flex-direction: row; "
                                            "column-gap: 10px; "
                                            "max-width: calc(95vw - 30px); "
                                            "max-height: calc(95vh - 65px); ");
    
    WImage *qrImage = new WImage( WLink(svgResource), window->contents() );
    qrImage->setInline( false );
    qrImage->resize( svg_size, svg_size );
    qrImage->setMargin( 5, Wt::Side::Bottom );
    if( title.empty() )
      qrImage->setMargin( 15, Wt::Side::Top );
    qrImage->setMargin( WLength::Auto, Wt::Side::Left );
    qrImage->setMargin( WLength::Auto, Wt::Side::Right );
    
    WContainerWidget *eclRow = new WContainerWidget( window->contents() );
    const char *ecl_style = is_phone
      ? "margin: 10px;"
        " display: flex;"
        " flex-direction: column;"
        " align-items: center;"
        " flex-wrap: nowrap;"
        " justify-content: space-between;"
      : "margin-bottom: 10px;"
        " display: flex;"
        " flex-wrap: nowrap;"
        " align-items: center;"
        " gap: 3px;"
        " margin-left: 10px;"
        " margin-right: 10px;";
    
    
    WContainerWidget *tol_sel = eclRow;
    if( is_phone )
    {
      tol_sel = new WContainerWidget( eclRow );
      tol_sel->setAttributeValue( "style",
                                 "display: flex;"
                                 " flex-direction: column;"
                                 " flex-wrap: nowrap;"
                                 " align-items: center;"
                                 " margin-top: 15px" );
    }
    
    eclRow->setAttributeValue( "style", ecl_style );
    WLabel *eclLabel = new WLabel( "Error Tolerance:", tol_sel );
    WComboBox *eclSelect = new WComboBox( tol_sel );
    eclLabel->setBuddy( eclSelect );
    eclSelect->setNoSelectionEnabled( false );
    eclSelect->addItem( "Approx. 7% Loss" );
    eclSelect->addItem( "Approx. 15% Loss" );
    eclSelect->addItem( "Approx. 25% Loss" );
    eclSelect->addItem( "Approx. 30% Loss" );
    eclSelect->setCurrentIndex( static_cast<int>(ecl) );
    
    if( !is_phone )
    {
      WText *spacer = new WText( "&nbsp;", eclRow );
      spacer->setAttributeValue( "style", "flex: 1 1;" );
    }
    const string sizeDesc = to_string(qr_size) + "x" + to_string(qr_size) + " elements";
    WText *sizeTxt = new WText( sizeDesc, eclRow );
    
    eclSelect->changed().connect( std::bind([eclSelect, sizeTxt, url, svgResource](){
      const int ecl = eclSelect->currentIndex();
      assert( ecl >= 0 && ecl <= 3 );
      if( ecl < 0 || ecl > 4 )
        return;
      
      try
      {
        const tuple<string,int,ErrorCorrLevel> qr_svg
                                = utf8_string_to_svg_qr( url, static_cast<ErrorCorrLevel>(ecl), 3 );
        
        const string &qr_svg_str = get<0>(qr_svg);
        if( qr_svg_str.empty() )
          throw runtime_error( "Error creating SVG" );
        
        const unsigned char *svg_begin = (unsigned char *) &(qr_svg_str[0]);
        const unsigned char *svg_end = svg_begin + qr_svg_str.size();
        const vector<unsigned char> svg_data( svg_begin, svg_end );
        
        svgResource->setData( svg_data );
        
        const int actualEcl = static_cast<int>( get<2>(qr_svg) );
        assert( actualEcl >= 0 && actualEcl <= 3 );
        eclSelect->setCurrentIndex( actualEcl );
        
        const int qr_size = get<1>(qr_svg);
        const string sizeDesc = to_string(qr_size) + "x" + to_string(qr_size) + " elements";
        sizeTxt->setText( sizeDesc );
      }catch( std::exception & )
      {
        // This shouldnt happen
        eclSelect->setNoSelectionEnabled( true );
        eclSelect->setCurrentIndex( -1 );
        passMessage( "Sorry, error setting tolerance level.", 3 );
      }//try / catch
    }) );
    
    
    if( description.length() )
    {
      WText *message = new WText( description );
      if( is_phone )
      {
        eclRow->insertWidget( 0, message );
      }else
      {
        window->contents()->addWidget( message );
        message->addStyleClass( "content" );
        message->setInline( false );
        // For a WText::setPadding() only supports left/right, so we will manually set bottom padding
        //  CSS, from the default 20px for .content, to 5px
        message->setAttributeValue( "style", "padding-bottom: 5px;" );
      }//if( is_phone ) / else
    }//if( description.length() )
    
    WContainerWidget *btndiv = new WContainerWidget( is_phone ? eclRow : window->contents() );
    if( is_phone )
      btndiv->setAttributeValue( "style",
                                "display: flex; "
                                "flex-direction: column; "
                                "flex-wrap: nowrap; "
                                "align-items: center; "
                                "row-gap: 10px;" );
    
#if( BUILD_AS_OSX_APP || IOS )
    WAnchor *svgDownload = new WAnchor( WLink(svgResource) );
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
    const string escaped_url = WString(url).jsStringLiteral();
    copyBtn->clicked().connect( "function(s,e){ "
                                 "Wt.WT.CopyUrlToClipboard(s,e,'" + copyBtn->id() + "','" + wApp->root()->id() + "'," + escaped_url + ");"
                                 "}" );
    if( is_phone )
    {
      btndiv->addWidget( svgDownload );
      btndiv->addWidget( copyBtn );
    }else
    {
      btndiv->addWidget( copyBtn );
      btndiv->addWidget( svgDownload );
    }
    
    return window;
  }catch( std::exception &e )
  {
    passMessage( "Sorry, error creating QR code: " + string(e.what()), 3 );
    cerr << "Error creating QR code window: " << e.what() << endl;
  }
  
  return nullptr;
}//displayTxtAsQrCode



}//namespace QrCode


