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

#include <Wt/WText.h>
#include <Wt/WImage.h>
#include <Wt/WLabel.h>
#include <Wt/WAnchor.h>
#include <Wt/WComboBox.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/WLink.h>
#include <Wt/WMemoryResource.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WJavaScriptPreamble.h>

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
  function( sender, event, id, text )
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
      Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'}, success_msg );
    }else
    {
      console.warn('Failed to copy to clipboard using textarea');
      Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'}, fail_msg );
    }
    
    return successful;
  }//if( !navigator.clipboard )
  
  navigator.clipboard.writeText(text).then(function() {
    console.log('Copying text to clipboard using async method.');
    Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'}, success_msg );
  }, function(err) {
    console.warn('Failed copying text to clipboard using async method.: ', err);
    Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'}, fail_msg );
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
                                 const Wt::WString &title,
                                 const Wt::WString &description )
{
  try
  {
    InterSpec *interspec = InterSpec::instance();
    const bool is_phone = (interspec && interspec->isPhone());
    
    int w = interspec ? interspec->renderedWidth() : 0;
    int h = interspec ? interspec->renderedHeight() : 0;
    
    if( (interspec && interspec->isMobile()) && (w < 100) )
    {
      w = wApp->environment().screenWidth();
      h = wApp->environment().screenHeight();
    }
    
    const bool narrow_layout = ((w > 100) && (w < 480));
    
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
    if( (w > 100) && (h > 100) )
    {
      const int otherVertSpace = (is_phone ? ((h/20) + 95) : 205) + (title.empty() ? 0 : 15) + (description.empty() ? 0 : 20);
      
      int wdim = w / (is_phone ? 2 : 3);
      if( narrow_layout )
        wdim = (3*w) / 4;
      
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
    
    SimpleDialog *window = SimpleDialog::make( title, "" );
    window->rejectWhenEscapePressed();
    window->addButton( "Close" );
    
    if( is_phone )
      window->setAttributeValue( "style", "max-width: 95vw; " + window->attributeValue("style") );
    
    const unsigned char *svg_begin = (unsigned char *) &(qr_svg_str[0]);
    const unsigned char *svg_end = svg_begin + qr_svg_str.size();
    const vector<unsigned char> svg_data( svg_begin, svg_end );
    auto svgResource = std::make_shared<WMemoryResource>( "image/svg+xml", svg_data );
    if( !title.empty() )
      svgResource->suggestFileName( title + ".svg", ContentDisposition::Attachment );
    else
      svgResource->suggestFileName( "qr.svg", ContentDisposition::Attachment );
    
    string contentStyle;
    if( is_phone && !narrow_layout )
      contentStyle = "display: flex; "
                     "flex-direction: row; "
                     "column-gap: 10px; "
                     "max-width: calc(95vw - 30px); "
                     "max-height: calc(95vh - 65px); ";
    if( narrow_layout )
      contentStyle += "font-size: small; ";
    
    if( !contentStyle.empty() )
      window->contents()->setAttributeValue( "style", contentStyle );
    
    WImage *qrImage = window->contents()->addNew<WImage>( WLink(std::static_pointer_cast<WResource>(svgResource)) );
    qrImage->setInline( false );
    qrImage->resize( svg_size, svg_size );
    qrImage->setMargin( 5, Wt::Side::Bottom );
    if( title.empty() )
      qrImage->setMargin( 15, Wt::Side::Top );
    qrImage->setMargin( WLength::Auto, Wt::Side::Left );
    qrImage->setMargin( WLength::Auto, Wt::Side::Right );
    
    WContainerWidget *sizeRow = nullptr;
    if( narrow_layout )
    {
      sizeRow = window->contents()->addNew<WContainerWidget>();
      
      const char *desc_style = "margin-bottom: 10px;"
        " display: flex;"
        " flex-wrap: nowrap;"
        "justify-content: center;"
        " align-items: center;";
      sizeRow->setAttributeValue( "style", desc_style );
    }//if( narrow_layout )
    
    
    WContainerWidget *eclRow = window->contents()->addNew<WContainerWidget>();
    if( !sizeRow )
      sizeRow = eclRow;
    
    const char *ecl_style = (is_phone && !narrow_layout)
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
    if( is_phone && !narrow_layout )
    {
      tol_sel = eclRow->addNew<WContainerWidget>();
      tol_sel->setAttributeValue( "style",
                                 "display: flex;"
                                 " flex-direction: column;"
                                 " flex-wrap: nowrap;"
                                 " align-items: center;"
                                 " margin-top: 15px" );
    }
    
    eclRow->setAttributeValue( "style", ecl_style );
    WLabel *eclLabel = tol_sel->addNew<WLabel>( "Error Tolerance:" );
    WComboBox *eclSelect = tol_sel->addNew<WComboBox>();
    eclLabel->setBuddy( eclSelect );
    eclSelect->setNoSelectionEnabled( false );
    eclSelect->addItem( "Approx. 7% Loss" );
    eclSelect->addItem( "Approx. 15% Loss" );
    eclSelect->addItem( "Approx. 25% Loss" );
    eclSelect->addItem( "Approx. 30% Loss" );
    eclSelect->setCurrentIndex( static_cast<int>(ecl) );
    
    if( !is_phone || narrow_layout )
    {
      WText *spacer = eclRow->addNew<WText>( "&nbsp;" );
      spacer->setAttributeValue( "style", "flex: 1 1;" );
    }
    const string sizeDesc = to_string(qr_size) + "x" + to_string(qr_size) + " elements";
    assert( sizeRow );
    WText *sizeTxt = sizeRow->addNew<WText>( sizeDesc );
    
    eclSelect->changed().connect( eclSelect, [eclSelect, sizeTxt, url, svgResource](){
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
    } );
    
    
    if( !description.empty() )
    {
      auto messageOwner = std::make_unique<WText>( description );
      WText *message = messageOwner.get();
      if( is_phone && !narrow_layout )
      {
        eclRow->insertWidget( 0, std::move(messageOwner) );
      }else if( narrow_layout )
      {
        window->contents()->insertBefore( std::move(messageOwner), qrImage );
        message->addStyleClass( "content" );
        message->setInline( false );
        message->setAttributeValue( "style", "padding-bottom: 5px;" );
      }else
      {
        window->contents()->addWidget( std::move(messageOwner) );
        message->addStyleClass( "content" );
        message->setInline( false );
        // For a WText::setPadding() only supports left/right, so we will manually set bottom padding
        //  CSS, from the default 20px for .content, to 5px
        message->setAttributeValue( "style", "padding-bottom: 5px;" );
      }//if( is_phone ) / else
    }//if( description.length() )

    WContainerWidget *btndiv = ((is_phone && !narrow_layout) ? eclRow : window->contents())->addNew<WContainerWidget>();
    if( is_phone && !narrow_layout )
      btndiv->setAttributeValue( "style",
                                "display: flex; "
                                "flex-direction: column; "
                                "flex-wrap: nowrap; "
                                "align-items: center; "
                                "row-gap: 10px;" );

#if( BUILD_AS_OSX_APP || IOS )
    {
      WLink lnk( std::static_pointer_cast<WResource>(svgResource) );
      lnk.setTarget( LinkTarget::NewWindow );
      WAnchor *svgDownload = btndiv->addNew<WAnchor>( lnk );
      svgDownload->setStyleClass( "LinkBtn DownloadLink DrfXmlDownload" );
      svgDownload->setText( "SVG" );
      svgDownload->setFloatSide( Wt::Side::Right );
    }
#else
    WPushButton *svgDownload = btndiv->addNew<WPushButton>();
    svgDownload->setIcon( "InterSpec_resources/images/download_small.svg" );
    {
      WLink lnk( std::static_pointer_cast<WResource>(svgResource) );
      lnk.setTarget( LinkTarget::NewWindow );
      svgDownload->setLink( lnk );
    }
    svgDownload->setStyleClass( "LinkBtn DownloadBtn" );

#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    svgDownload->clicked().connect( svgDownload, [svgResource](){
      android_download_workaround(svgResource.get(), "qr.svg");
    } );
#endif //ANDROID

    svgDownload->setText( "SVG" );
    svgDownload->setFloatSide( Wt::Side::Right );
#endif

    LOAD_JAVASCRIPT(wApp, "QrCode.cpp", "QrCode", wtjsCopyUrlToClipboard );

    WPushButton *copyBtn = btndiv->addNew<WPushButton>( "Copy url to clipboard" );
    copyBtn->setStyleClass( "LinkBtn" );
    copyBtn->setFloatSide( Wt::Side::Left );

    // TODO: there is probably a way to get wApp->root()->id() directly in the JS
    const string escaped_url = WString(url).jsStringLiteral();
    copyBtn->clicked().connect( "function(s,e){ "
                                 "Wt.WT.CopyUrlToClipboard(s,e,'" + copyBtn->id() + "'," + escaped_url + ");"
                                 "}" );
    // svgDownload and copyBtn are already children of btndiv from addNew<> above,
    // so no additional addWidget calls are needed here.
    
    return window;
  }catch( std::exception &e )
  {
    passMessage( WString::tr("app-qr-err").arg(e.what()), 3 );
    cerr << "Error creating QR code window: " << e.what() << endl;
  }
  
  return nullptr;
}//displayTxtAsQrCode



}//namespace QrCode


