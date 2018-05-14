#include "InterSpec/ChartToImageResource.h"

#include <fstream>

#include <Wt/WLength>
#include <Wt/WPainter>
#include <Wt/WResource>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/Chart/WCartesianChart>

#if( MAKE_PNG_RESOURCE )
#include <Wt/WRasterImage>
#else
#include <Wt/WPdfImage>
#endif


using namespace Wt;

ChartToImageResource::ChartToImageResource( Chart::WCartesianChart *chart,
                                            WObject *parentParam )
  : WResource(parentParam), m_chart( chart ), m_width( 1600 ), m_height( 1200 )

{
}


ChartToImageResource::~ChartToImageResource()
{

}

const char *ChartToImageResource::imageType() const
{
#if( MAKE_PNG_RESOURCE )
  return "Png";
#else
  return "Pdf";
#endif
}//const char *imageType() const

void ChartToImageResource::handleRequest( const Http::Request& request,
                                        Http::Response& response)
{
#if( MAKE_PNG_RESOURCE )
  Wt::WRasterImage pngImage( "png", m_width, m_height );

  {
    Wt::WPainter p( &pngImage );
    m_chart->paint( p );
  }

  response.setMimeType( "image/png" );
  pngImage.write( response.out() );
#else

  WPdfImage pdfImage( width, height );

  {
    Wt::WPainter p( &pdfImage );
    m_chart->paint( p );
  }

  response.setMimeType( "application/pdf" );
  pdfImage.write( response.out() );
#endif

}//void handleRequest( const Http::Request &, Http::Response & )
