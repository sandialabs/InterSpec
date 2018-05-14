//  Created by wcjohns on 20110505

#ifndef ChartToPdfResource_h
#define ChartToPdfResource_h

#include "InterSpec_config.h"

#include <sstream>

#include <Wt/WResource>
#include <Wt/Http/Request>
#include <Wt/Http/Response>

//By Default ChartToImageResource will generate a PDF of the report, if however
//  You would like to generate a PNG (require WRasterImage), then set switch
//  below
#define MAKE_PNG_RESOURCE 1


//Some forward declarations
namespace Wt
{
  class WObject;
  namespace Chart
  {
    class WCartesianChart;
  }//namespace WCartesianChart
}//namespace Wt

class ChartToImageResource : public Wt::WResource
{
  //A simple class to turn a WStandardItemModel into a CSV data file.
  //The csv string data is not produced untill data is requiested
  //  (eg, only when the user clicks the WAnchor object that this WResouce
  //       is assigned to)
  //Note that this class would probably be replaced by just WRasterImage
  //  or WPdfImage
public:
  ChartToImageResource( Wt::Chart::WCartesianChart *chart,
                        Wt::WObject *parentParam = NULL );
  ~ChartToImageResource();

  const char *imageType() const;

private:
  Wt::Chart::WCartesianChart *m_chart;

  virtual void handleRequest( const Wt::Http::Request& request,
                              Wt::Http::Response& response );
  
  int m_width, m_height;
};//class ChartToImageResource









#endif
