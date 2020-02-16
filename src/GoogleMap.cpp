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

// Disable streamsize <=> size_t warnings in boost
#pragma warning(disable:4244)

#include <vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <exception>

#include <Wt/WResource>
#include <Wt/WDateTime>
#include <Wt/WGoogleMap>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Http/Request>
#include <Wt/Http/Response>
#include <Wt/WContainerWidget>

#include "SpecUtils/SpecFile.h"
#include "InterSpec/GoogleMap.h"
#include "InterSpec/PhysicalUnits.h"

using namespace Wt;
using namespace std;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


struct PointInfo
{
  double latitude, longitude;
  Wt::WString html;
  
  PointInfo()
    : latitude(-999.9f),
      longitude(-999.9f),
      html()
  {}
  
  PointInfo( double lat, double lng, const WString &msg )
  : latitude(lat), longitude(lng), html(msg)
  {}
};//struct PointInfo


void wasClicked( const Wt::WGoogleMap::Coordinate &coord,
                 Wt::Signal<double /*Lat*/, double /*Lng*/> &signal )
{
  signal.emit( coord.latitude(), coord.longitude() );
}

class SrbGoogleMap: public Wt::WGoogleMap
{
public:
  SrbGoogleMap(  const bool trackMapExtent, WContainerWidget *parent = 0) :
    WGoogleMap( WGoogleMap::Version3, parent ),
    m_center( 37.675821, -121.709863 ),  //wcjohns office
    m_dontUseSandiaKey( false ),
    m_trackMapExtent( trackMapExtent )
  {
    if( m_trackMapExtent )
      m_mapGpsExtentFound.reset( new Wt::JSignal<double,double,double,double>( this, "mapGpsExtentRequested", true ) );
    
//    setMapTypeControl( WGoogleMap::DefaultControl );
  }
 
  void addPoint( const PointInfo &info )
  {
    m_points.push_back( info );
    const WGoogleMap::Coordinate coord(info.latitude,info.longitude);
    addMarker( coord );
    if( !info.html.empty() )
      openInfoWindow( coord, info.html );
  }//void addPoint( const PointInfo &info )
  
  
  void addInfoWindow( const PointInfo &info )
  {
    m_points.push_back( info );
    const WGoogleMap::Coordinate coord(info.latitude,info.longitude);
    addMarker( coord );
    
    if( info.html.empty() )
      return;
    
    string js;
    char lat[35], lng[35];
    snprintf(lat, sizeof(lat), "%.8f", info.latitude );
    snprintf(lng, sizeof(lng), "%.8f", info.longitude );
    
    js += "(function(){\n"
          "if( google===undefined || google.maps === undefined ) return;\n"
          "var infow = new google.maps.InfoWindow("
                                   "{content: '" + info.html.toUTF8() + "'});\n"
          "var latlng = new google.maps.LatLng(" + string(lat) + "," + string(lng) + ");\n"
          "var marker = new google.maps.Marker({position: latlng, map: " + jsRef() + ".map});\n"
          "google.maps.event.addListener(marker, 'click', function(){infow.open(map,marker);});\n"
          + jsRef() + ".map.overlays.push(marker);}\n)();";
    
    doGmJavaScript( js );
  }//void addInfoWindow( const PointInfo &info )

  
  
  //Right now the third argument, 'intensity' maps to a very crude color
  //  scale that should be improved in the future, as well as a legend added
  //  to the map widget, and even the ability for the user to customize the
  //  coloring.
  void addCircle( const double latitude, const double longitude,
                  const float intensity = -1.0f, const string &html = "" )
  {
    PointInfo info;
    info.longitude = longitude;
    info.latitude = latitude;
    info.html = html;
    m_points.push_back( info );
    
    
    //Should implement cubehelix, see: https://github.com/jradavenport/cubehelix/blob/master/cubehelix.py
    
    const WGoogleMap::Coordinate center(latitude,longitude);
    WColor strokeColor;
    if( intensity >= 0.0f && intensity <= 1.0f )
    {
      //generated at: http://www.perbang.dk/rgbgradient/
      //used FFFAFB to FB282D with 20 steps
      static const uint8_t reds[20]   = { 255, 254, 254, 254, 254, 253, 253, 253, 253, 253, 253, 252, 252, 252, 252, 251, 251, 251, 251, 251 };
      static const uint8_t greens[20] = { 250, 238, 227, 216, 205, 194, 183, 172, 161, 150, 139, 128, 117, 106,  95,  84,  73,  62,  51,  40 };
      static const uint8_t blues[20]  = { 251, 240, 229, 218, 207, 196, 185, 175, 164, 153, 142, 131, 120, 110,  99,  88,  77,  66,  55,  44 };
      static const size_t numGradients = sizeof(reds)/sizeof(reds[0]);
      
//      const int red = static_cast<int>( std::floor( 255*std::pow(intensity,3.0) + 0.5) );
//      const int greenblue = static_cast<int>( std::floor(255*(1.0-std::pow(intensity,1.0)) + 0.5) );
//      const int opacity = 75 + static_cast<int>(180.0*intensity);
//      strokeColor = WColor( red, greenblue, greenblue, opacity );
      
      const int opacity = 155 + static_cast<int>(100.0*intensity*intensity*intensity);
      const int index = std::min( static_cast<int>(numGradients-1), static_cast<int>( std::floor(numGradients*std::pow(intensity,4.0f)) ) );
      strokeColor = WColor( reds[index], greens[index], blues[index], opacity );
    }else
    {
      strokeColor = WColor(255,255,255,200);  //white
    }
    
    const WColor fillColor(255,255,255,0);
    
  
    std::stringstream strm;
    char lat[35], lng[35];
    snprintf(lat, sizeof(lat), "%.8f", latitude );
    snprintf(lng, sizeof(lng), "%.8f", longitude );
  
    const double radius = 15;
    const int strokeWidth = 6;
    
    const double strokeOpacity = strokeColor.alpha() / 255.0;
    const double fillOpacity = fillColor.alpha() / 255.0;
    
    strm << "(function(){"
            "if( google===undefined || google.maps === undefined ) return;\n"
            "var mapLocal = " << jsRef() + ".map;"
            "if(mapLocal===undefined) return;"
         << "var latLng = new google.maps.LatLng(" << lat << "," << lng << ");"
         << "var circle = new google.maps.Circle( "
            "{ "
            "  map: mapLocal, "
            "  radius: " << radius << ", "
            "  center: latLng,"
            "  position: latLng,"
            "  fillOpacity: \"" << fillOpacity << "\","
            "  fillColor: \"" << fillColor.cssText() << "\","
            "  strokeWeight: " << strokeWidth << ","
            "  strokeColor:\"" << strokeColor.cssText() << "\","
            "  strokeOpacity: " << strokeOpacity <<
            "} "
            ");";
    if( html.size() )
      strm << "var infow = new google.maps.InfoWindow({content: '" << html << "'});"
           << "google.maps.event.addListener(circle, 'click', function(){infow.open(map,circle);});";
    strm << "})();";
    
    
    doGmJavaScript(strm.str());

//    WGoogleMap::addCircle( center, 15, strokeColor, 6 , fillColor );
  }//void addCircle( double lat, double longitude )
  
  void clearPoints()
  {
    m_points.clear();
    clearOverlays();
  }//void clearPoints()
  
  
  void adjustPanAndZoom()
  {
    if( m_points.empty() )
      return;
    
    double minLat = 99999999.9, maxLat = -999999999.9;
    double minLon = 99999999.9, maxLon = -999999999.9;
    
    for( size_t i = 0; i < m_points.size(); ++i )
    {
      minLat = std::min( minLat, m_points[i].latitude );
      maxLat = std::max( maxLat, m_points[i].latitude );
      minLon = std::min( minLon, m_points[i].longitude );
      maxLon = std::max( maxLon, m_points[i].longitude );
    }//for( size_t i = 0; i < m_points.size(); ++i )
    
    if( m_points.size() > 1 )
    {
      const double lonDelta = 0.1*fabs(maxLon - minLon);
      const double latDelta = 0.1*fabs(maxLat - minLat);
      WGoogleMap::Coordinate ul( minLat-latDelta, maxLon+lonDelta );
      WGoogleMap::Coordinate lr( maxLat+latDelta, minLon-lonDelta );
      
      m_center.setLatitude( 0.5*(maxLat+minLat) );
      m_center.setLongitude( 0.5*(maxLon+minLon) );
      
      zoomWindow( ul, lr );
    }else
    {
      m_center.setLatitude( minLat );
      m_center.setLongitude( maxLon );
      
      panTo( WGoogleMap::Coordinate( minLat, maxLon ) );
      setZoom( 14 );
    }//if( m_points.size() > 1 )
  }//void adjustPanAndZoom()

  
protected:
  std::vector< PointInfo > m_points;
  WGoogleMap::Coordinate m_center;
  const bool m_dontUseSandiaKey;  //not tested for setting to true
  
  const bool m_trackMapExtent;
  std::unique_ptr<Wt::JSignal<double/*upper left latitude*/,double/*upper left longitude*/,
  double/*lower right latitude*/,double/*lower right longitude*/> >
      m_mapGpsExtentFound;
  friend class GoogleMap;
  
public:
  //We have to overide doGmJavaScript(...) so tht we can overide additions_
  //  with m_additions to allow us to execute this JS in render, since we
  //  are overiding it
  std::vector<std::string> m_additions;
  virtual void doGmJavaScript( const std::string& jscode )
  {
    if( m_dontUseSandiaKey )
    {
      WGoogleMap::doGmJavaScript( jscode );
      return;
    }//if( m_dontUseSandiaKey )
    
    if( isRendered() )
      doJavaScript( jscode );
    else
      m_additions.push_back( jscode );
  }//void doGmJavaScript( const std::string& jscode )
  
  
  void strmJSListener( const JSignal<Coordinate> &signal,
                                    std::string signalName,
                                    Wt::WStringStream &strm)
  {
      strm <<
      """if( google===undefined || google.maps === undefined ) return;"
      """google.maps.event.addListener(map, \"" << signalName << "\", "
      ""                              "function(event) {"
      ""  "if(event && event.latLng){"
      << signal.createCall("event.latLng.lat() +' '+ event.latLng.lng()")
      << ";}});";
  }//void strmJSListener(...)

  virtual void render( WFlags<RenderFlag> flags )
  {
    if( m_dontUseSandiaKey )
    {
      WGoogleMap::render( flags );
      return;
    }//const bool dontUseSandiaKey = true;

    
    if( flags & RenderFull )
    {
      WApplication *app = WApplication::instance();
      const string initFcn = app->javaScriptClass() + ".init_sandia_gmaps_" + id();
      
      // initialize the map
      WStringStream strm;
      strm <<"function() {"
      "var self = " << jsRef() << ";\n"
      "if(!self){setTimeout(" << initFcn << ", 0);return;}\n"
      "if( google===undefined || google.maps === undefined ) return;\n"
      "var myOptions = {"
        "" "zoom: 13,"
        "" "center: new google.maps.LatLng("
           << m_center.latitude() << ", " << m_center.longitude() << "),"
        "" "mapTypeId: google.maps.MapTypeId.ROADMAP"
      "};\n"
      "var map = new google.maps.Map(self, myOptions);\n"
      "map.overlays = [];\n"
      "map.infowindows = [];\n"
      "self.map = map;\n";
      
      if( m_trackMapExtent )
      {
        //overlay.getProjection() isnt defined until the map is drawn, so we
        //  must define it here, even if we wont ever use it... we could probably
        //  fix this by forcing a redraw, but whatever for now (since we;ll always
        //  have this capability enabled).
        strm << "var overlay = new google.maps.OverlayView();"
        "self.overlay = overlay;"
        "self.overlay.draw = function(){};"
        "self.overlay.setMap(map);"
        "var sendGpsExtent = function() {try{"
        "var ulpx = new google.maps.Point(0,0);"
        "var lrpx = new google.maps.Point($(self).width(),$(self).height());"
        "var ul = overlay.getProjection().fromContainerPixelToLatLng(ulpx);"
        "var lr = overlay.getProjection().fromContainerPixelToLatLng(lrpx);"
        "if(ul!==null&&lr!==null)"
          "Wt.emit( '" << id() << "', {name:'mapGpsExtentRequested'},"
          " ul.lat(), ul.lng(), lr.lat(), lr.lng() );"
        "}catch(e){console.log('Error in update map extent: ' + e );}};"
        "google.maps.event.addListener(map, 'idle', function(){sendGpsExtent();});";
      }//if( m_trackMapExtent )
      
      
      // eventhandling
      strmJSListener( clicked(), "click", strm );
      strmJSListener( doubleClicked(), "dblclick", strm );
      //      if (mouseMoved_)  //we'll never use this
      //        strmJSListener(*mouseMoved_, "mousemove", strm);
      
      // additional things
      for( const string &str : m_additions )
        strm << str;
      m_additions.clear();
    
      strm << "setTimeout(function(){delete " << initFcn << "; console.log('Deleted function!');},0);};";
      
      /*Add in a configuration option to wt_config to get the key; if there and non-empty
       Use it; else add in CMake option for key and if that is non empty, use it.
       Else, dont display the map.
       
       Also, add into the JavaScript:
       if ( google.maps == undefined ) { run the non-connected code }
      else { run the connected code }
    
       */
      
      //app->doJavaScript( strm.str(), true);  //started not working for some reason.
      app->declareJavaScriptFunction( "init_sandia_gmaps_" + id(), strm.str() );
      
      
      setJavaScriptMember(WT_RESIZE_JS,
                      "function(self, w, h){"
                      "if(w>=0)self.style.width=w+'px';"
                      "if(h>=0)self.style.height=h+'px';"
                      "if(self.map)google.maps.event.trigger(self.map,'resize');"
                      "}");
      
      
      string gmaps_key;
      
      //Check if a google maps key was specified at compile time.
#ifdef GOOGLE_MAPS_KEY
      if( GOOGLE_MAPS_KEY && strlen(GOOGLE_MAPS_KEY) > 1 )
        gmaps_key = GOOGLE_MAPS_KEY;
#endif
      
      //See if the configuration file contained a key we should use
      string configFileKey;
      if( Wt::WApplication::readConfigurationProperty("google_api_key", configFileKey) )
        gmaps_key = configFileKey;
      
      // if there is no google api key configured, use the one for
      // http://localhost:8080/
      if( gmaps_key.empty() )
        gmaps_key = "ABQIAAAAWqrN5o4-ISwj0Up_depYvhTwM0brOpm-All5BF6PoaKBxRWWERS-S9gPtCri-B6BZeXV8KpT4F80DQ";  //Taken from WGoogleMaps.C.. not sure if actually working.
      
        
      //Init the google javascript api
      // (will crash the client app if it doesnt load)
      const string jsapi = "https://maps.googleapis.com/maps/api/js?v=3&"
                           "key=" + gmaps_key + "&sensor=false&channel=InterSpec"
                           "&callback=" + initFcn;
      
      //const string jsapi = "https://maps.googleapis.com/maps/api/js?v=3&sensor=false&callback=" + initFcn;
      //If we use WApplication::require and we are not connected to the internet
      //  it will cause a fatal error and the user session will be ended.  So
      //  instead we will manually download the script, bypassing the caching
      //  mechanism of Wt, and instead use the browsers.
      //app->require(jsapi, "googlemapapi");
      
      string errorLoadingFcn = "function(jqxhr, textStatus, errorThrown){"
      "var self = " + jsRef() + ";\n"
      "console.log('loadin google maps api error thrown: ' + errorThrown + ', and text status: ' + textStatus);"
      "if(self)self.appendChild(document.createTextNode(\"Error loading google maps: \" + errorThrown + \". Check internet connection, or google maps key may be invalid\"));"
      "}";
    
      string successfcn = "function(data,textStatus,jqXHR){ $(window).data('HasLoadedGMaps',true); }";
      
      
      doJavaScript( "if(!$(window).data('HasLoadedGMaps'))"
                   " $.ajax({url: \"" + jsapi + "\","
                   + " dataType: \"script\""
                   + ", cache: true"
                   + ", timeout: 5000"
                   + ", error: " + errorLoadingFcn
                   + ", success: " + successfcn
                   + " });"
                   + " else " + initFcn + "();" );
      //doJavaScript( "$.getScript(\"" + jsapi + "\");" );
    }//if( flags & RenderFull )
    
    WCompositeWidget::render(flags);
  }//void render( WFlags<RenderFlag> flags )
};//class SrbGoogleMap



GoogleMap::GoogleMap( const bool trackMapExtent, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_map( 0 ),
    m_trackMapExtent( trackMapExtent ),
    m_clicked( this )
{
  init();
}
  
  
GoogleMap::~GoogleMap()
{
}


void GoogleMap::init()
{
  wApp->useStyleSheet( "InterSpec_resources/GoogleMap.css" );
  
  addStyleClass( "GoogleMapDiv" );
  m_map = new SrbGoogleMap( m_trackMapExtent );
  
  WGridLayout *layout = new WGridLayout();
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );
  layout->setContentsMargins( 0, 0, 0, 0 );
  setLayout( layout );
  
  layout->addWidget( m_map, 0, 0 );
  m_map->clicked().connect( boost::bind( &wasClicked, _1, boost::ref(m_clicked) ) );

  m_mapExtents[0] = m_mapExtents[1] = m_mapExtents[2] = m_mapExtents[3] = -999.9f;
  if( !!m_map->m_mapGpsExtentFound )
    m_map->m_mapGpsExtentFound->connect( boost::bind(&GoogleMap::updateMapGpsCoords, this, _1, _2, _3, _4 ) );
}//void init()

void GoogleMap::updateMapGpsCoords( const double upperLat, const double leftLng,
                                 const double lowerLat, const double rightLng )
{
  m_mapExtents[0] = static_cast<float>( upperLat );
  m_mapExtents[1] = static_cast<float>( leftLng );
  m_mapExtents[2] = static_cast<float>( lowerLat );
  m_mapExtents[3] = static_cast<float>( rightLng );
}


Wt::Signal<double /*Lat*/, double /*Lng*/> &GoogleMap::mapClicked()
{
  return m_clicked;
}


void GoogleMap::getMapExtent( float &upperLatitude, float &leftLongitude,
                              float &lowerLatitude, float &rightLongitude )
{
  if( !m_trackMapExtent )
    throw std::logic_error( "GoogleMap::getMapExtent(...): not initialized"
                            " with option to track map extent" );
  
  upperLatitude  = m_mapExtents[0];
  leftLongitude  = m_mapExtents[1];
  lowerLatitude  = m_mapExtents[2];
  rightLongitude = m_mapExtents[3];
}


void GoogleMap::addMarker( double latitude, double longitude )
{
  m_map->addPoint( PointInfo(latitude, longitude, "") );
}//void addMarker( double latitude, double longitude )


void GoogleMap::addInfoBox( double lat, double lng, const Wt::WString &html )
{
  m_map->addPoint( PointInfo(lat,lng,html) );
}//void addInfoBox(...)


void GoogleMap::adjustPanAndZoom()
{
  m_map->adjustPanAndZoom();
}

void GoogleMap::addMeasurment( std::shared_ptr<const SpecUtils::SpecFile> meas,
                               const WString &title,
                               const std::set<int> &displayed )
{
  if( !meas )
    return;

  typedef std::map< std::pair<float,float>, set<int> > GpsToMeas;
  GpsToMeas gpss;
  
  std::vector<float> cps;
  cps.reserve( meas->measurements().size() );
  
  for( const std::shared_ptr<const SpecUtils::Measurement> &m : meas->measurements() )
  {
    if( !m->has_gps_info() )
      continue;
    
    const pair<float,float> latlon( m->latitude(), m->longitude() );
    if( m->source_type() != SpecUtils::SourceType::IntrinsicActivity )
    {
      gpss[latlon].insert( m->sample_number() );
      cps.push_back( m->gamma_count_sum() / m->live_time() );
    }
  }//for( each measurments )
  
  std::sort( cps.begin(), cps.end() );
  
  char buffer[512];
  for( const GpsToMeas::value_type &vt : gpss )
  {
    const pair<float,float> &latlon = vt.first;
    const set<int> &samplenums = vt.second;
    
    string html = "<div clas=\"MapMeasPoint\">";
    
    const vector<std::shared_ptr<const SpecUtils::Measurement>> samples
                              = meas->sample_measurements( *vt.second.begin() );
    
    float lt = 0.0f, sumgammas = 0.0f;
    for( const std::shared_ptr<const SpecUtils::Measurement> &m : samples )
    {
      lt += m->live_time();
      sumgammas += m->gamma_count_sum();
    }
    
    html += "<div class=\"GpsTimeTxt\">";
    
    if( vt.second.size() < 2 )
    {
      // !samples[0]->position_time().is_special()
      if( samples.size() && !samples[0]->start_time().is_special() )
      {
        const WDateTime t = WDateTime::fromPosixTime( samples[0]->start_time());
        const string timeStr = t.toString(DATE_TIME_FORMAT_STR).toUTF8();
        html += timeStr + ", ";
      }//if( !m->start_time().is_special() )
    }else
    {
      snprintf( buffer, sizeof(buffer), "%i samples, ", int(vt.second.size()) );
      html += buffer;
    }//if( vt.second.size() < 2 ) / else
  
  
    const string ltstr = PhysicalUnits::printToBestTimeUnits( lt );
    html += "LT=" + ltstr + "</div>";
    
    
    if( samplenums.size() == 1 )
    {
      snprintf( buffer, sizeof(buffer), "<div>Sample %i</div>", *samplenums.begin() );
      html += buffer;
    }else if( samplenums.size() > 1 && samplenums.size() < 7 )
    {
      html += "<div>Sample Nums: ";
      for( set<int>::const_iterator i = samplenums.begin(); i != samplenums.end(); ++i )
      {
        if( i != samplenums.begin() )
          html += ", ";
        html += std::to_string(*i);
      }
      html += "</div>";
    }//if( samplenums.size() == 1 ) / else if there is more samples
    
    
    snprintf( buffer, sizeof(buffer),
              "<div class=\"GpsCoordsTxt\">(%.6f, %.6f)</div>",
              latlon.first, latlon.second );
    html += buffer;
    html += "</div>";
    
    
    //Lets make sure we dont load to much information to the client incase there
    //  are thousands of gps points.
    if( (gpss.size() < 10 || (gpss.size() < 20 && !meas->passthrough())) )
    {
      bool isdisplayed = false;
      for( const int s : samplenums )
        isdisplayed = (isdisplayed || displayed.count(s));
      
      if( isdisplayed && displayed.size() < 10 )
        m_map->addPoint( PointInfo(latlon.first, latlon.second, html) );
      else
        m_map->addInfoWindow( PointInfo(latlon.first, latlon.second, html) );
    }else
    {
      const std::vector<float>::const_iterator lb = std::upper_bound(cps.begin(),
                                                    cps.end(), (sumgammas/lt) );
      const float intensity = static_cast<float>(lb-cps.begin())/cps.size();
      m_map->addCircle( latlon.first, latlon.second, intensity, html );
    }
  }//for( const GpsToMeas::value_type &vt : gpss )
  
  m_map->adjustPanAndZoom();
}//void addMeasurment(...)


void GoogleMap::clearMeasurments()
{
  m_map->clearPoints();
}//void clearMeasurments()

