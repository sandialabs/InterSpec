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
#include <iostream>
#include <algorithm>

#include <Wt/WText>
#include <Wt/WCheckBox>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include <Wt/Json/Array>
#include <Wt/Json/Value>
#include <Wt/Json/Object>
#include <Wt/Json/Parser>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/ParseUtils.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/LeafletRadMap.h"

#if( IOS )
#include "InterSpec/InterSpecApp.h"
#endif

using namespace Wt;
using namespace std;

namespace
{
  void showMapWindow( const std::shared_ptr<const SpecMeas> meas,
                     double latitude, double longitude,
                     std::function<void(LeafletRadMap *)> on_create )
  {
    LeafletRadMapWindow *window = new LeafletRadMapWindow();
    
    LeafletRadMap *gpsmap = window->map();
    if( meas )
      gpsmap->displayMeasurementOnMap( meas );
    else if( SpecUtils::valid_latitude(latitude) && SpecUtils::valid_longitude(longitude) )
      gpsmap->displayCoordinate( latitude, longitude );
    
    if( on_create )
      on_create( gpsmap );
  }//void showMapWindow(...)
  
  
  SimpleDialog *startShowMap( const std::shared_ptr<const SpecMeas> meas,
                              double latitude, double longitude,
                              std::function<void(LeafletRadMap *)> on_create )
  {
    InterSpec *viewer = InterSpec::instance();
    
    const bool showWarning = InterSpecUser::preferenceValue<bool>( "ShowMapDataWarning", viewer );
    
    if( !showWarning )
    {
      showMapWindow( meas, latitude, longitude, on_create );
      return nullptr;
    }//if( !showWarning )
    
    // Show a warning dialog about requesting the data, before proceeding
    const char *title = "Before Proceeding";
    const char *msg =
    "Map tiles will be requested from <a href=\"https://arcgis.com\">https://arcgis.com</a>.<br />"
    "No radiation data will leave your device, but requests for<br />"
    "map tiles encompassing the measurements will be made to this service.";
    SimpleDialog *dialog = new SimpleDialog( title );
    
    WText *message = new WText( msg, dialog->contents() );
    message->addStyleClass( "content" );
    message->setInline( false );
    
    WCheckBox *cb = new WCheckBox( "Dont ask again", dialog->contents() );
    cb->setInline( false );
    cb->checked().connect( std::bind([cb](){
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        InterSpecUser::setPreferenceValue(viewer->m_user, "ShowMapDataWarning", !cb->isChecked(), viewer);
    }) );
    
    WPushButton *accept = dialog->addButton( "Proceed" );
    accept->clicked().connect( boost::bind( &showMapWindow, meas, latitude, longitude, on_create) );
    WPushButton *cancel = dialog->addButton( "Cancel" );
    cancel->clicked().connect( std::bind([](){
      InterSpec *viewer = InterSpec::instance();
      if( viewer )
        InterSpecUser::setPreferenceValue(viewer->m_user, "ShowMapDataWarning", true, viewer);
    }) );
    
    return dialog;
  }//static LeafletRadMap::showForMeasurement
}//namespace


SimpleDialog *LeafletRadMap::showForMeasurement( const std::shared_ptr<const SpecMeas> meas,
                                  std::function<void(LeafletRadMap *)> on_create )
{
  return startShowMap( meas, -999.9, -999.9, on_create );
}//


SimpleDialog *LeafletRadMap::showForCoordinate( double latitude, double longitude,
                                 std::function<void(LeafletRadMap *)> on_create )
{
  return startShowMap( nullptr, latitude, longitude, on_create );
}//static SimpleDialog *showForCoordinate


LeafletRadMapWindow::LeafletRadMapWindow()
: AuxWindow( "Map Tool",
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::EnableResize) ),
  m_map( nullptr )
{
  InterSpec *viewer = InterSpec::instance();
  
  m_map = new LeafletRadMap( contents() );
  m_map->setHeight( WLength(100,WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "leaflet-rad-map" );
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  resizeScaledWindow( 0.8, 0.8 );
  
  show();
  
  if( viewer && (viewer->renderedHeight() > 100) )
  {
    float safeAreas[4] = { 0.0f };
    
#if( IOS )
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    if( app )
    {
      InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
      app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
    }//if( app )
#endif
    setMaximumSize( WLength::Auto, viewer->renderedHeight() - std::max(0.5f*(safeAreas[0]+safeAreas[2]),6.0f) );
  }//if( we know the screens size )
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//LeafletRadMapWindow constructor
  
  
LeafletRadMapWindow::~LeafletRadMapWindow()
{
}


LeafletRadMap *LeafletRadMapWindow::map()
{
  return m_map;
}


LeafletRadMap::LeafletRadMap( Wt::WContainerWidget *parent )
 : Wt::WContainerWidget( parent ),
  m_meas(),
  m_jsmap( jsRef() + ".map" ),
  m_displaySamples( this, "loadSamples", false),
  m_loadSelected( this )
{
  addStyleClass( "LeafletRadMap" );
      
  // Load Leaflet 1.9.3 (20221118), see https://leafletjs.com/index.html
  wApp->useStyleSheet( "InterSpec_resources/assets/js/leaflet/leaflet_1.9.3/leaflet.css" );
  wApp->require( "InterSpec_resources/assets/js/leaflet/leaflet_1.9.3/leaflet.js" );

  // Load Esri Leaflet 3.1.10 (20230117), see https://github.com/Esri/esri-leaflet
  wApp->require( "InterSpec_resources/assets/js/leaflet/esri-leaflet_3.1.10/esri-leaflet.js" );
  
  // Load Esri Leaflet Vector Tile Plugin 4.0.0 (20220902), see https://github.com/Esri/esri-leaflet-vector
  wApp->require( "InterSpec_resources/assets/js/leaflet/esri-leaflet-vector_4.0.0/esri-leaflet-vector.js" );

  // Load Leaflet-Geoman 2.14.1 (20230118), see https://github.com/geoman-io/leaflet-geoman
  wApp->useStyleSheet( "InterSpec_resources/assets/js/leaflet/leaflet-geoman-free_2.14.1/leaflet-geoman.css" );
  wApp->require( "InterSpec_resources/assets/js/leaflet/leaflet-geoman-free_2.14.1/leaflet-geoman.min.js" );
      
  // Include d3-polygon 3.0.1 (20210605), see https://github.com/d3/d3-polygon
  wApp->require( "InterSpec_resources/assets/js/leaflet/d3-polygon/d3-polygon.min.js" );
    
  // Include Leaflet.heat 0.2.0 (20151026), see https://github.com/Leaflet/Leaflet.heat
  wApp->require( "InterSpec_resources/assets/js/leaflet/leaflet-heat_0.2.0/leaflet-heat.js" );

  // Include leaflet-markercluster 1.5.3 (20211017), see https://github.com/Leaflet/Leaflet.markercluster
  wApp->require( "InterSpec_resources/assets/js/leaflet/leaflet-markercluster_1.5.3/leaflet.markercluster-src.js" );
      
  // Include Leaflet.EasyButton 2.4.0 (20190110), see https://github.com/CliffCloud/Leaflet.EasyButton
  wApp->useStyleSheet( "InterSpec_resources/assets/js/leaflet/leaflet-EasyButton_2.4.0/easy-button.css" );
  wApp->require( "InterSpec_resources/assets/js/leaflet/leaflet-EasyButton_2.4.0/easy-button.js" );
  
  // Finally include our JS/CSS
  wApp->useStyleSheet( "InterSpec_resources/LeafletRadMap.css" );
  wApp->require( "InterSpec_resources/LeafletRadMap.js" );

  
  m_displaySamples.connect( boost::bind( &LeafletRadMap::handleLoadSamples, this,
                                        boost::placeholders::_1, boost::placeholders::_2 ) );
}//LeafletRadMap


LeafletRadMap::~LeafletRadMap()
{
}//~LeafletRadMap()
  

void LeafletRadMap::defineJavaScript()
{
  string options = "{"
    "apiKey: '" LEAFLET_MAPS_KEY "'"
  "}";
  
  setJavaScriptMember( "map", "new LeafletRadMap(" + jsRef() + "," + options + ");");
  
  //setJavaScriptMember( "resizeObserver",
  //  "new ResizeObserver(entries => {"
  //    "for (let entry of entries) {"
  //      "if( entry.target && (entry.target.id === '" + id() + "') )"
  //        + m_jsmap + ".handleResize();"
  //    "}"
  //  "});"
  //);
  //callJavaScriptMember( "resizeObserver.observe", m_chart->jsRef() );
  
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void LeafletRadMap::doJavaScript( const std::string& js )
{
  if( isRendered() )
    WContainerWidget::doJavaScript( js );
  else
    m_pendingJs.push_back( std::move(js) );
}//doJavaScript(...)


void LeafletRadMap::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


std::string LeafletRadMap::createGeoLocationJson( std::shared_ptr<const SpecMeas> meas )
{
  if( !meas || !meas->has_gps_info() )
  {
    cout << "createGeoLocationJson: no measurement, or no GPS info" << endl;
    return "null";
  }
  
  Wt::Json::Object json;//( {{"userLabel", Wt::Json::Value( Wt::WString::fromUTF8(label) )} } );
  json["fileName"] = Wt::WString::fromUTF8( meas->filename() );
  
  const vector<string> &det_names = meas->detector_names();
  Wt::Json::Array detectors;
  for( const auto det : det_names )
    detectors.push_back( Wt::WString::fromUTF8(det) );
  json["detectors"] = detectors;
  
  Wt::Json::Array samplesData;
  
  size_t numSamplesNoGps = 0;
  SpecUtils::time_point_t start_times_offset{};
  
  for( const int sample : meas->sample_numbers() )
  {
    vector<shared_ptr<const SpecUtils::Measurement>> meass = meas->sample_measurements(sample);
    
    meass.erase( std::remove_if(meass.begin(), meass.end(),
      [&](shared_ptr<const SpecUtils::Measurement> a){
        return !a || (end(det_names) == std::find(begin(det_names), end(det_names), a->detector_name()));
    }), meass.end() );
    
    
    bool hadGps = false;
    double latitude = -999.99, longitude = -999.99;
    int numGammaDet = 0, numNeutDet = 0;
    float realTime = 0.0f;
    SpecUtils::SourceType source_type = SpecUtils::SourceType::Unknown;
    double gammaRealTime = 0.0, gammaLiveTime = 0.0, gammaCounts = 0.0;
    double neutronRealTime = 0.0, neutronCounts = 0.0;
    
    SpecUtils::time_point_t meas_start_time{};
    
    for( const shared_ptr<const SpecUtils::Measurement> &m : meass )
    {
      if( m->has_gps_info() )
      {
        hadGps = true;
        // TODO: we could average or something here... not sure it matters much
        latitude = m->latitude();
        longitude = m->longitude();
      }
      
      if( SpecUtils::is_special(meas_start_time) && !SpecUtils::is_special(m->start_time()) )
      {
        meas_start_time = m->start_time();
        if( SpecUtils::is_special(start_times_offset) )
          start_times_offset = meas_start_time;
      }
      
      if( m->num_gamma_channels() >= 1 )
      {
        numGammaDet += 1;
        realTime = std::max( realTime, m->real_time() );
        gammaLiveTime += m->live_time();
        gammaRealTime += m->real_time();
        gammaCounts += m->gamma_count_sum();
      }
      
      if( m->contained_neutron() )
      {
        numNeutDet += 1;
        neutronRealTime += m->real_time();
        neutronCounts += m->neutron_counts_sum();
      }
      
      realTime = std::max( realTime, m->real_time() );
      
      if( m->source_type() != SpecUtils::SourceType::Unknown )
        source_type = m->source_type();
    }//for( const auto m : meass )
    
    if( !hadGps )
    {
      numSamplesNoGps += 1;
      continue;
    }
    
    if( !numNeutDet && !numGammaDet )
      continue;
    
    Wt::Json::Object sample_json;
    if( !SpecUtils::is_special(meas_start_time) )
    {
      const auto duration = start_times_offset - meas_start_time;
      const auto millisecs = chrono::duration_cast<chrono::milliseconds>(duration);
      sample_json["timeOffset"] = millisecs.count();
    }//if( SpecUtils::is_special(meas_start_time) )
    
    sample_json["nGDet"] = numGammaDet;
    if( numGammaDet )
    {
      sample_json["gSum"] = gammaCounts;
      sample_json["gRT"] = gammaRealTime;
      sample_json["gLT"] = gammaLiveTime;
    }
    
    sample_json["nNDet"] = numNeutDet;
    if( numNeutDet )
    {
      sample_json["nSum"] = neutronCounts;
      sample_json["nRT"] = neutronRealTime;
    }
    
    sample_json["rt"] = realTime;
    sample_json["src"] = static_cast<int>(source_type);
    
    Wt::Json::Array gps;
    gps.push_back(latitude);
    gps.push_back(longitude);
    sample_json["gps"] = gps;
    sample_json["sample"] = sample;
    
    samplesData.push_back( sample_json );
  }//for( const int sample : meas->sample_numbers() )
  
  if( !SpecUtils::is_special(start_times_offset) )
  {
    const auto dur = start_times_offset.time_since_epoch();
    const auto millisecs = chrono::duration_cast<chrono::milliseconds>(dur);
    json["startTime"] = millisecs.count();
  }//if( !SpecUtils::is_special(start_times_offset) )
  
  json["samples"] = samplesData;
  
  return Wt::Json::serialize(json);
}//void printGeoLocationJson( std::shared_ptr<SpecMeas> meas )



void LeafletRadMap::handleLoadSamples( const std::string &samples, std::string meas_type )
{
  cout << "Handle load samples: '" << samples << "', meas_type=" << meas_type << endl;
}//void handleLoadSamples( const std::vector<int> &samples, std::string meas_type );


void LeafletRadMap::displayMeasurementOnMap( const std::shared_ptr<const SpecMeas> meas )
{
  
  const string json = createGeoLocationJson( meas );
  
  cout << "displayMeasurementOnMap --> " << json << endl;
  
  doJavaScript( m_jsmap +  ".setData( " + json + " );" );
}//void displayMeasurementOnMap(...)


void LeafletRadMap::displayCoordinate( double latitude, double longitude )
{
  // blah blah blah
}//void displayCoordinate( double latitude, double longitude );



