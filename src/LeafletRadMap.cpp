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
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/LeafletRadMap.h"
#include "InterSpec/WarningWidget.h"

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
#include <Wt/Utils>
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecServer.h"
#endif


#if( IOS )
#include "InterSpec/InterSpecApp.h"
#endif

using namespace Wt;
using namespace std;

namespace
{
  void showMapWindow( const std::shared_ptr<const SpecMeas> meas,
                     const set<int> &sample_numbers,
                     const vector<string> &detector_names,
                     std::function<void(LeafletRadMapWindow *)> on_create )
  {
    LeafletRadMapWindow *window = new LeafletRadMapWindow();
    
    LeafletRadMap *gpsmap = window->map();
    gpsmap->displayMeasurementOnMap( meas, sample_numbers, detector_names );
    
    if( on_create )
      on_create( window );
  }//void showMapWindow(...)
}//namespace


SimpleDialog *LeafletRadMap::showForMeasurement( const std::shared_ptr<const SpecMeas> meas,
                                                 const set<int> &sample_numbers,
                                                 const vector<string> &detector_names,
                                                 function<void(LeafletRadMapWindow *)> on_create,
                                                 const bool forceNoWarning )
{
  InterSpec *viewer = InterSpec::instance();
  
  const bool showWarning = InterSpecUser::preferenceValue<bool>( "ShowMapDataWarning", viewer );
  
  if( forceNoWarning || !showWarning )
  {
    showMapWindow( meas, sample_numbers, detector_names, on_create );
    return nullptr;
  }//if( !showWarning )
  
  // Show a warning dialog about requesting the data, before proceeding
  const char *title = "Before Proceeding";
  string msg =
  "Map tiles will be requested from <a href=\"https://arcgis.com\">https://arcgis.com</a>."
  " No radiation data will leave your device, but requests for"
  " map tiles encompassing the measurements will be made to this service.";
  
  const string user_key = LeafletRadMap::get_user_arcgis_key();
  if( user_key.length() > 6 )
  {
    msg += "<p>Your custom arcgis key starting with '" + user_key.substr(0,6) + "' will be used"
    " to request map tiles.</p>";
  }else if( !user_key.empty() )
  {
    msg += "<p>An invalid arcgis key was specified in arcgis_key.txt, so the built-in key will"
    " be used.</p>";
  }
  
  SimpleDialog *dialog = new SimpleDialog( title );
  dialog->setWidth( 400 );
  
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
  accept->clicked().connect( boost::bind( &showMapWindow, meas, sample_numbers, detector_names, on_create) );
  WPushButton *cancel = dialog->addButton( "Cancel" );
  cancel->clicked().connect( std::bind([](){
    InterSpec *viewer = InterSpec::instance();
    if( viewer )
      InterSpecUser::setPreferenceValue(viewer->m_user, "ShowMapDataWarning", true, viewer);
  }) );
  accept->setFocus();
  
  return dialog;
}//


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
  
  // Include Leaflet.Measure 20e0e7d (20230426), see https://github.com/ptma/Leaflet.Measure/commits/master/
  wApp->useStyleSheet( "InterSpec_resources/assets/js/leaflet/Leaflet.Measure/src/leaflet.measure.css" );
  wApp->require( "InterSpec_resources/assets/js/leaflet/Leaflet.Measure/src/leaflet.measure.js" );
  
  // Finally include our JS/CSS
  wApp->useStyleSheet( "InterSpec_resources/LeafletRadMap.css" );
  wApp->require( "InterSpec_resources/LeafletRadMap.js" );

  
  m_displaySamples.connect( boost::bind( &LeafletRadMap::handleLoadSamples, this,
                                        boost::placeholders::_1, boost::placeholders::_2 ) );
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( viewer )
    viewer->displayedSpectrumChanged().connect( this, &LeafletRadMap::handleDisplayedSpectrumChanged );
}//LeafletRadMap


LeafletRadMap::~LeafletRadMap()
{
}//~LeafletRadMap()
  

#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
std::string LeafletRadMap::get_user_arcgis_key()
{
  try
  {
    // A user specified arcgis 
    const string user_data_dir = InterSpec::writableDataDirectory();
    const string user_key_file = SpecUtils::append_path(user_data_dir, "arcgis_key.txt" );
    
    string user_key;
    
    
    if( !SpecUtils::is_file( user_key_file ) )
    {
#if( BUILD_AS_ELECTRON_APP || BUILD_AS_OSX_APP || BUILD_AS_WX_WIDGETS_APP )
      try
      {
        const string app_data_dir = InterSpec::staticDataDirectory();
        
        const InterSpecServer::DesktopAppConfig app_config
                      = InterSpecServer::DesktopAppConfig::init( app_data_dir, user_data_dir );
        
        if( app_config.m_arcgis_key.empty() )
          return "";
        
        user_key = app_config.m_arcgis_key;
      }catch( std::exception & )
      {
        // This error would have already been reported to the user on application startup.
        return "";
      }
#else
      return "";
#endif
    }else
    {
      std::vector<char> data;
      SpecUtils::load_file_data( user_key_file.c_str(), data );
      user_key = string( begin(data), end(data) );
    }//if( "arcgis_key.txt" not present ) / else ( "arcgis_key.txt" is present )
    
    // Remove characters we dont expect
    while( !user_key.empty() )
    {
      const size_t pos = user_key.find_first_of( "\t\n \r\n{}'\"" );
      if( pos == std::string::npos )
        break;
      user_key.erase( begin(user_key) + pos );
    }//while( !user_key.empty() )
      
    SpecUtils::trim(user_key); //jic other white-space characters?
      
    // And to be super-careful
    auto wkey = WString::fromUTF8(user_key);
    Wt::Utils::removeScript(wkey);
    user_key = wkey.toUTF8();
      
    //We could test the key starts with "AAPK", but we wont right now
    //  Also, looks like it should just contain letters, numbers, _, and -  
    
    return user_key;
  }catch( std::exception &e )
  {
    
  }//
  
  return "empty";
}//std::string LeafletRadMap::get_user_arcgis_key()
#endif


void LeafletRadMap::defineJavaScript()
{
  string key = LEAFLET_MAPS_KEY;
  
  // Look for a user specified maps key.
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  key = LeafletRadMap::get_user_arcgis_key();
  
  if( key.length() > 6 )
  {
    passMessage( "Will use user-specified arcgis key, begining with '"
                + key.substr(0,6)
                + "' to request maps with.",
                WarningWidget::WarningMsgLevel::WarningMsgInfo );
  }else if( key.length() )
  {
    passMessage( "There was a arcgis_key.txt file, but it did not contain a valid key",
                WarningWidget::WarningMsgLevel::WarningMsgHigh );
  }
#endif
  
  if( key.length() < 6 )
    key = LEAFLET_MAPS_KEY;
  
  string options = "{apiKey: '" + key + "'}";
  
  setJavaScriptMember( "map", "new LeafletRadMap(" + jsRef() + "," + options + ");");
  
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



std::string LeafletRadMap::createGeoLocationJson( const std::shared_ptr<const SpecMeas> &meas,
                                                 const std::set<int> &sample_to_include,
                                                 const std::vector<std::string> &det_to_include,
                                                 const std::set<int> &foreground_samples,
                                                 const std::set<int> &background_samples,
                                                 const std::set<int> &secondary_samples )
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
    if( !sample_to_include.count(sample) && !sample_to_include.empty() )
      continue;
    
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
      const string &det_name = m->detector_name();
      if( !det_to_include.empty()
         && std::find( begin(det_to_include), end(det_to_include), det_name ) == end(det_to_include) )
        continue;
        
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
      const auto duration = meas_start_time - start_times_offset;
      const auto millisecs = chrono::duration_cast<chrono::milliseconds>(duration);
      sample_json["timeOffset"] = static_cast<long long int>( millisecs.count() );
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
    
    int displ_type = 3;
    if( foreground_samples.count(sample) )
      displ_type = 0;
    else if( secondary_samples.count(sample) )
      displ_type = 1;
    else if( background_samples.count(sample) )
      displ_type = 2;
    sample_json["disp"] = displ_type;
    
    
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
    json["startTime"] = static_cast<long long int>( millisecs.count() );
  }//if( !SpecUtils::is_special(start_times_offset) )
  
  json["samples"] = samplesData;
  
  return Wt::Json::serialize(json);
}//void createGeoLocationJson(...)


void LeafletRadMap::handleDisplayedSpectrumChanged()
{
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  
  if( !m_meas || !viewer )
    return;
  
  const shared_ptr<const SpecMeas> fore_meas = viewer->measurment( SpecUtils::SpectrumType::Foreground );
  set<int> foreground_samples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  if( fore_meas != m_meas )
    foreground_samples.clear();
  
  const shared_ptr<const SpecMeas> back_meas = viewer->measurment( SpecUtils::SpectrumType::Background );
  set<int> background_samples = viewer->displayedSamples(SpecUtils::SpectrumType::Background);
  if( back_meas != m_meas )
    background_samples.clear();
  
  const shared_ptr<const SpecMeas> second_meas = viewer->measurment( SpecUtils::SpectrumType::SecondForeground );
  set<int> secondary_samples = viewer->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
  if( second_meas != m_meas )
    secondary_samples.clear();
  
  vector<string> det_to_include = viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  if( fore_meas != m_meas )
    det_to_include = m_meas->detector_names();
  
  const string json = createGeoLocationJson( m_meas, m_meas->sample_numbers(),
                                             det_to_include, foreground_samples,
                                             background_samples, secondary_samples );
  
  doJavaScript( m_jsmap +  ".setData( " + json + ", true );" );
}//void handleDisplayedSpectrumChanged();


void LeafletRadMap::handleLoadSamples( const std::string &samples, const std::string &meas_type )
{
  try
  {
    SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
    
    if( meas_type == "foreground" )
      type = SpecUtils::SpectrumType::Foreground;
    else if( meas_type == "background" )
      type = SpecUtils::SpectrumType::Background;
    else if( meas_type == "secondary" )
      type = SpecUtils::SpectrumType::SecondForeground;
    else
      throw runtime_error( "couldnt decode spectrum type '" + meas_type + "'." );
    
    vector<int> sample_numbers;
    SpecUtils::split_to_ints( samples.c_str(), samples.length(), sample_numbers );
    
    if( sample_numbers.empty() )
    {
      string samples_str = samples;
      if( samples_str.length() > 12 )
        samples_str = samples_str.substr(0,12) + "...";
      throw runtime_error( "couldnt decode sample numbers, invalid format? ('" + samples_str + "')" );
    }//if( sample_numbers.empty() )
    
    const set<int> samplenums( begin(sample_numbers), end(sample_numbers) );
    
    InterSpec *viewer = InterSpec::instance();
    assert( viewer );
    if( !viewer )
      throw logic_error( "No valid InterSpec instance" );
    
    const std::shared_ptr<SpecMeas> prev_meas = viewer->measurment( type );
    
    if( prev_meas && (prev_meas == m_meas) )
    {
      viewer->changeDisplayedSampleNums( samplenums, type );
    }else
    {
      const std::shared_ptr<SpecMeas> fore_meas = viewer->measurment( SpecUtils::SpectrumType::Foreground );
      const std::shared_ptr<SpecMeas> back_meas = viewer->measurment( SpecUtils::SpectrumType::Background );
      const std::shared_ptr<SpecMeas> second_meas = viewer->measurment( SpecUtils::SpectrumType::SecondForeground );
      
      std::shared_ptr<SpecMeas> non_const_meas;
      if( m_meas == fore_meas )
        non_const_meas = fore_meas;
      else if( m_meas == back_meas )
        non_const_meas = back_meas;
      else if( m_meas == second_meas )
        non_const_meas = second_meas;
      
      if( !non_const_meas )
        throw runtime_error( "Spectrum file map was loaded from is not currently displayed." );
      
      viewer->setSpectrum( non_const_meas, samplenums, type );
    }
  }catch( std::exception &e )
  {
    passMessage( "Error loading map-selected samples: " + string(e.what()),
                WarningWidget::WarningMsgLevel:: WarningMsgHigh );
  }//try / catch
}//void handleLoadSamples( const std::vector<int> &samples, std::string meas_type );


void LeafletRadMap::displayMeasurementOnMap( const std::shared_ptr<const SpecMeas> &meas,
                                            std::set<int> sample_numbers,
                                            std::vector<std::string> detector_names )
{
  if( sample_numbers.empty() && meas )
    sample_numbers = meas->sample_numbers();
  if( detector_names.empty() && meas )
    detector_names = meas->detector_names();
  
  InterSpec *viewer = InterSpec::instance();
  assert( viewer );
  if( !viewer )
    return;
  
  m_meas = meas;
  
  const auto foremeas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
  const auto backmeas = viewer->measurment(SpecUtils::SpectrumType::Background);
  const auto secondmeas = viewer->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  set<int> fore_samples, back_samples, secon_samples;
  if( meas == foremeas )
    fore_samples = viewer->displayedSamples(SpecUtils::SpectrumType::Foreground);
  if( meas == backmeas )
    back_samples = viewer->displayedSamples(SpecUtils::SpectrumType::Background);
  if( meas == secondmeas )
    secon_samples = viewer->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
  
  const string json = createGeoLocationJson( meas, sample_numbers, detector_names,
                                            fore_samples, back_samples, secon_samples );
  
  //cout << "displayMeasurementOnMap --> " << json << endl;
  
  doJavaScript( m_jsmap +  ".setData( " + json + ", false );" );
}//void displayMeasurementOnMap(...)


